"""Wall layouts of the Soundblooms diffuser tiles, chosen with a far-field scattering model.

    python design/layout.py          v3 set (rule v3, qrd.RULE): copy the stage 1 wall layouts by the rule of
                                     qrd.layout_class() (no new search), evaluate them with the model on the v3 sizes
                                     -> out-v3/layouts.json and out-v3/model/*.png (out-v2 is the frozen v2 archive)
    python design/layout.py --stage1 DIR
                                     stage 1: optimise every wall (search, slow) -> DIR/layouts.json and DIR/model/*.png.
                                     design/out/ is the frozen canon of the generator test: never point DIR there.

Model (T. J. Cox, P. D'Antonio, "Acoustic Absorbers and Diffusers", 3rd ed., chapter on prediction models):
Fourier / Fraunhofer far field of a plane surface with a position-dependent reflection factor R(x).
    Deep tiles (wells with fins, 1D): R = exp(-2 i k d) per well (plane wave down the well and back), fins R = 1.
        p(theta) = sum_e w_e sinc(k w_e s / 2) exp(i k x_e s) R_e,   s = sin(theta) + sin(psi)
    Block tiles (2D, no fins): each square is a flat facet at height h, Kirchhoff phase exp(i k h (cos theta + cos psi)),
        p = sum_c px py sinc(k px u / 2) sinc(k py v / 2) exp(i k (x_c u + y_c v)) exp(i k h_c (cos theta + cos psi))
Diffusion coefficient (AES-4id-2001 / ISO 17497-2): d = [(sum E)^2 - sum E^2] / [(n - 1) sum E^2] over receivers with
equal weight (1D: semicircle every 5 degrees, 37 receivers; 2D: 400 equal-area points on the hemisphere), E averaged
over the frequencies of a third-octave band. Normalised against a flat plate of the same size:
d_n = (d - d_flat) / (1 - d_flat), per incidence, then averaged over incidences 0, +-30, +-60 degrees.
It is a prediction model for comparing layouts with each other, not a measurement: no number from it goes into a
public text.

Orders are in model (STL X) coordinates; seen from the room the STL X axis points left, so "room_left_to_right"
in layouts.json is the reversed order. Walls grow without moving tiles that are already stuck on: deep S 4 -> 7
columns, deep L 3 -> 5 -> 7 columns, blocks S 2x2 -> 2x3 -> 3x3, blocks L 2x2 -> 2x3. The optimiser maximises the
mean d_n of all wall sizes of a series; the deep score leaves out the third-octave band around N*f0.

Since 25.10.2026 the deep orders in out/layouts.json are not this optimiser's: they were replaced by the best nested
orders of the 2D BEM ranking (S 1 2 7 3 6 5 4, L 5 6 2 7 1 4 3; checks/bem/order_nested.py, docs/bem-check.md), with
the same score and the same growth rule. deep_layouts() and the stage 1 orders (S 6 3 5 4 2 7 1, L 6 3 4 2 5 7 1) stay
here as history; the block shifts are still the stage 1 ones.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qrd  # noqa: E402

STAGE1 = Path(__file__).resolve().parent / "out"          # frozen stage 1 canon (layouts.json), read only
OUT = Path(__file__).resolve().parent / "out-v3"
MODEL = OUT / "model"
C = qrd.C
BANDS = [200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000]
INC = [0.0, 30.0, -30.0, 60.0, -60.0]
RNG = np.random.default_rng(20261022)


def band_freqs(bands, per_band):
    """(n_bands * per_band,) frequencies, log-spaced inside each third-octave band."""
    off = 2 ** (np.linspace(-1 / 6, 1 / 6, per_band + 2)[1:-1]) if per_band > 1 else np.array([1.0])
    return np.concatenate([fc * off for fc in bands])


def diffusion(E):
    """AES-4id diffusion coefficient over the last axis (receivers)."""
    s, s2, n = E.sum(-1), (E ** 2).sum(-1), E.shape[-1]
    return (s ** 2 - s2) / ((n - 1) * s2)


def band_energy(p, per_band):
    E = np.abs(p) ** 2                                  # (F, I, R)
    return E.reshape(-1, per_band, *E.shape[1:]).mean(1)  # (B, I, R)


def dn(E, E_flat):
    d, d0 = diffusion(E), diffusion(E_flat)              # (B, I)
    return ((d - d0) / (1 - d0)).mean(-1)                # (B,)


# ---------------------------------------------------------------- 1D: deep tiles
def deep_elements(size, tile):
    """Strips across the tile: (x_center_m, width_m, depth_m), wells and the fins between them."""
    p = qrd.deep_params(size)
    N, pitch, fin, edge = p["N"], p["pitch"], p["fin"], p["edge"]
    d = qrd.deep_depths(size, tile)
    el = [(edge / 2, edge, 0.0), (N * pitch - edge / 2, edge, 0.0)]
    el += [(j * pitch, fin, 0.0) for j in range(1, N)]
    for n in range(N):
        xl = n * pitch + (edge if n == 0 else fin / 2)
        xr = (n + 1) * pitch - (edge if n == N - 1 else fin / 2)
        el.append(((xl + xr) / 2, xr - xl, d[n]))
    return np.array(el) / 1000.0


class Deep:
    def __init__(self, size, per_band):
        self.size, self.per_band = size, per_band
        self.W = qrd.deep_dims(size)[0] / 1000
        self.f = band_freqs(BANDS, per_band)
        self.theta = np.radians(np.arange(-90, 91, 5))
        k = 2 * np.pi * self.f / C
        s = np.sin(self.theta)[None, None, :] + np.sin(np.radians(INC))[None, :, None]   # (1, I, R)
        self.k, self.s = k[:, None, None], s
        self.T = {t: self._tile(t) for t in range(1, qrd.DEEP_TILES + 1)}

    def _tile(self, tile):
        p = 0
        for x, w, d in deep_elements(self.size, tile):
            p = p + w * np.sinc(self.k * w * self.s / (2 * np.pi)) * np.exp(1j * self.k * x * self.s) * np.exp(-2j * self.k * d)
        return p.astype(np.complex64)

    def shift(self, j):
        return np.exp(1j * self.k * (j * self.W) * self.s).astype(np.complex64)

    def flat(self, ncol):
        Wt = ncol * self.W
        return Wt * np.sinc(self.k * Wt * self.s / (2 * np.pi)) * np.exp(1j * self.k * (Wt / 2) * self.s)

    def wall(self, order):
        return sum(self.T[t] * self.shift(j) for j, t in enumerate(order))


# ---------------------------------------------------------------- 2D: block tiles
def hemisphere(n=400):
    i = np.arange(n) + 0.5
    ct = i / n                                          # equal area: cos(theta) uniform
    st = np.sqrt(1 - ct ** 2)
    ph = np.pi * (1 + 5 ** 0.5) * i
    return st * np.cos(ph), st * np.sin(ph), ct


class Blocks:
    def __init__(self, size, per_band, bands, n_rec=400):
        p = qrd.block_params(size)
        self.size, self.per_band, self.bands = size, per_band, bands
        self.Nx, self.Ny = p["Nx"], p["Ny"]
        self.W, self.H = p["W"] / 1000, p["H"] / 1000
        self.px, self.py = self.W / self.Nx, self.H / self.Ny
        self.f = band_freqs(bands, per_band)
        rx, ry, rz = hemisphere(n_rec)
        # incidences: normal, +-30 and +-60 in the horizontal plane (x), as for the deep tiles
        ix = np.sin(np.radians(INC)); iz = np.cos(np.radians(INC))
        k = 2 * np.pi * self.f / C
        self.k = k[:, None, None]
        self.u = rx[None, None, :] + ix[None, :, None]
        self.v = ry[None, None, :] + 0 * ix[None, :, None]
        self.cz = rz[None, None, :] + iz[None, :, None]
        self.cell = (self.px * self.py * np.sinc(self.k * self.px * self.u / (2 * np.pi))
                     * np.sinc(self.k * self.py * self.v / (2 * np.pi))).astype(np.complex64)
        self.ex = [np.exp(1j * self.k * ((n + 0.5) * self.px) * self.u).astype(np.complex64) for n in range(self.Nx)]
        self.ey = [np.exp(1j * self.k * ((m + 0.5) * self.py) * self.v).astype(np.complex64) for m in range(self.Ny)]
        self._T, self._at = {}, {}

    def T(self, kl):
        kl = tuple(kl)
        if kl not in self._T:
            h = np.asarray(qrd.block_heights(self.size, *kl)) / 1000
            levels, inv = np.unique(h.round(9), return_inverse=True)
            hp = [np.exp(1j * self.k * lv * self.cz).astype(np.complex64) for lv in levels]
            inv = inv.reshape(h.shape)
            p = 0
            for m in range(self.Ny):
                row = sum(self.ex[n] * hp[inv[m, n]] for n in range(self.Nx))
                p = p + row * self.ey[m]
            self._T[kl] = p * self.cell
        return self._T[kl]

    def at(self, r, c):
        if (r, c) not in self._at:
            self._at[r, c] = np.exp(1j * self.k * (c * self.W * self.u - r * self.H * self.v)).astype(np.complex64)
        return self._at[r, c]

    def flat(self, rows, cols):
        Wt, Ht = cols * self.W, rows * self.H
        return (Wt * Ht * np.sinc(self.k * Wt * self.u / (2 * np.pi)) * np.sinc(self.k * Ht * self.v / (2 * np.pi))
                * np.exp(1j * self.k * (Wt / 2 * self.u - (Ht / 2 - self.H) * self.v)))

    def wall(self, grid, shifts):
        """grid: rows of tile numbers (1-based, 0 = empty), shifts[tile-1] = (k, l)."""
        p = 0
        for r, row in enumerate(grid):
            for c, t in enumerate(row):
                if t:
                    p = p + self.T(shifts[t - 1]) * self.at(r, c)
        return p


# ---------------------------------------------------------------- optimisation
def deep_layouts(size):
    walls = {"S": [4, 7], "L": [3, 5, 7]}[size]
    band = qrd.deep_band(size)
    # design band without the third-octave band around N*f0, where every QRD reflects like a flat panel
    sel = [i for i, fc in enumerate(BANDS) if band["f0"] / 2 ** (1 / 6) <= fc <= band["f_flat"] / 2 ** (1 / 6)]
    m = Deep(size, per_band=3)
    pos = [m.shift(j) for j in range(7)]
    Tt = np.stack([m.T[t] for t in range(1, 8)])            # (7, F, I, R)
    flats = {n: band_energy(m.flat(n), m.per_band)[sel] for n in walls}
    orders = np.array(list(itertools.permutations(range(7))))
    score = np.zeros(len(orders))
    for n in walls:
        uniq, inv = np.unique(orders[:, :n], axis=0, return_inverse=True)
        sc = np.zeros(len(uniq))
        for b0 in range(0, len(uniq), 256):
            ob = uniq[b0:b0 + 256]
            p = sum(Tt[ob[:, j]] * pos[j][None] for j in range(n))         # (b, F, I, R)
            E = (np.abs(p) ** 2).reshape(len(ob), -1, m.per_band, *p.shape[2:]).mean(2)[:, sel]
            d, d0 = diffusion(E), diffusion(flats[n])[None]
            sc[b0:b0 + 256] = ((d - d0) / (1 - d0)).mean(axis=(1, 2))
        score += sc[inv.ravel()] / len(walls)
    best = orders[int(np.argmax(score))] + 1
    return [int(t) for t in best], float(score.max()), walls, sel


def block_layouts(size, iters=6000):
    p = qrd.BLOCKS[size]
    band = qrd.block_band(size)
    bands = [fc for fc in BANDS if band["f0"] / 2 ** (1 / 3) <= fc <= min(band["f_high"], 10000)]
    m = Blocks(size, per_band=2, bands=bands, n_rec=300)
    if size == "S":
        grids = {"2x2": [[1, 2], [3, 4]], "2x3": [[1, 2, 5], [3, 4, 6]], "3x3": [[1, 2, 5], [3, 4, 6], [7, 8, 9]]}
    else:
        grids = {"2x2": [[1, 2], [3, 4]], "2x3": [[1, 2, 5], [3, 4, 6]]}
    n_tiles = max(max(r) for g in grids.values() for r in g)
    flats = {g: band_energy(m.flat(len(v), len(v[0])), m.per_band) for g, v in grids.items()}
    all_shifts = [(k, l) for k in range(p["Nx"]) for l in range(p["Ny"])]

    def score(sh):
        return float(np.mean([dn(band_energy(m.wall(g, sh), m.per_band), flats[name]).mean() for name, g in grids.items()]))

    cur = [all_shifts[i] for i in RNG.choice(len(all_shifts), n_tiles, replace=False)]
    cur_s = score(cur)
    best, best_s = list(cur), cur_s
    temp0 = 0.02
    for it in range(iters):
        cand = list(cur)
        j = int(RNG.integers(n_tiles))
        free = [s for s in all_shifts if s not in cand]
        cand[j] = free[int(RNG.integers(len(free)))]
        s = score(cand)
        temp = temp0 * (1 - it / iters) + 1e-5
        if s > cur_s or RNG.random() < np.exp((s - cur_s) / temp):
            cur, cur_s = cand, s
            if s > best_s:
                best, best_s = list(cand), s
    return [list(map(int, s)) for s in best], best_s, grids


# ---------------------------------------------------------------- full evaluation and plots
def evaluate_deep(size, orders):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = Deep(size, per_band=6)
    res = {}
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for n, order in orders.items():
        Ef = band_energy(m.flat(n), m.per_band)
        per = dn(band_energy(m.wall([1] * n), m.per_band), Ef)
        ape = dn(band_energy(m.wall(order), m.per_band), Ef)
        res[str(n)] = dict(order=order, bands_hz=BANDS, dn_periodic=per.round(3).tolist(), dn_layout=ape.round(3).tolist())
        ax.plot(BANDS, per, "--", label=f"{n} columns, one tile repeated")
        ax.plot(BANDS, ape, "-", lw=2, label=f"{n} columns, layout {'-'.join(map(str, order))}")
    b = qrd.deep_band(size)
    ax.axvspan(b["f0"], b["f_high"], color="0.92", zorder=0, label="design band")
    ax.set_xscale("log"); ax.set_xticks(BANDS[::2]); ax.set_xticklabels([f"{f:g}" for f in BANDS[::2]])
    ax.set_xlabel("third-octave band (Hz)"); ax.set_ylabel("normalised diffusion (model)")
    ax.set_ylim(-0.1, 1); ax.legend(fontsize=7); ax.set_title(f"Deep {qrd.BED_SHORT.get(size, size)}: Fourier model, flat plate = 0")
    fig.savefig(MODEL / f"deep-{size}-diffusion.png", dpi=110, bbox_inches="tight")
    # polar plot at a mid-band frequency, normal incidence, widest wall
    n, order = max(orders.items())
    fc = 2000 if qrd.deep_params(size)["pitch"] < 30 else 1250
    th = np.radians(np.arange(-90, 90.5, 0.5))
    mm = Deep(size, per_band=1)
    mm.f = np.array([fc]); mm.theta = th
    mm.k = (2 * np.pi * mm.f / C)[:, None, None]; mm.s = (np.sin(th)[None, None, :] + 0 * np.sin(np.radians(INC[:1]))[None, :, None])
    mm.T = {t: mm._tile(t) for t in range(1, 8)}
    fig = plt.figure(figsize=(6, 4))
    a = fig.add_subplot(projection="polar")
    for lab, pp in (("flat plate", mm.flat(n)), ("one tile repeated", mm.wall([1] * n)), ("layout", mm.wall(order))):
        L = 20 * np.log10(np.abs(pp[0, 0]) + 1e-12); L -= L.max()
        a.plot(th, np.maximum(L, -40), label=lab)
    a.set_thetamin(-90); a.set_thetamax(90); a.set_theta_zero_location("N"); a.set_rlim(-40, 0)
    a.set_rticks([-30, -20, -10, 0]); a.set_yticklabels(["-30", "-20", "-10", "0 dB"], fontsize=7)
    a.set_title(f"Deep {qrd.BED_SHORT.get(size, size)}, {n} columns, {fc} Hz, normal incidence (dB, model)", fontsize=9)
    a.legend(fontsize=7, loc="lower left")
    fig.savefig(MODEL / f"deep-{size}-polar-{fc}.png", dpi=110, bbox_inches="tight")
    plt.close("all")
    return res


def evaluate_blocks(size, shifts, grids):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bands = [fc for fc in BANDS if fc >= 500]
    m = Blocks(size, per_band=3, bands=bands, n_rec=400)
    res = {}
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for name, g in grids.items():
        Ef = band_energy(m.flat(len(g), len(g[0])), m.per_band)
        per = dn(band_energy(m.wall(g, [(0, 0)] * 9), m.per_band), Ef)
        ape = dn(band_energy(m.wall(g, shifts), m.per_band), Ef)
        res[name] = dict(grid=g, bands_hz=bands, dn_periodic=per.round(3).tolist(), dn_layout=ape.round(3).tolist())
        ax.plot(bands, per, "--", label=f"{name}, one tile repeated")
        ax.plot(bands, ape, "-", lw=2, label=f"{name}, numbered layout")
    b = qrd.block_band(size)
    ax.axvspan(b["f0"], min(b["f_high"], 10000), color="0.92", zorder=0, label="design band")
    ax.set_xscale("log"); ax.set_xticks(bands[::2]); ax.set_xticklabels([f"{f:g}" for f in bands[::2]])
    ax.set_xlabel("third-octave band (Hz)"); ax.set_ylabel("normalised diffusion (model)")
    ax.set_ylim(-0.1, 1); ax.legend(fontsize=7); ax.set_title(f"Blocks {qrd.BED_SHORT.get(size, size)}: Fourier/Kirchhoff model, flat plate = 0")
    fig.savefig(MODEL / f"blocks-{size}-diffusion.png", dpi=110, bbox_inches="tight")
    plt.close("all")
    return res


def v2_layouts() -> dict:
    """v2 / v3 layouts: values copied from the frozen stage 1 layouts.json by the rule of qrd.layout_class(); no search.
    The orders do not depend on the wall thickness, so v3 has the same ones as v2 (the model is run again on v3 sizes)."""
    s1 = json.loads((STAGE1 / "layouts.json").read_text(encoding="utf-8"))
    out = dict(deep={}, blocks={})
    for key, p in qrd.deep_sizes().items():
        cls = qrd.layout_class("deep", p)
        src = s1["deep"][cls]
        out["deep"][key] = dict(order=src["order"], walls=src["walls"], room_left_to_right=src["room_left_to_right"],
                                from_stage1=cls)
    for key, p in qrd.BLOCKS_V2.items():
        cls = qrd.layout_class("blocks", p)
        src = s1["blocks"][cls]
        out["blocks"][key] = dict(shifts=src["shifts"], grids=src["grids"], from_stage1=cls)
    return out


def main_v2():
    """Evaluation only (model, no search) of the copied layouts on the parameters of the current rule (v3)."""
    MODEL.mkdir(parents=True, exist_ok=True)
    out = v2_layouts()
    target = OUT / "layouts.json"
    target.write_text(json.dumps(out, indent=1), encoding="utf-8")        # the layouts first, evaluation can take minutes
    for key, d in out["deep"].items():
        d["evaluation"] = evaluate_deep(key, {int(n): w for n, w in d["walls"].items()})
        print("deep", key, d["room_left_to_right"], flush=True)
    for key, d in out["blocks"].items():
        d["evaluation"] = evaluate_blocks(key, d["shifts"], d["grids"])
        print("blocks", key, d["grids"], flush=True)
    target.write_text(json.dumps(out, indent=1), encoding="utf-8")


def main():
    MODEL.mkdir(parents=True, exist_ok=True)
    out = dict(model=__doc__.split("\n\n")[2].strip(), deep={}, blocks={})
    for size in qrd.DEEP:
        order, sc, walls, _ = deep_layouts(size)
        orders = {n: order[:n] for n in walls}
        ev = evaluate_deep(size, orders)
        out["deep"][size] = dict(order=order, walls={str(n): order[:n] for n in walls},
                                 room_left_to_right={str(n): order[:n][::-1] for n in walls},
                                 score=round(sc, 3), evaluation=ev)
        print("deep", size, order, round(sc, 3), flush=True)
    for size in qrd.BLOCKS:
        shifts, sc, grids = block_layouts(size)
        ev = evaluate_blocks(size, shifts, grids)
        out["blocks"][size] = dict(shifts=shifts, grids=grids, score=round(sc, 3), evaluation=ev)
        print("blocks", size, shifts, round(sc, 3), flush=True)
    (OUT / "layouts.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    if sys.argv[1:2] == ["--stage1"]:
        if len(sys.argv) < 3:
            sys.exit("--stage1 needs a target folder")
        OUT = Path(sys.argv[2]).resolve()
        if OUT == STAGE1.resolve():
            sys.exit("design/out/ is the frozen stage 1 canon; give another folder")
        MODEL = OUT / "model"
        main()
    else:
        main_v2()
