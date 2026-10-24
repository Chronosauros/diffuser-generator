"""Every order of a 3-tile P1S / X1 deep wall in BEM and in the Fourier model, ranked by the optimiser's score.

    python checks/bem/order_ranking.py [--workers 8]

Score as in layout.deep_layouts: mean d_n over the incidences 0, +-30, +-60 degrees and the third-octave bands from f0
to below N f0 (800 Hz - 4 kHz on the P1S / X1), 3 frequencies per band, 37 receivers, d_n against the flat reference of
the same width (BEM: rigid box W x T; Fourier: its own flat strip). There are 7 * 6 * 5 = 210 orders. Turning a tile
upside down gives tile 9 - t (tile 1 stays), so the order a b c seen in a mirror is m(c) m(b) m(a): the same wall
mirrored, with the same score over symmetric incidences and receivers. BEM runs one order of every mirror pair (and both
orders of two pairs, to check that); the Fourier model runs all 210.
Band energies per order go to out/ranking-P1S-3-bem.jsonl as they come (a stopped run resumes), the summary to
out/ranking-P1S-3.json and the figure to out/fig-P1S-3-ranking.png.
"""
from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import bem2d  # noqa: E402
import geometry as geo  # noqa: E402

L = geo.layout
SIZE, NCOL = "P1S", 3
PER_BAND = 3
INC = [0.0, 30.0, -30.0, 60.0, -60.0]
TH = np.arange(-90, 90.01, 5.0)
JSONL = OUT / "ranking-P1S-3-bem.jsonl"
PUBLISHED = tuple(json.loads((geo.REF / "out" / "layouts.json").read_text(encoding="utf-8"))["deep"]["L"]["walls"]["3"])
DEADLINE_S = 3 * 3600


def bands():
    b = geo.qrd.deep_band(SIZE)
    sel = [i for i, fc in enumerate(L.BANDS) if b["f0"] / 2 ** (1 / 6) <= fc <= b["f_flat"] / 2 ** (1 / 6)]
    return [L.BANDS[i] for i in sel]


def freqs():
    return L.band_freqs(bands(), PER_BAND)


def mirror(order):
    m = {1: 1, 2: 7, 3: 6, 4: 5, 5: 4, 6: 3, 7: 2}
    return tuple(m[t] for t in order[::-1])


def canonical(order):
    return min(tuple(order), mirror(order))


def band_energy(F):
    """(F, R, I) far fields -> (B, I, R) band energies."""
    E = np.abs(F).transpose(0, 2, 1) ** 2
    return E.reshape(-1, PER_BAND, *E.shape[1:]).mean(1)


def score(E, E0):
    d, d0 = L.diffusion(E), L.diffusion(E0)
    return (d - d0) / (1 - d0)                       # (B, I)


def task(order):
    poly = geo.centre_x(geo.flat(SIZE, NCOL) if order == "flat" else geo.wall(SIZE, list(order)))
    rows = []
    for f in freqs():
        k = 2 * np.pi * f / bem2d.C
        ms = bem2d.mesh(poly, bem2d.element_size(f))
        rows.append(bem2d.farfield(ms, k, bem2d.solve(ms, k, INC), TH))
    return order, band_energy(np.stack(rows)).astype(np.float64)


def run_bem(workers):
    done = {}
    if JSONL.exists():
        for line in JSONL.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            done[r["order"] if r["order"] == "flat" else tuple(r["order"])] = np.array(r["E"])
    reps = sorted({canonical(o) for o in itertools.permutations(range(1, 8), NCOL)})
    extra = [mirror(reps[0]), mirror(reps[-1])]                 # both sides of two pairs: the mirror check
    if PUBLISHED not in reps:
        extra.append(PUBLISHED)
    todo = [o for o in ["flat"] + reps + extra if o not in done]
    print(f"BEM: {len(todo)} runs to do ({len(done)} done), {len(freqs())} frequencies each, {workers} workers", flush=True)
    t0 = time.time()
    pool = get_context("spawn").Pool(workers)
    try:
        with JSONL.open("a", encoding="utf-8") as fh:
            for i, (order, E) in enumerate(pool.imap_unordered(task, todo)):
                done[order] = E
                fh.write(json.dumps(dict(order=order if order == "flat" else list(order), E=E.round(9).tolist())) + "\n")
                fh.flush()
                if (i + 1) % 10 == 0 or i + 1 == len(todo):
                    print(f"{time.time() - t0:7.0f}s  {i + 1}/{len(todo)}", flush=True)
                if time.time() - t0 > DEADLINE_S:
                    raise TimeoutError("ranking deadline")
    finally:
        pool.terminate()
        pool.join()
    return done


def run_fourier():
    m = geo.fourier(SIZE, freqs(), TH, INC)
    E0 = band_energy(m.flat(NCOL).transpose(0, 2, 1))
    out = {}
    for order in itertools.permutations(range(1, 8), NCOL):
        out[order] = score(band_energy(m.wall(list(order)).transpose(0, 2, 1)), E0)
    return out


def ranks(values):
    """Rank 1 = best (highest), ties share the mean rank."""
    v = np.asarray(values)
    order = np.argsort(-v, kind="stable")
    r = np.empty(len(v))
    r[order] = np.arange(1, len(v) + 1)
    for x in np.unique(v):
        idx = v == x
        r[idx] = r[idx].mean()
    return r


def kendall(a, b):
    a, b = np.asarray(a), np.asarray(b)
    s = np.sign(a[:, None] - a[None]) * np.sign(b[:, None] - b[None])
    n = len(a)
    return float(s[np.triu_indices(n, 1)].sum() / (n * (n - 1) / 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    bem = run_bem(a.workers)
    four = run_fourier()
    bnd = bands()
    b = geo.qrd.deep_band(SIZE)
    low = [i for i, fc in enumerate(bnd) if fc < b["f_period"]]  # f0 .. 2 f0: 800, 1000, 1250 Hz
    E0 = bem["flat"]
    mirror_check = []
    for o in [o for o in bem if o != "flat" and mirror(o) in bem and o < mirror(o)]:
        s1, s2 = score(bem[o], E0).mean(), score(bem[mirror(o)], E0).mean()
        mirror_check.append(dict(order=list(o), mirror=list(mirror(o)), score=round(float(s1), 5), score_mirror=round(float(s2), 5)))
    orders = list(four)
    sb = np.array([score(bem[canonical(o)], E0).mean() for o in orders])
    sf = np.array([four[o].mean() for o in orders])
    sb_low = np.array([score(bem[canonical(o)], E0)[low].mean() for o in orders])
    sf_low = np.array([four[o][low].mean() for o in orders])
    # one entry per mirror pair for the statistics (the 210 orders hold every pair twice)
    pair = sorted({canonical(o) for o in orders})
    ip = [orders.index(o) for o in pair]
    rb, rf = ranks(sb[ip]), ranks(sf[ip])
    rb_low, rf_low = ranks(sb_low[ip]), ranks(sf_low[ip])
    spearman = float(np.corrcoef(rb, rf)[0, 1])
    spearman_low = float(np.corrcoef(rb_low, rf_low)[0, 1])
    pub = pair.index(canonical(PUBLISHED))
    top = lambda r, k: {pair[i] for i in np.argsort(r)[:k]}
    res = dict(
        method=__doc__.strip(), size=SIZE, columns=NCOL, bands_hz=bnd, low_bands_hz=[bnd[i] for i in low], incidences=INC,
        published=list(PUBLISHED), pairs=len(pair), mirror_check=mirror_check,
        spearman_rank=spearman, kendall_tau=kendall(sb[ip], sf[ip]),
        spearman_rank_f0_2f0=spearman_low, kendall_tau_f0_2f0=kendall(sb_low[ip], sf_low[ip]),
        published_rank_bem=float(rb[pub]), published_rank_fourier=float(rf[pub]),
        published_rank_bem_f0_2f0=float(rb_low[pub]), published_rank_fourier_f0_2f0=float(rf_low[pub]),
        top10_overlap=len(top(rb, 10) & top(rf, 10)),
        pairs_table=[dict(order=list(pair[j]), mirror=list(mirror(pair[j])), bem=round(float(sb[ip][j]), 4),
                          fourier=round(float(sf[ip][j]), 4), bem_f0_2f0=round(float(sb_low[ip][j]), 4),
                          fourier_f0_2f0=round(float(sf_low[ip][j]), 4), rank_bem=float(rb[j]), rank_fourier=float(rf[j]))
                     for j in np.argsort(rb)],
    )
    (OUT / "ranking-P1S-3.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    figure(res)
    print(json.dumps({k: res[k] for k in ("pairs", "spearman_rank", "kendall_tau", "spearman_rank_f0_2f0",
                                         "published_rank_bem", "published_rank_fourier", "top10_overlap")}))


def figure(res):
    import style
    from style import plt
    rows = res["pairs_table"]
    pub = res["published"]
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 6))
    for ax, kb, kf, title in ((axs[0], "bem", "fourier", "800 Hz - 4 kHz (the optimiser's score)"),
                              (axs[1], "bem_f0_2f0", "fourier_f0_2f0", "f0 - 2 f0 only (800 - 1250 Hz)")):
        x = np.array([r[kf] for r in rows]); y = np.array([r[kb] for r in rows])
        ax.scatter(x, y, s=22, color=style.BLUE, alpha=0.55, lw=0, label="one order and its mirror image")
        best = int(np.argmax(y))
        ip = [i for i, r in enumerate(rows) if pub in (r["order"], r["mirror"])][0]
        ax.scatter([x[ip]], [y[ip]], s=90, color=style.ORANGE, zorder=3, label="published order")
        ax.scatter([x[best]], [y[best]], s=90, facecolor="none", edgecolor=style.INK, lw=1.6, zorder=3,
                   label="best in BEM for this score")
        ax.annotate(" ".join(map(str, pub)), (x[ip], y[ip]), xytext=(0, -16), textcoords="offset points", ha="center",
                    fontsize=11, color=style.INK)
        bo = rows[best]
        ax.annotate(f"{' '.join(map(str, bo['order']))} / {' '.join(map(str, bo['mirror']))}", (x[best], y[best]),
                    xytext=(0, 10), textcoords="offset points", ha="center", fontsize=11, color=style.INK)
        lo, hi = min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02
        ax.plot([lo, hi], [lo, hi], color=style.FLAT, lw=1, ls=(0, (2, 2)), zorder=0)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.grid(True)
        ax.set_xlabel("Fourier model, mean d_n"); ax.set_ylabel("BEM, mean d_n")
        ax.set_title(title, pad=10)
    axs[0].legend(loc="upper left", fontsize=10.5)
    fig.suptitle(f"P1S / X1, 3 tiles: all {res['pairs']} orders (mirror pairs merged), Spearman rank correlation "
                 f"{res['spearman_rank']:.2f} / {res['spearman_rank_f0_2f0']:.2f}", fontsize=14)
    fig.savefig(OUT / "fig-P1S-3-ranking.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
