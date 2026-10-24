"""BEM against the Fourier model of layout.py: diffusion coefficients, region summaries, figures, results.json.

    python checks/bem/analysis.py          (after sweep.py)

Diffusion coefficient AES-4id / ISO 17497-2 over 37 receivers (5 degrees), energy averaged over the 6 frequencies of
each third-octave band (layout.band_energy, layout.diffusion), normalised per incidence to the flat reference of the
same width: d_n = (d - d_flat) / (1 - d_flat). BEM reference: rigid box W x T (wells filled). Fourier reference: its
own flat strip (layout.Deep.flat). "opt" = mean over 0, +-30, +-60 degrees, the score the optimiser maximised.
Polar plots: dB re the maximum of the flat reference of the same model at that frequency and incidence.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import geometry as geo  # noqa: E402
import style  # noqa: E402
import sweep  # noqa: E402
from style import plt  # noqa: E402

L = geo.layout
BANDS = np.array(L.BANDS, float)
INC = sweep.INC
I0, I45 = INC.index(0.0), INC.index(45.0)
IOPT = [INC.index(a) for a in (0.0, 30.0, -30.0, 60.0, -60.0)]
NB = len(BANDS)


def load(name):
    z = np.load(sweep.RAW / f"{name}.npz")
    return dict(cfg=json.loads(str(z["cfg"])), f=z["f"], F37=z["F37"], FP=z["FP"], fpolar=z["fpolar"],
                elements=z["elements"])


def band_E(F37):
    """(F, R, I) far fields of the band frequencies -> (B, I, R) band energies."""
    E = np.abs(F37[:NB * sweep.PER_BAND]).transpose(0, 2, 1) ** 2
    return E.reshape(NB, sweep.PER_BAND, *E.shape[1:]).mean(1)


def dn_of(E, E0):
    d, d0 = L.diffusion(E), L.diffusion(E0)
    return (d - d0) / (1 - d0)                               # (B, I)


def fourier_fields(size, order, ncol, f, th):
    m = geo.fourier(size, f, th, INC)
    return m.wall(order), m.flat(ncol)                       # (F, I, R)


def size_info(size):
    b = geo.qrd.deep_band(size)
    W = geo.qrd.deep_dims(size)[0]
    return dict(f0=b["f0"], f2=2 * b["f0"], f7=7 * b["f0"], W=W, f_w=b["f_width"])


def regions(size, ncol):
    """Band indices of the frequency regions the question is about."""
    s = size_info(size)
    fw = geo.qrd.C / (ncol * s["W"] / 1000)                  # c / wall width
    lo, hi = 2 ** (-1 / 6), 2 ** (1 / 6)
    in7 = lambda fc: fc * lo <= s["f7"] < fc * hi
    by = {
        "below c/W": lambda fc: fc < min(fw, s["f0"]),
        "c/W .. f0": lambda fc: fw <= fc < s["f0"],
        "f0 .. 2 f0": lambda fc: s["f0"] <= fc < s["f2"],            # a single tile: this is still below c/W
        "2 f0 .. 7 f0": lambda fc: s["f2"] <= fc < s["f7"] and not in7(fc),
        "7 f0 band": in7,
        "above 7 f0": lambda fc: fc * lo > s["f7"],
    }
    return {r: [i for i, fc in enumerate(BANDS) if test(fc)] for r, test in by.items()}, fw


def main():
    cfgs = sweep.configs()
    res = dict(method=__doc__.strip(), bands_hz=BANDS.tolist(), incidences=INC, sizes={})
    for size in sweep.SIZES:
        walls = sweep.walls_of(size)
        s = size_info(size)
        fs_all = sweep.frequencies(size)
        fb = fs_all[:NB * sweep.PER_BAND]
        out = dict(info=s, walls=walls, cases={})
        for n in [1] + sorted(walls):
            flat = load(f"{size}-flat-{n}")
            E0 = band_E(flat["F37"])
            kinds = [("tile", [1])] if n == 1 else [("layout", walls[n]), ("periodic", [1] * n)]
            case = dict(ncol=n, width_mm=n * s["W"], c_over_W=geo.qrd.C / (n * s["W"] / 1000),
                        d_flat_bem=L.diffusion(E0).round(4).tolist())
            reg, fw = regions(size, n)
            case["regions"] = reg
            for kind, order in kinds:
                name = f"{size}-tile-1" if kind == "tile" else f"{size}-{'wall' if kind == 'layout' else kind}-{n}"
                bem = load(name)
                assert np.allclose(bem["f"], fs_all)
                dnb = dn_of(band_E(bem["F37"]), E0)
                Fw, Ff = fourier_fields(size, order, n, fb, sweep.TH37)
                Ew = (np.abs(Fw) ** 2).reshape(NB, sweep.PER_BAND, *Fw.shape[1:]).mean(1)
                Ef = (np.abs(Ff) ** 2).reshape(NB, sweep.PER_BAND, *Ff.shape[1:]).mean(1)
                dnf = dn_of(Ew, Ef)
                c = dict(order=order, elements=int(bem["elements"].max()))
                for lab, dn in (("bem", dnb), ("fourier", dnf)):
                    c[f"dn_{lab}_0"] = dn[:, I0].round(3).tolist()
                    c[f"dn_{lab}_45"] = dn[:, I45].round(3).tolist()
                    c[f"dn_{lab}_m45"] = dn[:, INC.index(-45.0)].round(3).tolist()
                    c[f"dn_{lab}_opt"] = dn[:, IOPT].mean(1).round(3).tolist()
                # bands where the BEM flat box is already diffuse (d_flat > 0.9) are left out of the means: d_n
                # divides by 1 - d_flat ~ 0 there and swings by whole units on tiny differences
                d0 = np.array(case["d_flat_bem"])
                ok = {"0": d0[:, I0] <= 0.9, "45": d0[:, I45] <= 0.9, "m45": d0[:, INC.index(-45.0)] <= 0.9,
                      "opt": d0[:, IOPT].mean(1) <= 0.9}
                c["region_means"] = {}
                for r, idx in reg.items():
                    row = {}
                    for k in [k for k in c if k.startswith("dn_")]:
                        sel = [i for i in idx if ok[k.rsplit("_", 1)[1]][i]]
                        row[k] = float(np.mean(np.array(c[k])[sel])) if sel else None
                    c["region_means"][r] = row
                # lowest band where d_n >= 0.3 in two consecutive bands
                # (bands with a diffuse flat box do not count)
                c["first_band_dn_0.3"] = {}
                for k in [k for k in c if k.startswith("dn_")]:
                    v = (np.array(c[k]) >= 0.3) & ok[k.rsplit("_", 1)[1]]
                    hit = [i for i in range(NB - 1) if v[i] and v[i + 1]]
                    c["first_band_dn_0.3"][k] = float(BANDS[hit[0]]) if hit else None
                case[kind] = c
            out["cases"][str(n)] = case
        res["sizes"][size] = out
        figures(size, out)
    (OUT / "results.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


# ---------------------------------------------------------------- figures
def mark_freqs(ax, s, fw=None):
    marks = [(s["f0"], "f0"), (s["f2"], "2 f0"), (s["f7"], "7 f0")]
    if fw:
        if abs(fw / s["f2"] - 1) < 0.05:
            marks[1] = (s["f2"], "2 f0 = c/W")
        elif abs(fw / s["f0"] - 1) < 0.12:
            marks[0] = (s["f0"], "f0, c/W")
        else:
            marks.append((fw, "c/W"))
    for f, lab in marks:
        ax.axvline(f, color=style.FLAT, lw=0.9, ls=(0, (1, 2)), zorder=0)
        ax.text(f, 1.02, lab, rotation=0, ha="center", va="bottom", fontsize=10, color=style.INK2,
                transform=ax.get_xaxis_transform())


def shade_diffuse_flat(ax, d_flat):
    """Bands where the flat box itself already scatters almost evenly (d_flat > 0.9): d_n divides by ~0 there."""
    lo, hi = 2 ** (-1 / 6), 2 ** (1 / 6)
    first = True
    for fc, d0 in zip(BANDS, d_flat):
        if d0 > 0.9:
            ax.axvspan(fc * lo, fc * hi, color=style.BAND, zorder=0, lw=0,
                       label="flat box already diffuse (d > 0.9)" if first else None)
            first = False


def dn_axes(ax):
    ax.set_xscale("log")
    ax.set_xlim(180, 11200)
    ticks = [200, 500, 1000, 2000, 5000, 10000]
    ax.set_xticks(ticks); ax.set_xticklabels(["200", "500", "1k", "2k", "5k", "10k"])
    ax.minorticks_off()
    ax.set_ylim(-0.25, 1.0)
    ax.axhline(0, color=style.INK2, lw=0.8)
    ax.grid(True, axis="y")
    ax.set_xlabel("third-octave band (Hz)")


def figures(size, out):
    s = out["info"]
    label = geo.qrd.BED_SHORT[size]
    fp = sweep.polar_freqs(size)
    fplab = ["200 Hz", "500 Hz", f"f0 {fp[2]:.0f} Hz", f"2 f0 {fp[3]:.0f} Hz", f"7 f0 {fp[4]:.0f} Hz"]
    # summary: d_n averaged over the incidences of the optimiser, layout vs one tile repeated, per wall width
    ns = [n for n in sorted(map(int, out["cases"])) if n > 1]
    fig, axs = plt.subplots(1, max(1, len(ns)), figsize=(6.2 * max(1, len(ns)), 5.2), sharey=True)
    for ax, n in zip(np.atleast_1d(axs), ns):
        case = out["cases"][str(n)]
        dn_axes(ax)
        shade_diffuse_flat(ax, np.array(case["d_flat_bem"])[:, IOPT].mean(1))
        ax.plot(BANDS, case["layout"]["dn_bem_opt"], color=style.BLUE, label="BEM, wall order")
        ax.plot(BANDS, case["layout"]["dn_fourier_opt"], color=style.BLUE, ls=(0, (4, 2)), lw=1.6,
                label="Fourier, wall order")
        ax.plot(BANDS, case["periodic"]["dn_bem_opt"], color=style.ORANGE, label="BEM, one tile repeated")
        ax.plot(BANDS, case["periodic"]["dn_fourier_opt"], color=style.ORANGE, ls=(0, (4, 2)), lw=1.6,
                label="Fourier, one tile repeated")
        mark_freqs(ax, s, case["c_over_W"])
        ax.set_title(f"{n} tiles, {n * s['W'] / 10:.0f} cm", pad=22)
    np.atleast_1d(axs)[0].set_ylabel("normalised diffusion d_n")
    np.atleast_1d(axs)[0].legend(loc="upper left", fontsize=11)
    fig.suptitle(f"Deep {label}: diffusion against a flat box of the same width, mean of incidences 0, ±30, ±60°",
                 y=1.02, fontsize=15)
    fig.savefig(OUT / f"fig-{size}-summary.png", dpi=110, bbox_inches="tight")
    plt.close(fig)

    for n in sorted(map(int, out["cases"])):
        case = out["cases"][str(n)]
        tile = n == 1
        flat = load(f"{size}-flat-{n}")
        kinds = ["tile"] if tile else ["layout", "periodic"]
        bems = {k: load(f"{size}-tile-1" if k == "tile" else f"{size}-{'wall' if k == 'layout' else k}-{n}")
                for k in kinds}
        fig = plt.figure(figsize=(22, 13))
        gs = fig.add_gridspec(3, 5, height_ratios=[1.15, 1.0, 1.0], hspace=0.32, wspace=0.25)
        for j, (ii, lab) in enumerate(((I0, "0°"), (I45, "45°"))):
            ax = fig.add_subplot(gs[0, 2 * j: 2 * j + 2])
            dn_axes(ax)
            shade_diffuse_flat(ax, np.array(case["d_flat_bem"])[:, ii])
            key = "0" if ii == I0 else "45"
            main = "tile" if tile else "layout"
            ax.plot(BANDS, case[main][f"dn_bem_{key}"], color=style.BLUE,
                    label="BEM, one tile" if tile else "BEM, wall order")
            ax.plot(BANDS, case[main][f"dn_fourier_{key}"], color=style.BLUE, ls=(0, (4, 2)), lw=1.6,
                    label="Fourier, one tile" if tile else "Fourier, wall order")
            if not tile:
                ax.plot(BANDS, case["periodic"][f"dn_bem_{key}"], color=style.ORANGE, label="BEM, one tile repeated")
                ax.plot(BANDS, case["periodic"][f"dn_fourier_{key}"], color=style.ORANGE, ls=(0, (4, 2)), lw=1.6,
                        label="Fourier, one tile repeated")
            mark_freqs(ax, s, case["c_over_W"])
            ax.set_title(f"Diffusion d_n, incidence {lab}", pad=22)
            if j == 0:
                ax.set_ylabel("normalised diffusion d_n")
                ax.legend(loc="upper left", fontsize=11, ncol=1)
        ax = fig.add_subplot(gs[0, 4]); ax.axis("off")
        order = "tile 1" if tile else " ".join(map(str, case["layout"]["order"]))
        ax.text(0, 0.95, f"Deep {label}\n{n} x {s['W']:.0f} mm = {n * s['W']:.0f} mm\n\n"
                f"model order (STL X): {order}\n\nf0 = {s['f0']:.0f} Hz\nc / period = 2 f0 = {s['f2']:.0f} Hz\n"
                f"c / wall width = {case['c_over_W']:.0f} Hz\n7 f0 (flat plate) = {s['f7']:.0f} Hz\n\n"
                "Polar: dB re the flat-plate peak\nof the same model.\nGrey: flat box, BEM.",
                va="top", fontsize=13, color=style.INK)
        for r, (ii, lab) in enumerate(((I0, "0°"), (I45, "45°"))):
            th = sweep.THP
            Fl = flat["FP"][:, :, ii]                           # (P, T)
            Fw, Ff = fourier_fields(size, case[kinds[0]]["order"], n, fp, th)
            for c in range(5):
                ax = style.polar_axes(fig, gs[r + 1, c])
                refb = np.abs(Fl[c]).max()
                reff = np.abs(Ff[c, ii]).max()
                lv = lambda F, ref: np.maximum(20 * np.log10(np.abs(F) / ref + 1e-12), -40)
                ax.plot(np.radians(th), lv(Fl[c], refb), color=style.FLAT, lw=1.3, label="flat box, BEM")
                if not tile:
                    ax.plot(np.radians(th), lv(bems["periodic"]["FP"][c, :, ii], refb), color=style.ORANGE, lw=1.5,
                            label="BEM, one tile repeated")
                ax.plot(np.radians(th), lv(Fw[c, ii], reff), color=style.BLUE, lw=1.5, ls=(0, (4, 2)),
                        label="Fourier, " + ("one tile" if tile else "wall order"))
                ax.plot(np.radians(th), lv(bems[kinds[0]]["FP"][c, :, ii], refb), color=style.BLUE, lw=1.8,
                        label="BEM, " + ("one tile" if tile else "wall order"))
                ax.set_title(f"{fplab[c]}, incidence {lab}", pad=10, fontsize=13)
                if r == 1 and c == 0:
                    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.05), fontsize=11, ncol=2)
        fig.savefig(OUT / f"fig-{size}-{'tile' if tile else f'wall-{n}'}.png", dpi=100, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
