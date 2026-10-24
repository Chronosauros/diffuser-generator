"""Checks of the BEM solver before it is used on the tiles. Writes out/validation.json and out/fig-validation-*.png.

    python checks/bem/validate.py

1. Rigid circular cylinder against the exact Bessel series (far field), incl. an interior Dirichlet resonance ka = 2.405.
2. Fictitious resonance of a closed box: the P1S flat box 231 x 134 mm has its first interior Dirichlet eigenfrequency
   at 1479.6 Hz. CBIE alone against Burton-Miller, both against Burton-Miller on a twice finer mesh.
3. Flat box against Kirchhoff (asymptotic, not exact: valid for k W >> 1 near the specular direction).
4. Standard N = 7 QRD of Cox, Angus and D'Antonio (f0 490 Hz, unit 50 mm, wells 75 mm, period 525 mm), 6 periods:
   energy of the grating lobes at m f0 (Schroeder: equal) and the flat-plate frequency 7 f0 = 3430 Hz.
5. Mesh convergence on P1S: element min(lambda/10, 2 mm), 4 Gauss points, against min(lambda/20, 1 mm), 8 points;
   the first run with lambda/8 capped at 8 mm was up to 13 % off in the diffusion coefficient (thin fins).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import h1vp, jvp

HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "out"
sys.path.insert(0, str(HERE))
import bem2d  # noqa: E402
import geometry as geo  # noqa: E402
import style  # noqa: E402
from style import plt  # noqa: E402

C = bem2d.C
OUT = {}


def db(x, ref):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-12) / ref)


def solve_ff(poly, f, psi, theta, h=None, method="bm"):
    k = 2 * np.pi * f / C
    ms = bem2d.mesh(poly, h or bem2d.element_size(f))
    p = bem2d.solve(ms, k, psi, method=method)
    return bem2d.farfield(ms, k, p, theta), len(ms)


# ---------------------------------------------------------------- 1. cylinder
def cylinder_exact(ka, a, theta, psi):
    obs = np.pi / 2 - np.radians(theta)
    prop = np.radians(-90 - psi)
    F = 0
    for n in range(int(ka + 40)):
        F = F + (1 if n == 0 else 2) * jvp(n, ka) / h1vp(n, ka) * np.cos(n * (obs - prop))
    return 4j * F


def check_cylinder():
    a = 0.1
    th = np.arange(-180, 180, 1.0)
    rows = []
    for ka in [1.0, 2.404826, 10.0, 30.0]:
        k = ka / a
        nseg = max(256, int(np.ceil(2 * np.pi * a / (2 * np.pi / k / 12))))
        ms = bem2d.mesh(geo.cylinder(a, nseg), 1.0)
        Fe = cylinder_exact(ka, a, th, 0.0)
        row = dict(ka=ka, elements=len(ms))
        for method in ("cbie", "bm"):
            F = bem2d.farfield(ms, k, bem2d.solve(ms, k, [0.0], method=method), th)[:, 0]
            row[f"relerr_{method}"] = float(np.linalg.norm(F - Fe) / np.linalg.norm(Fe))
            if method == "bm" and ka == 10.0:
                keep = (F, Fe)
        rows.append(row)
    OUT["cylinder"] = rows
    return th, keep


# ---------------------------------------------------------------- 2. fictitious resonance of the flat P1S box
def check_resonance():
    poly = geo.flat("P1S", 1)
    th = np.arange(-90, 91, 5.0)
    fD = C / 2 * np.hypot(1 / 0.231, 1 / 0.134)
    fs = np.arange(1440.0, 1521.0, 2.0)
    res = dict(f_dirichlet=fD, f=fs.tolist(), err_cbie=[], err_bm=[])
    for f in fs:
        ref, _ = solve_ff(poly, f, [0.0], th, h=bem2d.element_size(f) / 2)
        for m in ("cbie", "bm"):
            F, n = solve_ff(poly, f, [0.0], th, method=m)
            res[f"err_{m}"].append(float(np.linalg.norm(F - ref) / np.linalg.norm(ref)))
    i = int(np.argmax(res["err_cbie"]))
    res["worst_cbie"] = dict(f=float(fs[i]), err=res["err_cbie"][i], err_bm_there=res["err_bm"][i])
    res["max_err_bm"] = float(max(res["err_bm"]))
    OUT["resonance_P1S_box"] = res
    return res


# ---------------------------------------------------------------- 3. flat box against Kirchhoff
def kirchhoff_flat(W, f, theta, psi):
    k = 2 * np.pi * f / C
    s = np.sin(np.radians(theta)) + np.sin(np.radians(psi))
    return 2 * k * np.cos(np.radians(theta)) * W * np.abs(np.sinc(k * W * s / (2 * np.pi)))


def check_flat():
    th = np.arange(-90, 90.01, 0.25)
    rows, curves = [], {}
    for ncol in (1, 7):
        W = 0.231 * ncol
        poly = geo.centre_x(geo.flat("P1S", ncol))
        for f in (200.0, 1000.0, 5000.0):
            F, n = solve_ff(poly, f, [0.0, 45.0], th)
            for j, psi in enumerate((0.0, 45.0)):
                K = kirchhoff_flat(W, f, th, psi)
                pk = 2 * 2 * np.pi * f / C * W * np.cos(np.radians(psi))
                spec = int(np.argmin(np.abs(th + psi)))
                LB, LK = db(F[:, j], pk), db(K, pk)
                lobe = (np.abs(np.sin(np.radians(th)) + np.sin(np.radians(psi))) < 0.5 * C / f / W) & (LK > -20)
                rows.append(dict(width_mm=round(W * 1000), f=f, psi=psi, kW=2 * np.pi * f / C * W, elements=n,
                                 bem_specular_re_kirchhoff_db=float(LB[spec]),
                                 max_abs_diff_half_main_lobe_db=float(np.max(np.abs(LB - LK)[lobe])) if lobe.sum() else None))
                curves[(ncol, f, psi)] = (LB, LK)
    OUT["flat_vs_kirchhoff"] = rows
    return th, curves


# ---------------------------------------------------------------- 4. standard QRD of Cox, Angus and D'Antonio
def lobe_energies(th, F, f, P, model, psi=0.0):
    """|R_q|^2 of every propagating grating order q (|sin theta_q| <= 0.95, grazing orders left out), from the energy
    between the half-way points to the neighbouring orders. Fourier has no obliquity factor: |R_q|^2 ~ int |F|^2 dsin.
    BEM is physical: a lobe carries power ~ |R_q|^2 cos theta_q, so |R_q|^2 ~ int |F|^2 dtheta / cos theta_q."""
    lam = C / f
    s = np.sin(np.radians(th)) + np.sin(np.radians(psi))
    E = np.abs(F) ** 2
    out = []
    for q in range(-int(2 / (lam / P)) - 1, int(2 / (lam / P)) + 2):
        sq = q * lam / P
        st = sq - np.sin(np.radians(psi))
        if abs(st) > 0.95:
            continue
        win = np.abs(s - sq) < 0.5 * lam / P
        if model == "fourier":
            e = np.trapezoid(E[win], np.sin(np.radians(th[win])))
        else:
            e = np.trapezoid(E[win], np.radians(th[win])) / np.sqrt(1 - st * st)
        out.append((q, float(e)))
    return out


def check_cox():
    nper = 6
    poly = geo.centre_x(geo.wall("COX", [geo.COX_TILE] * nper))
    P = 0.525
    f0 = geo.qrd.deep_band("COX")["f0"]
    th = np.arange(-90, 90.01, 0.1)
    res, curves = [], {}
    for m in range(1, 8):
        f = m * f0
        F, n = solve_ff(poly, f, [0.0], th)
        fm = geo.fourier("COX", [f], th, [0.0])
        Ff = fm.wall([geo.COX_TILE] * nper)[0, 0]
        Fl = fm.flat(nper)[0, 0]
        Fb_flat, _ = solve_ff(geo.centre_x(geo.flat("COX", nper)), f, [0.0], th)
        row = dict(m=m, f=f, elements=n)
        for name, FF in (("bem", F[:, 0]), ("fourier", Ff)):
            le = lobe_energies(th, FF, f, P, name)
            e = np.array([v for _, v in le])
            row[f"lobes_{name}"] = [q for q, _ in le]
            row[f"lobe_db_re_mean_{name}"] = (10 * np.log10(e / e.mean())).round(2).tolist()
            row[f"lobe_spread_db_{name}"] = float(10 * np.log10(e.max() / e.min()))
            spec = [v for q, v in le if q == 0][0]
            row[f"specular_share_{name}"] = float(spec / e.sum())
        res.append(row)
        curves[m] = (db(F[:, 0], np.abs(Fb_flat[:, 0]).max()), db(Ff, np.abs(Fl).max()))
    OUT["cox_qrd"] = dict(periods=nper, f0=f0, rows=res)
    return th, curves


# ---------------------------------------------------------------- 5. mesh convergence
def check_convergence():
    th = np.arange(-90, 91, 5.0)
    rows = []
    b = geo.qrd.deep_band("P1S")
    cases = [("tile 1", geo.wall("P1S", [1]), f) for f in (b["f0"], 2 * b["f0"], 7 * b["f0"], 10000.0)]
    cases += [("flat box 1 tile", geo.flat("P1S", 1), f) for f in (b["f0"], 10000.0)]
    cases += [("wall 7", geo.wall("P1S", [6, 3, 4, 2, 5, 7, 1]), f) for f in (b["f0"], 10000.0)]
    for name, poly, f in cases:
            poly = geo.centre_x(poly)
            bem2d.GX, bem2d.GW = np.polynomial.legendre.leggauss(4)
            F1, n1 = solve_ff(poly, f, [0.0, 45.0], th)
            bem2d.GX, bem2d.GW = np.polynomial.legendre.leggauss(8)
            F2, n2 = solve_ff(poly, f, [0.0, 45.0], th, h=min(C / f / 20, 0.001))
            bem2d.GX, bem2d.GW = np.polynomial.legendre.leggauss(4)
            E1, E2 = np.abs(F1) ** 2, np.abs(F2) ** 2
            L1, L2 = 10 * np.log10(E1 / E1.max(0)), 10 * np.log10(E2 / E2.max(0))
            mask = L2 > -20
            d = lambda E: ((E.sum(0) ** 2 - (E ** 2).sum(0)) / ((len(E) - 1) * (E ** 2).sum(0)))
            rows.append(dict(case=name, f=round(f, 1), elements=[n1, n2],
                             max_abs_diff_db_above_minus20=float(np.max(np.abs(L1 - L2)[mask])),
                             diffusion_coef=[d(E1).round(4).tolist(), d(E2).round(4).tolist()]))
    OUT["convergence_P1S"] = rows


# ---------------------------------------------------------------- figures
def figures(cyl, flat, cox):
    th, (Fb, Fe) = cyl
    fig = plt.figure(figsize=(15, 5.4))
    ax = fig.add_subplot(131, projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ref = np.abs(Fe).max()
    ax.plot(np.radians(th), np.maximum(db(Fe, ref), -40), color=style.INK2, lw=4, alpha=0.35, label="exact series")
    ax.plot(np.radians(th), np.maximum(db(Fb, ref), -40), color=style.BLUE, lw=1.6, label="BEM")
    ax.set_rlim(-40, 3); ax.set_rticks([-30, -20, -10, 0]); ax.set_yticklabels(["-30", "-20", "-10", "0 dB"], fontsize=10)
    ax.set_title("Rigid cylinder, ka = 10\nscattered far field", pad=14)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.15, -0.18))
    th2, curves = flat
    for i, (f, psi) in enumerate(((200.0, 0.0), (1000.0, 45.0))):
        ax = style.polar_axes(fig, 132 + i)
        LB, LK = curves[(7, f, psi)]
        ax.plot(np.radians(th2), np.maximum(LK, -40), color=style.INK2, lw=4, alpha=0.35, label="Kirchhoff (asymptotic)")
        ax.plot(np.radians(th2), np.maximum(LB, -40), color=style.BLUE, lw=1.6, label="BEM")
        ax.set_title(f"Flat box 1617 x 134 mm\n{f:g} Hz, incidence {psi:g}°", pad=14)
        if i == 1:
            ax.legend(loc="lower left", bbox_to_anchor=(-0.15, -0.18))
    fig.savefig(OUTDIR / "fig-validation-solver.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    th, curves = cox
    fig = plt.figure(figsize=(15, 9.5))
    for i, m in enumerate((1, 2, 3, 4, 6, 7)):
        ax = style.polar_axes(fig, 231 + i)
        LB, LF = curves[m]
        ax.plot(np.radians(th), np.maximum(LF, -40), color=style.BLUE, lw=1.4, ls=(0, (4, 2)), label="Fourier model")
        ax.plot(np.radians(th), np.maximum(LB, -40), color=style.BLUE, lw=1.8, label="BEM")
        f = m * OUT["cox_qrd"]["f0"]
        ax.set_title(f"{m} x f0 = {f:.0f} Hz", pad=12)
        if i == 0:
            ax.legend(loc="lower left", bbox_to_anchor=(-0.2, -0.2))
    fig.suptitle("Standard N = 7 QRD (Cox, Angus, D'Antonio), 6 periods = 3.15 m, normal incidence\n"
                 "dB re the peak of a flat box of the same size, each model against its own flat plate", fontsize=14)
    fig.savefig(OUTDIR / "fig-validation-cox-qrd.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    cyl = check_cylinder(); print("cylinder", OUT["cylinder"], flush=True)
    r = check_resonance(); print("resonance", r["worst_cbie"], r["max_err_bm"], flush=True)
    flat = check_flat(); print("flat done", flush=True)
    cox = check_cox(); print("cox done", flush=True)
    check_convergence(); print("convergence done", flush=True)
    (OUTDIR / "validation.json").write_text(json.dumps(OUT, indent=1), encoding="utf-8")
    figures(cyl, flat, cox)
