"""Cross-sections (2D polygons, metres, counter-clockwise) of the deep tiles, the walls made of them and flat references.

Geometry follows reference/diffuser.py (rule v3): x across the wells, y from the wall (0) to the front plane (T).
Outer walls `edge`, fins `fin`, two tiles side by side make a 2 * edge joint fin; 45 degree chamfers of `chamfer` mm in
the bottom corners of every well; the body is solid from the back to the well floors (the closed air chambers inside
are not part of the exterior problem). 2D means the wells are infinitely long.
Tile parameters and depths come from qrd.py; a temporary key put into qrd.DEEP (as its docstring allows) adds the
standard QRD of Cox, Angus and D'Antonio.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REF = Path(__file__).resolve().parents[2] / "reference"
sys.path.insert(0, str(REF))
import qrd  # noqa: E402
import layout  # noqa: E402

# Standard N = 7 QRD (Cox, Angus, D'Antonio, Proc. IOA 25(5) 2003): f0 490 Hz, depth unit 50 mm, wells 75 mm,
# period 525 mm. Thin fins: 1.6 mm like ours (edge 0.8, so the joint between periods is 1.6 too), no chamfers.
qrd.DEEP["COX"] = dict(N=7, pitch=75.0, length=1000.0, dmax=200.0, fin=1.6, edge=0.8, floor=1.6, back=2.0, cap=1.2,
                       chamfer=0.0, bed="Cox, Angus, D'Antonio 2003")
COX_TILE = 4          # shift 0: depths 0 1 4 2 2 4 1


def front_profile(size, order):
    """Front polyline from (0, T) to (n W, T) in mm for tiles `order` placed left to right in STL X."""
    p = qrd.deep_params(size)
    N, pitch, fin, edge, c = p["N"], p["pitch"], p["fin"], p["edge"], p.get("chamfer", 0.0)
    W, _, T = qrd.deep_dims(size)
    pts = [(0.0, T)]
    for j, tile in enumerate(order):
        X0 = j * W
        d = qrd.deep_depths(size, tile)
        for n in range(N):
            xl = X0 + n * pitch + (edge if n == 0 else fin / 2)
            xr = X0 + (n + 1) * pitch - (edge if n == N - 1 else fin / 2)
            if d[n] <= 0:
                continue
            yf = T - d[n]
            cc = min(c, d[n], (xr - xl) / 2)
            pts += [(xl, T), (xl, yf + cc), (xl + cc, yf)] if cc > 0 else [(xl, T), (xl, yf)]
            pts += [(xr - cc, yf), (xr, yf + cc), (xr, T)] if cc > 0 else [(xr, yf), (xr, T)]
    pts.append((len(order) * W, T))
    return np.array(pts)


def wall(size, order):
    """Closed counter-clockwise polygon (m) of a wall of tiles."""
    W, _, T = qrd.deep_dims(size)
    Wt = len(order) * W
    prof = front_profile(size, order)[::-1]                   # right to left along the front
    poly = np.vstack([[[0.0, 0.0], [Wt, 0.0]], prof])
    return poly / 1000.0


def flat(size, ncol):
    """Rigid box of the same envelope (wells filled): ncol * W wide, T deep."""
    W, _, T = qrd.deep_dims(size)
    return np.array([[0, 0], [ncol * W, 0], [ncol * W, T], [0, T]], float) / 1000.0


def plate(width_mm, thick_mm):
    return np.array([[0, 0], [width_mm, 0], [width_mm, thick_mm], [0, thick_mm]], float) / 1000.0


def cylinder(radius, nseg):
    a = np.linspace(0, 2 * np.pi, nseg, endpoint=False)
    return np.stack([radius * np.cos(a), radius * np.sin(a)], 1)


def centre_x(poly):
    """Shift x so the wall is centred on x = 0 (only the phase reference of the far field moves)."""
    q = poly.copy()
    q[:, 0] -= 0.5 * (q[:, 0].min() + q[:, 0].max())
    return q


# ---------------------------------------------------------------- the Fourier model of layout.py on any frequencies
def fourier(size, f, theta_deg, psi_deg):
    """A layout.Deep instance evaluated at frequencies f (F,), receivers theta (R,), incidences psi (I,).
    Its _tile / shift / wall / flat methods are the exact code that chose the wall orders."""
    m = layout.Deep.__new__(layout.Deep)
    m.size, m.per_band = size, 1
    m.W = qrd.deep_dims(size)[0] / 1000
    m.f = np.asarray(f, float)
    m.theta = np.radians(np.asarray(theta_deg, float))
    k = 2 * np.pi * m.f / qrd.C
    m.k = k[:, None, None]
    m.s = np.sin(m.theta)[None, None, :] + np.sin(np.radians(np.asarray(psi_deg, float)))[None, :, None]
    m.T = {t: m._tile(t) for t in range(1, qrd.DEEP_TILES + 1)}
    return m
