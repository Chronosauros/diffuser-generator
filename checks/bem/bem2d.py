"""2D boundary element solver for sound scattering by a rigid body (exterior Helmholtz, Neumann, plane wave).

Time convention exp(-i w t), free-space Green's function G = (i/4) H0^(1)(k r), lengths in metres.
Boundary: closed polygon, vertices counter-clockwise, outward normal n = (t_y, -t_x) points into the air.
Total field p on the surface from the Burton-Miller equation (no fictitious interior resonances):
    (1/2 I - K + alpha H) p = p_inc - alpha dp_inc/dn,   alpha = i / k
    K p(x) = PV int p dG/dn_y ds,   H p(x) = f.p. int p d2G/dn_x dn_y ds
Constant elements on straight segments, collocation at element midpoints. Near-singular and singular parts are exact:
the Laplace kernels (-1/2pi ln r and its normal derivative) are integrated analytically over every segment, only the
smooth Helmholtz-minus-Laplace remainder goes through Gauss quadrature. H uses the Maue form
    H p(x) = d/ds_x int G dp/ds_y ds_y + k^2 n_x . int G n_y p ds_y,
which for piecewise constant p turns into point values of grad G at the element ends.
Far field: p_s ~ (i/4) sqrt(2 / (pi k R)) exp(i (k R - pi/4)) F(theta), F = sum_e p_e (-i k xhat.n_e) int_e exp(-i k xhat.y) ds.
Angles: theta from the front normal (+y), positive towards +x; a plane wave from direction psi travels along
-(sin psi, cos psi), its specular reflection off a plane y = const goes to theta = -psi.
"""
from __future__ import annotations

import numpy as np
from scipy.special import j0, j1, y0, y1

C = 343.0
NG = 4
GX, GW = np.polynomial.legendre.leggauss(NG)


class Mesh:
    def __init__(self, a, b):
        self.a, self.b = np.asarray(a, float), np.asarray(b, float)
        d = self.b - self.a
        self.L = np.hypot(d[:, 0], d[:, 1])
        self.t = d / self.L[:, None]
        self.n = np.stack([self.t[:, 1], -self.t[:, 0]], 1)
        self.m = 0.5 * (self.a + self.b)

    def __len__(self):
        return len(self.L)


def simplify(poly, tol=1e-12):
    """Drop repeated and collinear vertices of a closed polygon."""
    P = [np.asarray(p, float) for p in poly]
    out = []
    for i, p in enumerate(P):
        if out and np.hypot(*(p - out[-1])) < tol:
            continue
        out.append(p)
    if np.hypot(*(out[0] - out[-1])) < tol:
        out.pop()
    changed = True
    while changed:
        changed = False
        for i in range(len(out)):
            p0, p1, p2 = out[i - 1], out[i], out[(i + 1) % len(out)]
            u, v = p1 - p0, p2 - p1
            if abs(u[0] * v[1] - u[1] * v[0]) < tol * (np.hypot(*u) + np.hypot(*v)) and u @ v > 0:
                out.pop(i)
                changed = True
                break
    return np.array(out)


def area(poly):
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def mesh(poly, h):
    """Split every polygon side into ceil(length / h) equal elements."""
    poly = simplify(poly)
    if area(poly) < 0:
        poly = poly[::-1]
    A, B = [], []
    for i in range(len(poly)):
        p, q = poly[i], poly[(i + 1) % len(poly)]
        n = max(1, int(np.ceil(np.hypot(*(q - p)) / h - 1e-9)))
        s = np.linspace(0, 1, n + 1)[:, None]
        pts = p + (q - p) * s
        A.append(pts[:-1]); B.append(pts[1:])
    return Mesh(np.concatenate(A), np.concatenate(B))


def _Gp(k, r):
    """dG/dr of G = (i/4) H0(k r)."""
    kr = k * r
    return -0.25j * k * j1(kr) + 0.25 * k * y1(kr)


def assemble(ms: Mesh, k: float, rows=None, block=256, alpha=None):
    """K, S = int G ds and H (Maue) for collocation rows (default all), each (rows, N) complex.
    With alpha: only the Burton-Miller matrix 1/2 I - K + alpha H (saves memory on big meshes)."""
    N = len(ms)
    rows = np.arange(N) if rows is None else np.asarray(rows)
    if alpha is not None:
        A = np.empty((len(rows), N), complex)
    else:
        K = np.empty((len(rows), N), complex)
        S = np.empty((len(rows), N), complex)
        H = np.empty((len(rows), N), complex)
    a, b, t, n, L = ms.a, ms.b, ms.t, ms.n, ms.L
    g = 0.5 * (GX + 1)                                         # Gauss points on [0, 1]
    for r0 in range(0, len(rows), block):
        ri = rows[r0:r0 + block]
        X = ms.m[ri]                                           # (B, 2)
        da = a[None] - X[:, None]                              # a_e - x
        ua = np.einsum("bnk,nk->bn", da, t)                    # along the element
        h = -np.einsum("bnk,nk->bn", da, n)                    # (x - a_e) . n_e, signed distance
        self_ = ri[:, None] == np.arange(N)[None]
        h[self_] = 0.0
        ub = ua + L[None]
        ah = np.abs(h)
        hs = np.where(ah > 0, ah, 1.0)
        # Laplace: int ln r ds and int d(ln r)/dn_y ds over the segment, exactly
        def F(u):
            r2 = u * u + h * h
            lg = np.where(r2 > 0, 0.5 * np.log(np.where(r2 > 0, r2, 1.0)), 0.0)
            return u * lg - u + np.where(ah > 0, ah * np.arctan(u / hs), 0.0)
        ilog = F(ub) - F(ua)
        ang = np.where(ah > 0, np.sign(h) * (np.arctan(ub / hs) - np.arctan(ua / hs)), 0.0)
        S_L = -ilog / (2 * np.pi)
        K_L = ang / (2 * np.pi)
        # Helmholtz minus Laplace, smooth, by Gauss
        u = ua[..., None] + g[None, None] * L[None, :, None]   # (B, N, NG)
        r = np.sqrt(u * u + (h * h)[..., None])
        kr = k * r
        R0 = 0.25j * j0(kr) - 0.25 * y0(kr) + np.log(r) / (2 * np.pi)
        R1 = -0.25j * k * j1(kr) + 0.25 * k * y1(kr) + 1 / (2 * np.pi * r)
        w = 0.5 * GW[None, None] * L[None, :, None]
        S_blk = S_L + np.sum(w * R0, -1)
        K_blk = K_L + np.sum(w * R1 * (-h[..., None]) / r, -1)
        # Maue: tangential derivative of G at the element ends
        db = b[None] - X[:, None]
        ra = np.hypot(da[..., 0], da[..., 1]); rb = np.hypot(db[..., 0], db[..., 1])
        tx = t[ri][:, None, :]
        ca = -np.sum(da * tx, -1) / ra                          # t_x . (x - a) / r
        cb = -np.sum(db * tx, -1) / rb
        nn = n[ri] @ n.T
        H_blk = _Gp(k, ra) * ca - _Gp(k, rb) * cb + k * k * nn * S_blk
        if alpha is not None:
            A[r0:r0 + block] = 0.5 * self_ - K_blk + alpha * H_blk
        else:
            K[r0:r0 + block], S[r0:r0 + block], H[r0:r0 + block] = K_blk, S_blk, H_blk
    return A if alpha is not None else (K, S, H)


def incident(ms: Mesh, k, psi_deg):
    """Plane wave from direction psi (degrees from the front normal): values and normal derivatives at midpoints."""
    psi = np.radians(np.atleast_1d(psi_deg))
    d = -np.stack([np.sin(psi), np.cos(psi)], 1)              # propagation directions (P, 2)
    ph = np.exp(1j * k * ms.m @ d.T)                           # (N, P)
    dn = 1j * k * (ms.n @ d.T) * ph
    return ph, dn


def solve(ms: Mesh, k, psi_deg, method="bm", alpha=None):
    """Surface pressure (N, P) for incidences psi_deg. method 'bm' (Burton-Miller) or 'cbie'."""
    pi, dpi = incident(ms, k, psi_deg)
    if method == "bm":
        al = 1j / k if alpha is None else alpha
        return np.linalg.solve(assemble(ms, k, alpha=al), pi - al * dpi)
    K, S, H = assemble(ms, k)
    return np.linalg.solve(0.5 * np.eye(len(ms)) - K, pi)


def farfield(ms: Mesh, k, p, theta_deg):
    """F(theta) (T, P) of the scattered field for surface pressures p (N, P)."""
    th = np.radians(np.atleast_1d(theta_deg))
    xh = np.stack([np.sin(th), np.cos(th)], 1)                 # (T, 2)
    ph = np.exp(-1j * k * xh @ ms.m.T)                         # (T, N)
    sc = np.sinc(k * ms.L[None] * (xh @ ms.t.T) / (2 * np.pi))
    kern = (-1j * k) * (xh @ ms.n.T) * ms.L[None] * sc * ph
    return kern @ p


def element_size(f, per_wavelength=10, hmax=0.002):
    """Elements must not be longer than the fins are thick: with 8 mm elements on 1.6 mm fins the diffusion coefficient
    of a P1S tile was 13 % off at f0, from 2 mm down it converges to 0.1 % (validate.py, convergence_P1S)."""
    return min(C / f / per_wavelength, hmax)
