"""Bass traps, physics: room modes, cell impedance, the single-mode trap model, flags, peak drop. Python twin of the
`bass-physics` core of index.html (createBassPhysics).

Source (read-only research package, copied and trimmed on 25.10.2026, session `diffuser-bass-traps-publish`):
  - mode-traps/room/roomtrap.py (v3: radiation mass with the square root of the number of images): zs_qwr, Cell, room_modes,
    psi, base_decay, sabine_decay, ANCHORS, cluster, cells_psi, loading_mass, solve_single, single_mode_effect, Trap, t60;
  - mode-traps/room/run_room.py: the flags `orient`, `overlap_at`, `no_effect` (base decay > cell gamma), `nodal`;
  - absorbers/tmm.py (constants, slit_eff, dl_slit_array, absorption), absorbers/absorbers.py (qwr_Zin),
    absorbers/design_Q.py (first_res and solve_L on the 30 - 520 Hz grid of 0.25 Hz).
  (all under audio/diffuser/research/2026-10-25-*/ in the private design project).
Changes against the sources: no plotting, no multi-mode model, no brentq cell tuning (the cell length always comes from
solve_L, as the geometry does), anchors have English ids, and everything below "NEW" is new (peak drop in dB, evaluate_*,
best anchor, anchor texts, uncertainty reasons). Convention exp(+j w t); decay <=> Im(w) > 0.
Calculated, not measured: every number here is a model result, see docs/bass-traps.md.

All lengths in mm unless the name says _m; frequencies in Hz; decay rates in 1/s (pressure amplitude).
"""
from __future__ import annotations

import math

import numpy as np

# ---------------------------------------------------------------------------------------------------- tmm.py constants
C0 = 343.2
RHO0 = 1.204
MU = 1.81e-5
KTH = 0.0257
CP = 1005.0
GAM = 1.4
P0 = 101325.0
Z0 = RHO0 * C0
MM = 1e-3
CELL = 0.12                 # cell face, m
CELL_AREA = CELL * CELL
LN1000 = math.log(1000.0)   # 60 dB of pressure amplitude
Y_CELL = 120.0              # cell face width, mm (design_Q.Y)
FG = np.arange(30.0, 520.0, 0.25)      # design_Q.FG: the grid of first_res


def dl_slit_array(h, d, nterms=4000):
    """End correction of a slit array (tmm.dl_slit_array, eq. 15 of arXiv:1708.03343); h slit, d period, both in m."""
    sg = h / d
    m = np.arange(1, nterms + 1)
    return sg * h * np.sum(np.sin(m * np.pi * sg) ** 2 / (m * np.pi * sg) ** 3)


def absorption(zs):
    """Normal-incidence absorption of a surface with specific impedance zs [Pa s/m] (tmm.absorption)."""
    r = (zs - Z0) / (zs + Z0)
    return 1 - np.abs(r) ** 2


_DL_CACHE = {}


def _dl(w_mm, y_mm):
    k = (round(w_mm, 6), round(y_mm, 6))
    if k not in _DL_CACHE:
        _DL_CACHE[k] = dl_slit_array(w_mm * MM, y_mm * MM)
    return _DL_CACHE[k]


# ---------------------------------------------------------------------------------------------------- cell impedance
def slit_eff(w, f, rough=1.0):
    """Effective density and bulk modulus of a slit of width w [m] at real frequencies f (tmm.slit_eff)."""
    f = np.atleast_1d(np.asarray(f, float))
    om = 2 * np.pi * f
    mu = MU * rough
    g = (w / 2) * np.sqrt(1j * om * RHO0 / mu)
    gt = (w / 2) * np.sqrt(1j * om * RHO0 * CP / KTH)
    xi = 1 - np.tanh(g) / g
    psi = 1 - np.tanh(gt) / gt
    return RHO0 / xi, GAM * P0 / (GAM - (GAM - 1) * psi)


def qwr_zin(f, y_mm, w_mm, l_mm, rough=1.0):
    """absorbers.qwr_Zin: acoustic impedance of the straightened slit channel of a QWR cell per unit width at real f."""
    f = np.atleast_1d(np.asarray(f, float))
    rho, k_eff = slit_eff(w_mm * MM, f, rough)
    s = w_mm * MM
    om = 2 * np.pi * f
    k = om * np.sqrt(rho / k_eff)
    zc = np.sqrt(rho * k_eff) / s
    zin = -1j * zc / np.tan(k * (l_mm * MM))
    return zin + 1j * om * RHO0 * dl_slit_array(w_mm * MM, y_mm * MM) / s


def first_res(w, L, rough=1.0):
    """design_Q.first_res: first resonance [Hz] of a straight channel of acoustic length L [mm] and slit w [mm];
    first upward zero crossing of Im Zs on the grid FG, linear interpolation inside the grid cell (NaN when none)."""
    zs = qwr_zin(FG, Y_CELL, w, L, rough) * Y_CELL * MM
    k = np.where(np.diff(np.sign(zs.imag)) > 0)[0]
    if not len(k):
        return float("nan")
    i = k[0]
    y0, y1 = zs.imag[i], zs.imag[i + 1]
    return float(FG[i] + (FG[i + 1] - FG[i]) * (-y0) / (y1 - y0))


def solve_L(f0, w, rough=1.0):
    """design_Q.solve_L: acoustic channel length L [mm] whose first_res is f0 (40 bisection steps)."""
    lo, hi = 0.8 * C0 / (4 * f0) * 1000, 1.02 * C0 / (4 * f0) * 1000
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        fm = first_res(w, mid, rough)
        if fm > f0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def zs_qwr(om, w, L, Y=Y_CELL, rough=1.0):
    """Specific surface impedance [Pa s/m] of the QWR cell at ANY (also complex) angular frequency om (roomtrap.zs_qwr)."""
    om = np.asarray(om, complex)
    mu = MU * rough
    ws = w * MM
    g = (ws / 2) * np.sqrt(1j * om * RHO0 / mu)
    gt = (ws / 2) * np.sqrt(1j * om * RHO0 * CP / KTH)
    xi = 1 - np.tanh(g) / g
    psi = 1 - np.tanh(gt) / gt
    rho = RHO0 / xi
    k_eff = GAM * P0 / (GAM - (GAM - 1) * psi)
    k = om * np.sqrt(rho / k_eff)
    zc = np.sqrt(rho * k_eff) / ws
    zin = -1j * zc / np.tan(k * L * MM)
    zin = zin + 1j * om * RHO0 * _dl(w, Y) / ws
    return zin * Y * MM


class Cell:
    """120 x 120 mm QWR cell, slit w, acoustic channel length L [mm] (from solve_L of the geometry design).
    extra_mass [kg/m2]: radiation mass of the patch, Zs += j w extra_mass. rough = viscosity multiplier of the channel."""

    def __init__(self, f_design, w=6.0, Y=Y_CELL, rough=1.0, extra_mass=0.0, L=None):
        self.f_design, self.w, self.Y, self.rough = f_design, w, Y, rough
        self.L = solve_L(f_design, w, rough) if L is None else float(L)
        self.extra_mass = float(extra_mass)

    def zs(self, om):
        om = np.asarray(om, complex)
        return zs_qwr(om, self.w, self.L, self.Y, self.rough) + 1j * om * self.extra_mass

    def zs_f(self, f):
        return self.zs(2 * np.pi * np.asarray(f, float))

    def dzs(self, om):
        h = 1e-5 * np.abs(om)
        return (self.zs(om + h) - self.zs(om - h)) / (2 * h)

    def alpha_f(self, f):
        return absorption(self.zs_f(f))

    def gamma(self, f=None):
        """Re Zs / |dZs/dw| at f (default the design frequency) [1/s]: half width of the cell resonance (band())."""
        om = 2 * np.pi * (self.f_design if f is None else f)
        z = self.zs(om)
        return float(z.real / abs(self.dzs(om)))


# ---------------------------------------------------------------------------------------------------- room modes
def room_modes(L, fmax=120.0, c=C0):
    """Rigid-wall modes up to fmax: dict(n, f, kind, Lam), sorted by f. L = (Lx, Ly, Lz) in m."""
    Lx, Ly, Lz = L
    V = Lx * Ly * Lz
    nm = [int(math.floor(2 * fmax * Li / c)) for Li in L]
    out = []
    for nx in range(nm[0] + 1):
        for ny in range(nm[1] + 1):
            for nz in range(nm[2] + 1):
                if nx == ny == nz == 0:
                    continue
                f = c / 2 * math.sqrt((nx / Lx) ** 2 + (ny / Ly) ** 2 + (nz / Lz) ** 2)
                if f <= fmax:
                    k = sum(1 for x in (nx, ny, nz) if x > 0)
                    out.append(dict(n=(nx, ny, nz), f=f, kind=("axial", "tangential", "oblique")[k - 1], Lam=V / 2 ** k))
    out.sort(key=lambda m: (m["f"], m["n"]))
    return out


def psi(n, L, pts):
    pts = np.atleast_2d(pts)
    r = np.ones(pts.shape[0])
    for i in range(3):
        r = r * np.cos(n[i] * np.pi * pts[:, i] / L[i])
    return r


def wall_psi2_integral(n, L, wall):
    """Integral of psi_n^2 over a wall: wall area x 2^(-number of nonzero n in the wall plane). Walls 0:x0 1:xL 2:y0 3:yL 4:z0 5:zL."""
    ax = wall // 2
    others = [i for i in range(3) if i != ax]
    area = L[others[0]] * L[others[1]]
    return area * 2.0 ** (-sum(1 for i in others if n[i] > 0))


def beta_wall_from_alpha(alpha):
    """Re(beta) of a wall from its Sabine coefficient (small beta, diffuse field, local reaction): alpha = 8 Re(beta)."""
    return np.asarray(alpha, float) / 8.0


def base_decay(mode, L, alpha_walls, c=C0):
    """Decay rate [1/s] of the pressure amplitude of a mode from the walls alone. alpha_walls: one value or six."""
    aw = np.broadcast_to(np.asarray(alpha_walls, float), (6,))
    bw = beta_wall_from_alpha(aw)
    s = sum(bw[w] * wall_psi2_integral(mode["n"], L, w) for w in range(6))
    return c * s / (2 * mode["Lam"])


def sabine_decay(L, alpha_walls, c=C0):
    aw = np.broadcast_to(np.asarray(alpha_walls, float), (6,))
    sx, sy, sz = L[1] * L[2], L[0] * L[2], L[0] * L[1]
    s = np.array([sx, sx, sy, sy, sz, sz])
    v = L[0] * L[1] * L[2]
    return c * float(np.sum(aw * s)) / (8 * v)


def t60(delta):
    return LN1000 / delta if delta > 0 else float("inf")


# ---------------------------------------------------------------------------------------------------- cell clusters
def wall_axes(wall):
    ax = wall // 2
    others = [i for i in range(3) if i != ax]
    return ax, others[0], others[1]


def _face_pts(center, wall, n=3, size=CELL):
    ax, i1, i2 = wall_axes(wall)
    offs = (np.arange(n) - (n - 1) / 2) * size / n
    pts = []
    for a in offs:
        for b in offs:
            p = np.array(center, float)
            p[i1] += a
            p[i2] += b
            pts.append(p)
    return np.array(pts)


def make_cell(L, wall, u, v):
    ax, i1, i2 = wall_axes(wall)
    c = np.zeros(3)
    c[ax] = 0.0 if wall % 2 == 0 else L[ax]
    c[i1], c[i2] = u, v
    return dict(wall=wall, center=c)


# id: (anchor point as fractions of the dimensions, adjoining walls); research names in the comment
ANCHORS = {
    "corner":  ((0, 0, 0), (0, 2, 4)),       # narożnik
    "edge-z":  ((0, 0, .5), (0, 2)),         # krawędź-z
    "edge-y":  ((0, .5, 0), (0, 4)),         # krawędź-y
    "edge-x":  ((.5, 0, 0), (2, 4)),         # krawędź-x
    "wall-x":  ((0, .5, .5), (0,)),          # środek ściany-x
    "wall-y":  ((.5, 0, .5), (2,)),          # środek ściany-y
    "wall-z":  ((.5, .5, 0), (4,)),          # środek ściany-z
}
IMG_FACTOR = {"corner": 2.0, "edge-x": math.sqrt(2.0), "edge-y": math.sqrt(2.0), "edge-z": math.sqrt(2.0)}   # others 1


def cluster(L, anchor, nmax):
    """The nmax cells closest to the anchor (a prefix of N cells = a trap of N cells), roomtrap.cluster. Grid 0.12 m."""
    P = np.array(ANCHORS[anchor][0]) * np.array(L)
    walls = ANCHORS[anchor][1]
    cands = []
    for wi, wall in enumerate(walls):
        ax, i1, i2 = wall_axes(wall)
        coords = []
        for ii in (i1, i2):
            if P[ii] == 0.0:
                coords.append([CELL * (k + 0.5) for k in range(int(L[ii] / CELL))])
            else:
                kk = int(L[ii] / CELL / 2) - 1
                coords.append([P[ii] + CELL * k for k in range(-kk, kk + 1)])
        for u in coords[0]:
            for v in coords[1]:
                cell = make_cell(L, wall, u, v)
                d = float(np.linalg.norm(cell["center"] - P))
                cands.append((round(d, 9), wi, v, u, cell))
    cands.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
    return [c[4] for c in cands[:nmax]]


def cells_psi(modes, L, cells):
    """Matrix <psi_m>_s (mean over the cell face, 3 x 3 points) [S x M]."""
    P = np.zeros((len(cells), len(modes)))
    for s, cl in enumerate(cells):
        pts = _face_pts(cl["center"], cl["wall"])
        for m, md in enumerate(modes):
            P[s, m] = psi(md["n"], L, pts).mean()
    return P


def loading_mass(N, anchor):
    """Radiation mass [kg/m2] of a patch of N cells (piston of radius a_eq, images enlarge it by sqrt(images)): roomtrap v3."""
    a = math.sqrt(N * CELL_AREA / math.pi)
    return IMG_FACTOR.get(anchor, 1.0) * RHO0 * 8 * a / (3 * math.pi)


# ---------------------------------------------------------------------------------------------------- single-mode equation
def _newton(F, x0, tol=1e-10, maxit=60):
    x = complex(x0)
    h = 1e-6 * max(abs(x), 1.0)
    for _ in range(maxit):
        f0 = F(x)
        d = (F(x + h) - F(x - h)) / (2 * h)
        if d == 0:
            return None
        dx = f0 / d
        x = x - dx
        if abs(dx) < tol * abs(x):
            return x
    return None


def solve_single(mode_w, delta_w, G, cell, c=C0):
    """Roots near the mode of (w^2 - wn^2 - 2 j w dw) Zs(w) = j w c Z0 G; sorted by Im, from the cubic with Zs linearised
    around wn, polished by Newton on the exact Zs."""
    wn = mode_w

    def F(w):
        return (w * w - wn * wn - 2j * w * delta_w) * cell.zs(w) - 1j * w * c * Z0 * G

    z = complex(cell.zs(wn))
    z1 = complex(cell.dzs(wn))
    z0 = z - z1 * wn
    p = np.polymul([1, -2j * delta_w, -wn * wn], [z1, z0])
    p = np.polyadd(p, [-1j * c * Z0 * G, 0])
    roots = np.roots(p)
    out = []
    for r in roots:
        if r.real <= 0 or abs(r.real - wn) > 0.2 * wn or r.imag < -0.2 * wn:
            continue
        rr = _newton(F, r)
        if rr is not None and rr.real > 0 and abs(rr.real - wn) < 0.25 * wn:
            if not any(abs(rr - o) < 1e-6 * wn for o in out):
                out.append(rr)
    out.sort(key=lambda x: x.imag)
    return out


def single_mode_effect(mode, delta_w_eff, G, cell, c=C0):
    """Exact single-mode result (Trap.single): roots, the slower of the two roots nearest the mode as delta_eff."""
    wn = 2 * np.pi * mode["f"]
    roots = solve_single(wn, delta_w_eff, G, cell, c)
    z = complex(cell.zs(wn))
    z1 = complex(cell.dzs(wn))
    omz = wn - z / z1
    g2 = 1j * c * Z0 * G / (2 * z1)
    dz = omz - wn
    lam = 4 * abs(g2) / abs(dz - 1j * delta_w_eff) ** 2
    res = dict(roots=[(r.real / (2 * np.pi), r.imag) for r in roots], g2=g2, Dz=dz, lam=float(lam),
               gamma=float(omz.imag), detune_Hz=float(omz.real / (2 * np.pi) - mode["f"]))
    if roots:
        rr = sorted(roots, key=lambda x: abs(x.real - wn))[:2]
        slow = min(rr, key=lambda x: x.imag)
        res["delta_eff"] = float(slow.imag)
        res["f_eff"] = float(slow.real / (2 * np.pi))
        res["split_Hz"] = float(abs(rr[0].real - rr[1].real) / (2 * np.pi)) if len(rr) > 1 else 0.0
    return res


class Trap:
    """A trap of the first N cells of `cells` (all of the cell type `cell`) in room L, for the modes in `modes`."""

    def __init__(self, L, modes, cells, cell, alpha_walls, P=None, anchor=None):
        self.L, self.modes, self.cells, self.cell, self.aw = L, modes, cells, cell, alpha_walls
        self.anchor, self._cn = anchor, {}
        self.P = cells_psi(modes, L, cells) if P is None else P
        aw = np.broadcast_to(np.asarray(alpha_walls, float), (6,))
        self.bw = np.array([beta_wall_from_alpha(aw[cl["wall"]]) for cl in cells])
        self.dw = np.array([base_decay(m, L, alpha_walls) for m in modes])
        self.lam = np.array([m["Lam"] for m in modes])

    def cell_N(self, N):
        """The cell with the radiation mass of an N-cell patch."""
        if self.anchor is None:
            return self.cell
        if N not in self._cn:
            c = self.cell
            self._cn[N] = Cell(c.f_design, c.w, c.Y, c.rough, L=c.L, extra_mass=loading_mass(N, self.anchor))
        return self._cn[N]

    def Gn(self, n_idx, N):
        return CELL_AREA * float(np.sum(self.P[:N, n_idx] ** 2)) / self.lam[n_idx]

    def delta_w_eff(self, n_idx, N):
        gc = CELL_AREA * self.P[:N, n_idx] ** 2 / self.lam[n_idx]
        return float(self.dw[n_idx] - 0.5 * C0 * np.sum(gc * self.bw[:N]))

    def single(self, n_idx, N):
        return single_mode_effect(self.modes[n_idx], self.delta_w_eff(n_idx, N), self.Gn(n_idx, N), self.cell_N(N))


# ====================================================================================================== NEW
DETUNE_PCT = 0.5            # f cell = f mode x 1.005
SLOT_W = 6.0                # the only public slit width, mm (slot-width-decision.md)
N_LIST = (1, 2, 4, 8, 16)
PEAK_POINTS = 3001
PEAK_SPAN = (0.85, 1.15)


def cell_for_mode(f_mode, L0=None, w=SLOT_W):
    """The cell of the geometry design for a mode: f cell = f mode x 1.005, L0 = solve_L(f cell, w) unless L0 is given."""
    f_cell = f_mode * (1 + DETUNE_PCT / 100.0)
    return Cell(f_cell, w, L=solve_L(f_cell, w) if L0 is None else L0)


def peak_drop_db(mode, delta_base, delta_w_eff, G, cell, c=C0):
    """Drop of the resonance peak of the mode in dB for a force on the mode:
      H_base(w) = 1 / (w^2 - wn^2 - 2 j w d_n),  H_trap(w) = Zs / ((w^2 - wn^2 - 2 j w d') Zs - j w c Z0 G),
    peak = max |H| on w = wn x linspace(0.85, 1.15, 3001); result 20 log10(peak_trap / peak_base), negative = lower peak."""
    wn = 2 * np.pi * mode["f"]
    w = wn * np.linspace(PEAK_SPAN[0], PEAK_SPAN[1], PEAK_POINTS)
    hb = 1.0 / (w * w - wn * wn - 2j * w * delta_base)
    zs = cell.zs(w)
    ht = zs / ((w * w - wn * wn - 2j * w * delta_w_eff) * zs - 1j * w * c * Z0 * G)
    return float(20 * np.log10(np.max(np.abs(ht)) / np.max(np.abs(hb))))


def nearest_other(modes, i):
    """Index of the mode nearest in frequency to mode i (ties: lower index), None for a single mode."""
    best, bi = None, None
    for k, m in enumerate(modes):
        if k == i:
            continue
        d = abs(m["f"] - modes[i]["f"])
        if best is None or d < best:
            best, bi = d, k
    return bi


def overlap_flag(modes, i, L, alpha):
    """run_room.overlap_at: the nearest mode is closer than the sum of the two mode widths (delta / pi each)."""
    o = nearest_other(modes, i)
    if o is None:
        return False
    return bool(abs(modes[o]["f"] - modes[i]["f"]) < (base_decay(modes[i], L, alpha) + base_decay(modes[o], L, alpha)) / np.pi)


def orient(kind, overlap, N):
    """run_room.orient: True = approximate (single-mode model off the multi-mode / FEM results)."""
    return bool(overlap or N > 8 or (kind != "axial" and N >= 4))


def approximate_reasons(kind, overlap, N, anchor, neighbour_in_window):
    r = []
    if overlap:
        r.append("overlapping modes")
    if N > 8:
        r.append("more than 8 cells")
    if kind != "axial" and N >= 4:
        r.append("tangential or oblique mode with 4 or more cells")
    if anchor == "corner" and neighbour_in_window:
        r.append("corner with a neighbouring mode within 10 % of the cell frequency")
    return r


def neighbour_in_window(modes, i, f_cell, window=0.10):
    return any(k != i and abs(m["f"] - f_cell) < window * f_cell for k, m in enumerate(modes))


def uncertain_reasons(kind, alpha_scalar, delta_base, gamma):
    """Reasons the result is shown as uncertain (REVIEW.md 3, slot-width-decision.md): near the no-effect limit, oblique
    modes at alpha >= 0.10, alpha >= 0.15 (the sign of the effect then depends on the unmeasured cell resistance)."""
    r = []
    if delta_base > 0.8 * gamma and delta_base <= gamma:
        r.append("base decay above 0.8 of the cell width")
    if kind == "oblique" and alpha_scalar >= 0.10 - 1e-9:      # mean of six 0.1 is 0.0999...: compare with a tolerance
        r.append("oblique mode at alpha >= 0.10")
    if alpha_scalar >= 0.15 - 1e-9:
        r.append("alpha >= 0.15")
    return r


def fmt_m(x):
    s = repr(round(x * 100) / 100)
    return s if "." in s else s + ".0"


def anchor_text(anchor, L):
    """English description of an anchor position in the room L = (Lx, Ly, Lz) in m."""
    lx, ly, lz = (fmt_m(v) for v in L)
    return {
        "corner": "floor corner (the three walls meet)",
        "edge-z": f"vertical corner edge, half height ({lz} m)",
        "edge-y": f"floor edge under the {ly} x {lz} m wall, middle of the {ly} m edge",
        "edge-x": f"floor edge under the {lx} x {lz} m wall, middle of the {lx} m edge",
        "wall-x": f"middle of the {ly} x {lz} m wall",
        "wall-y": f"middle of the {lx} x {lz} m wall",
        "wall-z": f"middle of the floor ({lx} x {ly} m)",
    }[anchor]


def alpha_scalar(alpha_walls):
    return float(np.mean(np.broadcast_to(np.asarray(alpha_walls, float), (6,))))


def analyze_room(L, alpha_walls, fmax=120.0):
    """All modes up to fmax with base decay, T60 before and the overlap flag (the table behind the mode picker)."""
    modes = room_modes(L, fmax)
    rows = []
    for i, m in enumerate(modes):
        d = base_decay(m, L, alpha_walls)
        o = nearest_other(modes, i)
        rows.append(dict(index=i, n=list(m["n"]), f=m["f"], kind=m["kind"], Lam=m["Lam"], delta_base=d, t60_before=t60(d),
                         spacing_Hz=abs(modes[o]["f"] - m["f"]) if o is not None else None,
                         overlap=overlap_flag(modes, i, L, alpha_walls)))
    return modes, rows


def evaluate_mode(L, modes, i, alpha_walls, ns=N_LIST, anchors=None, L0=None, w=SLOT_W):
    """Everything for mode i of `modes` (room L): cell, base decay, flags, and for every anchor and N the exact single-mode
    result with T60 after, decay factor and peak drop. `best[N]` = anchor with the largest delta_eff (ties: first)."""
    mode = modes[i]
    anchors = list(ANCHORS) if anchors is None else list(anchors)
    cell = cell_for_mode(mode["f"], L0, w)
    d_base = float(base_decay(mode, L, alpha_walls))
    gamma = cell.gamma()
    a_s = alpha_scalar(alpha_walls)
    ov = overlap_flag(modes, i, L, alpha_walls)
    nbw = neighbour_in_window(modes, i, cell.f_design)
    out = dict(mode=dict(n=list(mode["n"]), f=mode["f"], kind=mode["kind"], Lam=mode["Lam"]), delta_base=d_base,
               t60_before=t60(d_base), gamma_cell=gamma, f_cell=cell.f_design, L0=cell.L, no_effect=bool(d_base > gamma),
               overlap=ov, uncertain=uncertain_reasons(mode["kind"], a_s, d_base, gamma), anchors={}, best={})
    nmax = max(ns)
    for an in anchors:
        cl = cluster(L, an, nmax)
        P = cells_psi([mode], L, cl)
        tr = Trap(L, [mode], cl, cell, alpha_walls, P=P, anchor=an)
        rec = dict(text=anchor_text(an, L), psi2=float(P[0, 0] ** 2), nodal=bool(P[0, 0] ** 2 < 0.05), cells_available=len(cl), by_n={})
        for N in ns:
            if N > len(cl):
                rec["by_n"][str(N)] = dict(available=False)
                continue
            s = tr.single(0, N)
            dw_eff, G = tr.delta_w_eff(0, N), tr.Gn(0, N)
            pk = peak_drop_db(mode, d_base, dw_eff, G, tr.cell_N(N))
            reasons = approximate_reasons(mode["kind"], ov, N, an, nbw)
            de = s.get("delta_eff")
            rec["by_n"][str(N)] = dict(
                available=True, solved=de is not None, delta_eff=de, f_eff=s.get("f_eff"),
                t60_after=t60(de) if de is not None else None,
                decay_factor=de / d_base if de is not None else None,
                peak_drop_db=pk, delta_walls_eff=dw_eff, G=G,
                approximate=bool(orient(mode["kind"], ov, N) or reasons), approximate_reasons=reasons)
        out["anchors"][an] = rec
    for N in ns:
        best, bd = None, None
        for an in anchors:
            r = out["anchors"][an]["by_n"].get(str(N), {})
            if r.get("available") and r.get("solved") and (bd is None or r["delta_eff"] > bd):
                best, bd = an, r["delta_eff"]
        out["best"][str(N)] = best
    return out
