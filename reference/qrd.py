"""Sequences, well depths and design band of the Soundblooms diffuser tiles. Pure functions, no geometry.

Source of the formulas: T. J. Cox, P. D'Antonio, "Acoustic Absorbers and Diffusers", 3rd ed., CRC Press 2016,
chapter on Schroeder diffusers (quadratic residue diffuser, QRD). Notation:
    s_n   = n^2 mod N                      quadratic residue sequence, N prime
    d_n   = s_n * lambda0 / (2N)           well depth, lambda0 = c / f0
    f0    = c * s_max / (2 * N * d_max)    design frequency from the deepest well
    f_w   = c / (2 w)                      well width w stops being small against the wavelength
    f_N   = N * f0                         every depth is a whole number of half wavelengths: flat-plate frequency
2D tiles add the phases of two 1D sequences (their reflection factors multiply): phi = frac(s_n / Nx + s_m / Ny),
for Nx = Ny this is the classic (n^2 + m^2) mod N.
All lengths in mm, frequencies in Hz. The Stage 2 generator must port these functions 1:1.
"""
from __future__ import annotations

C = 343.0          # speed of sound, m/s at 20 C
SITE = "Soundblooms"


def residues(N: int) -> list[int]:
    return [(n * n) % N for n in range(N)]


def shifted(N: int, k: int) -> list[int]:
    s = residues(N)
    return [s[(n + k) % N] for n in range(N)]


# ---------------------------------------------------------------- deep series: 1D QRD, N = 7, wells with fins
# pitch = well + one fin; the tile is exactly one period (N wells) wide; tile "length" runs along the wells.
DEEP = {
    "S": dict(N=7, pitch=25.0, length=175.0, dmax=100.0, fin=1.2, edge=0.8, floor=1.6, back=2.0, cap=1.2,
              chamfer=3.0, bed="Bambu Lab A1 mini (180 x 180 x 180)"),
    "L": dict(N=7, pitch=35.0, length=205.0, dmax=140.0, fin=1.2, edge=0.8, floor=1.6, back=2.0, cap=1.2,
              chamfer=4.0, bed="Prusa MK4 (250 x 210 x 220), Core One, Bambu P1S/X1"),
}
DEEP_TILES = 7     # every cyclic shift of the N = 7 sequence is one tile


def deep_shift(tile: int) -> int:
    """Tile 1 is the symmetric period centred on the flat well (2 4 1 0 1 4 2); tile t uses shift (t + 3) mod 7.
    Turning tile t upside down gives tile 9 - t (tile 1 stays tile 1)."""
    return (tile + 3) % 7


def deep_depths(size: str, tile: int) -> list[float]:
    p = deep_params(size)
    s = shifted(p["N"], deep_shift(tile))
    smax = max(residues(p["N"]))
    return [si * p["dmax"] / smax for si in s]


def deep_band(size: str) -> dict:
    p = deep_params(size)
    N, smax = p["N"], max(residues(p["N"]))
    f0 = C * smax / (2 * N * p["dmax"] / 1000)
    w = p["pitch"] - p["fin"]
    f_w = C / (2 * w / 1000)
    f_n = N * f0
    period = N * p["pitch"]
    return dict(f0=f0, f_width=f_w, f_flat=f_n, f_high=min(f_w, f_n), well_width=w, period=period,
                f_period=C / (period / 1000), lambda0=C / f0 * 1000, unit_depth=p["dmax"] / smax)


def deep_dims(size: str) -> tuple[float, float, float]:
    """(width across wells, length along wells, depth from the wall) in mm."""
    p = deep_params(size)
    return p["N"] * p["pitch"], p["length"], p["dmax"] + p["back"]


# ---------------------------------------------------------------- block series: 2D QRD, square blocks, no fins
BLOCKS = {
    "S": dict(Nx=13, Ny=13, W=175.0, H=175.0, hmax=48.0, base=3.0, bed="Bambu Lab A1 mini (180 x 180 x 180)"),
    "L": dict(Nx=13, Ny=11, W=245.0, H=205.0, hmax=60.0, base=3.0,
              bed="Prusa MK4 (250 x 210 x 220), Core One, Bambu P1S/X1"),
}


def phase_fraction(Nx: int, Ny: int, n: int, m: int, k: int = 0, l: int = 0) -> float:
    a = ((n + k) ** 2) % Nx
    b = ((m + l) ** 2) % Ny
    if Nx == Ny:                       # exact integers for the square case
        return ((a + b) % Nx) / Nx
    v = a / Nx + b / Ny
    return v - int(v)


def phase_max(Nx: int, Ny: int) -> float:
    return max(phase_fraction(Nx, Ny, n, m) for n in range(Nx) for m in range(Ny))


def block_heights(size: str, k: int, l: int) -> list[list[float]]:
    """Block heights above the base plate, [m][n] = row m (along H), column n (along W).
    The tallest block is the zero-depth well; depth below the top plane = hmax - height."""
    p = block_params(size)
    Nx, Ny, hmax = p["Nx"], p["Ny"], p["hmax"]
    pm = phase_max(Nx, Ny)
    return [[hmax * (1 - phase_fraction(Nx, Ny, n, m, k, l) / pm) for n in range(Nx)] for m in range(Ny)]


def block_band(size: str) -> dict:
    p = block_params(size)
    pm = phase_max(p["Nx"], p["Ny"])
    f0 = C * pm / (2 * p["hmax"] / 1000)
    w = min(p["W"] / p["Nx"], p["H"] / p["Ny"])
    f_w = C / (2 * w / 1000)
    return dict(f0=f0, f_width=f_w, f_high=f_w, block=w, phase_max=pm, lambda0=C / f0 * 1000,
                f_period=C / (max(p["W"], p["H"]) / 1000))


def block_dims(size: str) -> tuple[float, float, float]:
    p = block_params(size)
    return p["W"], p["H"], p["base"] + p["hmax"]


# ---------------------------------------------------------------- v2 / v3 sizes (owner decisions 22.10.2026)
# Tiles at the full bed size lift off the bed: rules v2 and v3 keep a margin of 25 mm on every bed axis (a tile is at
# least 20 mm smaller than its stage 1 size). Rule v3 (owner decision 22.10.2026, variant A) only thickens the deep walls
# (see derive). derive() is a port of `derive` in generator/CONTRACT.md (the generator still has rule v2 until it is ported
# back). DEEP and BLOCKS above stay the frozen stage 1 canon (generator test, scripts inject parameters into them under a
# temporary key); the deep parameters live in DEEP_V3 (current) and DEEP_V2 (frozen archive, design/out-v2/), the block
# parameters in BLOCKS_V2 (same for both rules), keyed like BEDS. Every function above takes a stage 1 key ("S", "L"), a
# v2 / v3 key ("A1mini", "MK4", "P1S"; the deep ones resolve through RULE, default 3) or a temporary key that a script put
# into DEEP / BLOCKS.
MARGIN = 25.0
BEDS = {                                   # key: (label, x, y, z) in mm
    "A1mini": ("Bambu Lab A1 mini", 180, 180, 180),
    "MK4": ("Prusa MK4 / Core One", 250, 210, 220),
    "P1S": ("Bambu Lab P1S / X1", 256, 256, 256),
}
BED_SHORT = {"A1mini": "A1 mini", "MK4": "MK4 / Core One", "P1S": "P1S / X1"}


def bed_label(key: str) -> str:
    label, x, y, z = BEDS[key]
    return f"{label} ({x} x {y} x {z})"


def ceil01(v: float) -> float:
    """Round up to 0.1 mm (the epsilon keeps an exact multiple, e.g. 2.0, from jumping to 2.1)."""
    import math
    return math.ceil(v * 10 - 1e-9) / 10


def derive(series: str, bed: tuple, depth: float | None = None, label: str | None = None, rule: int = 3) -> dict:
    """Parameters of rule v3 (rule=2: rule v2) for a bed (x, y, z) in mm and a depth in mm (None = default depth).
    The rules differ only in the wall thickness of the deep series: v2 fin 1.2 / edge 0.8 (thin: the tiles swayed when
    printed standing), v3 fin = max(1.6, ceil01(dmax / 82.5)), edge = max(2.0, ceil01(dmax / 66)): thickness grows with
    the depth so that the free walls of a deeper tile (depth slider of the generator) stay as stiff. At the default
    depths (88 / 128 / 132 mm) v3 is always fin 1.6 / edge 2.0. Blocks are the same in both rules."""
    import math
    x, y, z = bed
    label = label or f"custom {x} x {y} x {z}"
    if series == "deep":
        pitch = float((x - MARGIN) // 7)
        back = 2.0
        hi = y - MARGIN - back
        d = min(4 * pitch, hi) if depth is None else depth
        fin, edge = (1.2, 0.8) if rule == 2 else (max(1.6, ceil01(d / 82.5)), max(2.0, ceil01(d / 66)))
        return dict(N=7, pitch=pitch, length=float(min(y - MARGIN, z - MARGIN)), dmax=float(d), fin=fin, edge=edge,
                    floor=1.6, back=back, cap=1.2, chamfer=3.0 if pitch < 30 else 4.0, bed=label)
    W = x - MARGIN
    H = min(y - MARGIN, W)
    Ny = 13 if H / W >= 12 / 13 else 11
    if H / W < 10 / 13:
        W = math.floor(H * 13 / 10 + 0.5)          # JavaScript Math.round: half up
    base = 3.0
    hi = min(z - MARGIN - base, 4 * min(W / 13, H / Ny))
    h = min(48 if Ny == 13 else 60, math.floor(hi)) if depth is None else depth
    return dict(Nx=13, Ny=Ny, W=float(W), H=float(H), hmax=float(h), base=base, bed=label)


DEEP_V2 = {k: derive("deep", b[1:], label=bed_label(k), rule=2) for k, b in BEDS.items()}   # frozen: design/out-v2/
DEEP_V3 = {k: derive("deep", b[1:], label=bed_label(k), rule=3) for k, b in BEDS.items()}   # current: design/out-v3/
BLOCKS_V2 = {k: derive("blocks", b[1:], label=bed_label(k)) for k, b in BEDS.items()}      # blocks: same in v2 and v3

# Rule of the whole toolchain (diffuser.py, layout.py, plan.py, tools, tests): v3. set_rule(2) switches the deep sizes
# back to DEEP_V2 so that the v2 set can be regenerated byte for byte (`diffuser.py --rule 2 DIR`, the regression gate).
RULE = 3


def set_rule(rule: int) -> None:
    global RULE
    if rule not in (2, 3):
        raise ValueError("rule must be 2 or 3")
    RULE = rule


def deep_sizes() -> dict:
    """The v2 / v3 deep sizes of the current rule, keyed A1mini / MK4 / P1S."""
    return DEEP_V3 if RULE == 3 else DEEP_V2


def deep_params(size: str) -> dict:
    return DEEP[size] if size in DEEP else deep_sizes()[size]


def block_params(size: str) -> dict:
    return BLOCKS[size] if size in BLOCKS else BLOCKS_V2[size]


def layout_class(series: str, params: dict) -> str:
    """Which stage 1 wall layout a v2 size copies: deep pitch < 30 -> "S" else "L"; blocks Ny 13 -> "S" else "L"."""
    if series == "deep":
        return "S" if params["pitch"] < 30 else "L"
    return "S" if params["Ny"] == 13 else "L"


def name(series: str, size: str, tile: int) -> str:
    dims = deep_dims(size) if series == "deep" else block_dims(size)
    tag = "qrd-deep" if series == "deep" else "qrd-blocks"
    return f"soundblooms-diffuser-{tag}-{'x'.join(f'{v:g}' for v in dims)}-tile{tile}"


if __name__ == "__main__":
    print(f"v{RULE} sizes (default depth):")
    for sz in deep_sizes():
        print("  deep", sz, deep_dims(sz), {k: round(v, 1) for k, v in deep_band(sz).items()})
        print("  blocks", sz, block_dims(sz), {k: round(v, 3) for k, v in block_band(sz).items()})
    print("stage 1 sizes (archive):")
    for sz in DEEP:
        b = deep_band(sz)
        print("deep", sz, deep_dims(sz), {k: round(v, 1) for k, v in b.items()})
        for t in range(1, DEEP_TILES + 1):
            print("   tile", t, "shift", deep_shift(t), deep_depths(sz, t))
    for sz in BLOCKS:
        b = block_band(sz)
        print("blocks", sz, block_dims(sz), {k: round(v, 3) for k, v in b.items()})
