"""STL generator of the Soundblooms diffuser tiles (millimetres, print orientation). Depths and bands: qrd.py.

    python design/diffuser.py                 v3 set, 45 tiles (3 beds) -> design/out-v3/*.stl + *.json, manifest.json
    python design/diffuser.py deep A1mini 1   one v3 tile -> design/out-v3/ (sizes: A1mini, MK4, P1S; see qrd.py)
    python design/diffuser.py --rule 2 DIR    the v2 set (deep fins 1.2 / edge 0.8) -> DIR, byte for byte the files of
                                              design/out-v2/ (frozen archive: DIR must be another folder)
    python design/diffuser.py --stage1 DIR    the 29 tiles of stage 1 (full-bed size, keys S / L) -> DIR (an empty
                                              folder of your own; design/out/ is the frozen canon of the generator test)

Rule v3 (current, owner decision 22.10.2026) thickens the deep walls: fin 1.6, edge 2.0 at the default depths (qrd.derive);
the face of the tile, the depths and every file name are the same as in v2. Blocks are identical in v2 and v3.

Shifts of the block tiles are always the stage 1 ones in design/out/layouts.json: a v2 size copies the layout of the
class qrd.layout_class() names (no new search).

Deep series (1D QRD, wells with fins), printed standing on its end: X across the wells (0..W), Y from the wall (0) to the
front plane (T), Z along the wells (0 on the bed). Every layer is the same cross-section, so there are no overhangs and
no bridges in the wells. Under each shallow well floor there is a closed air chamber: a solid cap on the bed and a
ridge roof at the top (two 45 degree slopes meeting in a line), so the chamber closes without a bridge
(tools/overhangs.py measures it). Tile number and an up arrow are sunk 0.4 mm into the back.

Block series (2D QRD, square blocks, no fins), printed lying on its back: X (0..W), Y (0..H, up on the wall),
Z from the wall. Blocks are solid in the STL; the slicer prints them with two walls and lightning infill.
The tile number is sunk into the bed face, mirrored so it reads from behind.

Construction is plain: rectangles and triangles -> CrossSection union -> extrude -> boolean with boxes. The Stage 2
generator ports these functions 1:1 (manifold-3d); equality is checked on heights, bounding box and volume.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import trimesh
from manifold3d import CrossSection, JoinType, Manifold, OpType

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qrd  # noqa: E402

STAGE1_LAYOUTS = Path(__file__).resolve().parent / "out" / "layouts.json"   # frozen stage 1 canon, read only
OUT = Path(__file__).resolve().parent / "out-v3"                              # save() writes here; scripts repoint it
LABEL_DEPTH = 0.4
SEGMENTS = {1: "bc", 2: "abged", 3: "abgcd", 4: "fgbc", 5: "afgcd", 6: "afgecd", 7: "abc", 8: "abcdefg", 9: "abcdfg"}


def rect(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def xs(polys):
    return CrossSection(polys) if polys else CrossSection()


def union2d(polys):
    return CrossSection.batch_boolean([CrossSection([p]) for p in polys], OpType.Add)


# ---------------------------------------------------------------- label: 7-segment digit with an up arrow
def label_polys(digit: int, cx: float, cy: float, h: float = 24.0, w: float = 13.0, s: float = 3.0):
    """Polygons in a 2D frame (u right, v up) centred at (cx, cy): arrow above, digit below."""
    seg = {"a": (0, h - s, w, h), "d": (0, 0, w, s), "g": (0, h / 2 - s / 2, w, h / 2 + s / 2),
           "f": (0, h / 2, s, h), "e": (0, 0, s, h / 2), "b": (w - s, h / 2, w, h), "c": (w - s, 0, w, h / 2)}
    x0, y0 = cx - w / 2, cy - (h + 10) / 2
    polys = [rect(x0 + a, y0 + b, x0 + c, y0 + d) for a, b, c, d in (seg[k] for k in SEGMENTS[digit])]
    ay = y0 + h + 4                                  # arrow 4 mm above the digit, 10 wide, 6 tall
    polys.append([(cx - 5, ay), (cx + 5, ay), (cx, ay + 6)])
    return polys


# ---------------------------------------------------------------- deep 1D QRD tile
def deep_tile(size: str, tile: int) -> tuple[Manifold, dict]:
    p = qrd.deep_params(size)
    N, pitch, Lw = p["N"], p["pitch"], p["length"]
    fin, edge, floor, back, cap = p["fin"], p["edge"], p["floor"], p["back"], p["cap"]
    W, _, T = qrd.deep_dims(size)
    d = qrd.deep_depths(size, tile)
    dmax = p["dmax"]

    def well_x(n):          # inner faces of the walls around well n
        xl = n * pitch + (edge if n == 0 else fin / 2)
        xr = (n + 1) * pitch - (edge if n == N - 1 else fin / 2)
        return xl, xr

    polys = [rect(0, 0, W, back), rect(0, 0, edge, T), rect(W - edge, 0, W, T)]
    polys += [rect(j * pitch - fin / 2, 0, j * pitch + fin / 2, T) for j in range(1, N)]
    for n in range(N):
        if d[n] < dmax:
            polys.append(rect(n * pitch, 0, (n + 1) * pitch, T - d[n]))
        if d[n] > 0:                                   # 45 degree chamfers in both bottom corners of the well
            xl, xr = well_x(n)
            yf, c = T - d[n], p["chamfer"]
            polys.append([(xl, yf), (xl + c, yf), (xl, yf + c)])
            polys.append([(xr, yf), (xr, yf + c), (xr - c, yf)])
    body = union2d(polys).extrude(Lw)

    # Optional horizontal dividers (parameters "dividers" = count, "divider" = thickness; absent or 0 = none, the output
    # stays byte for byte the same). A divider is a flat plate in the XY plane inside every well with depth > 0: in X
    # between the inner faces of the two walls, in Y from the well floor to the front plane, so its front is flush with
    # the fin fronts. Centres at (i + 0.5) * length / count: half a gap at both tile ends, so stacked tiles keep the
    # rhythm. Sides and bottom grow into the walls and the floor (overlap, no shared planes), the front plane is shared
    # with the fins on purpose (flush, same outward normal).
    ndiv, tdiv = int(p.get("dividers", 0)), float(p.get("divider", 1.2))
    if ndiv > 0:
        ov = 0.3
        plates = [rect(well_x(n)[0] - ov, T - d[n] - 0.5, well_x(n)[1] + ov, T) for n in range(N) if d[n] > 0]
        slab = union2d(plates).extrude(tdiv)
        divs = [slab.translate([0, 0, (i + 0.5) * Lw / ndiv - tdiv / 2]) for i in range(ndiv)]
        body = Manifold.batch_boolean([body] + divs, OpType.Add)

    chambers = []
    for n in range(N):
        y0, y1 = back, T - d[n] - floor
        if d[n] >= dmax or y1 - y0 < 4:
            continue
        xl, xr = well_x(n)
        cf = (xr - xl) / 2
        z0, z1 = cap, Lw - cap                         # ridge roof: two 45 degree slopes meet, nothing to bridge
        prof = [(xl, z0), (xr, z0), (xr, z1 - cf), ((xl + xr) / 2, z1), (xl, z1 - cf)]
        # profile in the (x, z) plane, extruded along y: frame (u, v, w) -> rotate X 90 -> (u, -w, v) -> shift y
        chambers.append(CrossSection([prof]).extrude(y1 - y0).rotate([90, 0, 0]).translate([0, y1, 0]))
    lab = union2d(label_polys(tile, W / 2, Lw / 2)).extrude(LABEL_DEPTH + 0.1).rotate([90, 0, 0]).translate(
        [0, LABEL_DEPTH, 0])
    solid = body - Manifold.batch_boolean(chambers + [lab], OpType.Add)

    band = qrd.deep_band(size)
    info = dict(series="deep", size=size, tile=tile, shift=qrd.deep_shift(tile), sequence=qrd.shifted(N, qrd.deep_shift(tile)),
                depths=d, dims_whd=[W, Lw, T], stl_axes="X = width across wells, Y = depth from wall, Z = along wells",
                stl_extent_xyz=[W, T, Lw], front_plane_y=T,
                wells=[dict(n=n, x_center=round((well_x(n)[0] + well_x(n)[1]) / 2, 4), floor_y=T - d[n]) for n in range(N)],
                band=band, bed=p["bed"], params=p, upside_down_is_tile=(9 - tile) if tile != 1 else 1)
    if ndiv > 0:
        info.update(dividers=ndiv, divider_z=[round((i + 0.5) * Lw / ndiv, 3) for i in range(ndiv)])
    return solid, info


# ---------------------------------------------------------------- 2D QRD block tile
def block_tile(size: str, tile: int, k: int, l: int) -> tuple[Manifold, dict]:
    p = qrd.block_params(size)
    Nx, Ny, W, H, base = p["Nx"], p["Ny"], p["W"], p["H"], p["base"]
    px, py = W / Nx, H / Ny
    h = qrd.block_heights(size, k, l)
    # No two solids may touch along a plane (coplanar unions leave zero-volume sheets in the STL): cells of equal
    # height are merged in 2D first, every region grows by EPS so neighbours overlap, the base and the regions start
    # below/inside each other, and one clip box cuts the exact outline at the end.
    eps = 0.01
    parts = [Manifold.cube([W + 2, H + 2, base + 1]).translate([-1, -1, -1])]
    for lv in sorted({round(v, 6) for row in h for v in row if v > 1e-6}):
        cells = [rect(n * px, m * py, (n + 1) * px, (m + 1) * py)
                 for m in range(Ny) for n in range(Nx) if abs(h[m][n] - lv) < 1e-6]
        parts.append(union2d(cells).offset(eps, JoinType.Miter).extrude(base + lv - 1).translate([0, 0, 1]))
    body = Manifold.batch_boolean(parts, OpType.Add) ^ Manifold.cube([W, H, base + p["hmax"] + 1])
    # label in the bed face, seen from behind (viewer below, up = +Y): u -> -x, mirrored
    polys = [[(W - u, v) for u, v in poly][::-1] for poly in label_polys(tile, W / 2, H / 2)]
    lab = union2d(polys).extrude(LABEL_DEPTH + 0.1).translate([0, 0, -0.1])
    solid = body - lab
    band = qrd.block_band(size)
    info = dict(series="blocks", size=size, tile=tile, shift=[k, l], dims_whd=[W, H, base + p["hmax"]],
                stl_axes="X = width, Y = height on the wall, Z = depth from wall", stl_extent_xyz=[W, H, base + p["hmax"]],
                base=base, cells=[dict(n=n, m=m, x_center=(n + 0.5) * px, y_center=(m + 0.5) * py, top_z=base + h[m][n])
                                  for m in range(Ny) for n in range(Nx)],
                band=band, bed=p["bed"], params=p)
    return solid, info


def save(solid: Manifold, info: dict, stem: str):
    m = solid.to_mesh()
    tm = trimesh.Trimesh(vertices=np.asarray(m.vert_properties)[:, :3], faces=np.asarray(m.tri_verts), process=False)
    OUT.mkdir(parents=True, exist_ok=True)
    tm.export(OUT / f"{stem}.stl")
    info.update(file=f"{stem}.stl", volume_cm3=round(solid.volume() / 1000, 2), triangles=len(tm.faces))
    (OUT / f"{stem}.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    return info


def block_shifts(size: str) -> list:
    """Shifts of the block tiles of a size: stage 1 sizes use their own entry of the frozen stage 1 layouts.json
    (chosen by layout.py), a v2 size copies the entry of the class qrd.layout_class() names."""
    d = json.loads(STAGE1_LAYOUTS.read_text(encoding="utf-8"))["blocks"]
    key = size if size in d else qrd.layout_class("blocks", qrd.block_params(size))
    return d[key]["shifts"]


def sets(stage1: bool) -> list[tuple[str, str, int]]:
    """(series, size key, tile number) of a whole set: stage 1 = 29 tiles, v2 = 45 tiles."""
    deep, blocks = (qrd.DEEP, qrd.BLOCKS) if stage1 else (qrd.deep_sizes(), qrd.BLOCKS_V2)
    return ([("deep", sz, t) for sz in deep for t in range(1, qrd.DEEP_TILES + 1)]
            + [("blocks", sz, t) for sz in blocks for t in range(1, len(block_shifts(sz)) + 1)])


def build(series: str, size: str, tile: int) -> dict:
    if series == "deep":
        s, i = deep_tile(size, tile)
    else:
        k, l = block_shifts(size)[tile - 1]
        s, i = block_tile(size, tile, k, l)
    return save(s, i, qrd.name(series, size, tile))


def main(argv):
    global OUT
    if argv and argv[0] in ("deep", "blocks"):                       # one tile: series size tile
        print(build(argv[0], argv[1], int(argv[2]))["volume_cm3"])
        return
    stage1 = bool(argv) and argv[0] == "--stage1"
    if stage1:
        if len(argv) < 2:
            sys.exit("--stage1 needs a target folder")
        OUT = Path(argv[1]).resolve()
        if OUT == STAGE1_LAYOUTS.parent.resolve():
            sys.exit("design/out/ is the frozen stage 1 canon; give another folder")
    if argv and argv[0] == "--rule":
        if argv[1:2] != ["2"] or len(argv) < 3:
            sys.exit("--rule 2 DIR regenerates the v2 set into DIR")
        OUT = Path(argv[2]).resolve()
        if OUT in (STAGE1_LAYOUTS.parent.resolve(), STAGE1_LAYOUTS.parent.parent.resolve() / "out-v2"):
            sys.exit("design/out/ and design/out-v2/ are frozen; give another folder")
        qrd.set_rule(2)
    rows = [build(*a) for a in sets(stage1)]
    (OUT / "manifest.json").write_text(json.dumps([dict(file=r["file"], series=r["series"], size=r["size"], tile=r["tile"],
                                                        volume_cm3=r["volume_cm3"]) for r in rows], indent=1), encoding="utf-8")
    for r in rows:
        print(r["file"], r["volume_cm3"], "cm3")


if __name__ == "__main__":
    main(sys.argv[1:])
