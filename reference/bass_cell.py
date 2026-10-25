"""Bass traps, geometry: the folded quarter-wave cell (QWR), blocks of cells, tuning inserts and the inserts plate.
Python twin of the `bass-geometry` core of index.html (createBassGeometry).

Source (read-only research package, copied and trimmed on 25.10.2026, session `diffuser-bass-traps-publish`):
  mode-traps/cell/qwr_cell.py (v2 defaults: bend_dL = -0.245 w, f cell = f mode x 1.005, inserts -1 ... +4 % in 0.5 % steps),
  channel length from design_Q.solve_L / first_res (here in bass_room.py), decision w = 6 mm only: mode-traps/slot-width-decision.md.
Changes against the source: no PNG, no TMM summary, no CLI; the CAD part (build123d) is imported inside the export functions only,
so everything else runs with numpy alone; `plan_block` returns the rectangle plan of a block without building any solid.

Coordinates (mm): x = depth from the face (face plane x = 0), y = face width, z = height = print axis (the channel is extruded
along z, so there are no supports). The insert closes the channel (nominal insert k = 0 is a mandatory end cap); other k are
optional tuning parts, to be used only with a decay measurement before and after.
"""
from __future__ import annotations

import json
import math
import os

import bass_room as br

DEF = dict(
    w=6.0,            # slit width [mm]
    wall=1.2,         # wall and floor thickness [mm]
    width=120.0,      # face width Y [mm]
    height=120.0,     # face height H = extrusion [mm]
    bend_dL=None,     # length correction per 90 degree bend; None = BEND_K * w
    detune_pct=0.5,   # f cell = f mode * (1 + detune/100)
    steps_down=2,     # inserts below the nominal
    steps_up=8,       # inserts above the nominal (cut when the last leg is too short)
    step_pct=0.5,     # tuning step [% of f]
    max_depth=150.0,  # warning when the cell is deeper [mm]
    insert_min=10.0,  # shortest insert [mm]
    ins_cx=0.25,      # insert clearance on the slit side [mm]
    ins_cz=0.30,      # insert clearance in height [mm] (per side)
    flange_t=1.6,     # insert flange thickness [mm]
    flange_front=8.0, # the flange reaches this far in front of the slit [mm]
)
BEND_K = -0.245       # length correction per 90 degree bend = BEND_K * w [mm] (2D FEM)


def default_bend_dL(w):
    return BEND_K * w


# ---------------------------------------------------------------------------------------------------- axis and plan
def axis_length(n, y_in, p, w, TW, Y):
    """Length of the axis (base orientation: leg 0 goes +y) from the face to the end of an insert of length p. Returns (L, y_end)."""
    xc0 = TW + w / 2
    yT, yB = Y - TW - w / 2, TW + w / 2
    yc_in = y_in + w / 2
    pitch = w + TW
    if n == 1:
        y_end = Y - p
        return xc0 + (y_end - yc_in), y_end
    L = xc0 + (yT - yc_in) + (n - 1) * pitch + (n - 2) * (yT - yB)
    if (n - 1) % 2 == 0:
        y_end = Y - p
        L += y_end - yB
    else:
        y_end = p
        L += yT - y_end
    return L, y_end


def plan_cell(L_ac, w, TW, Y, bend_dL, p_nom):
    """Number of legs n, inlet position y_in and nominal insert length p (>= p_nom) so that L_acoustic = L_axis + bends * bend_dL = L_ac."""
    for n in range(1, 80):
        bends = 1 + 2 * (n - 1)
        L_geom = L_ac - bends * bend_dL
        needed = L_geom + p_nom
        Lmax = axis_length(n, TW, 0.0, w, TW, Y)[0]
        if Lmax < needed:
            continue
        y_in = TW + (Lmax - needed)
        p = p_nom
        ymax = Y - TW - w
        if y_in > ymax:
            y_in = ymax
            Lmin = axis_length(n, ymax, 0.0, w, TW, Y)[0]
            p = p_nom + (Lmin - needed)
        return dict(n=n, y_in=y_in, p_nom=p, bends=bends, L_geom=L_geom)
    raise ValueError("L too long")


def design(f0=None, L=None, w=DEF["w"], wall=DEF["wall"], width=DEF["width"], height=DEF["height"],
           bend_dL=DEF["bend_dL"], detune_pct=DEF["detune_pct"], steps_down=DEF["steps_down"], steps_up=DEF["steps_up"],
           step_pct=DEF["step_pct"], insert_min=DEF["insert_min"], max_depth=DEF["max_depth"]):
    """L nominal (acoustic), f cell, inserts k -> (L_k, p_k), leg plan. f0 = MODE frequency: f cell = f0 * (1 + detune_pct / 100)."""
    if bend_dL is None:
        bend_dL = default_bend_dL(w)
    f_mode = f0
    if L is None:
        if f0 is None:
            raise ValueError("give f0 or L")
        f0n = f0 * (1 + detune_pct / 100.0)
        L0 = br.solve_L(f0n, w)
    else:
        L0 = float(L)
        f0n = br.first_res(w, L0)
        f_mode = f_mode or f0n
    p_limit = width - wall - w - 3.0
    warn = []
    Lk = {k: (L0 if k == 0 else br.solve_L(f0n * (1 + k * step_pct / 100.0), w)) for k in range(-steps_down, steps_up + 1)}
    p_nom_req = insert_min + (Lk[-steps_down] - L0)
    pl = plan_cell(L0, w, wall, width, bend_dL, p_nom_req)
    p_nom = pl["p_nom"]
    pk = {k: p_nom + (L0 - Lk[k]) for k in Lk}
    up = steps_up
    while up > 0 and pk[up] > p_limit:
        up -= 1
    if up < steps_up:
        warn.append(f"upper tuning range cut from +{steps_up * step_pct:g} % to +{up * step_pct:g} % ({up} steps of {step_pct:g} %): "
                    f"the +{steps_up * step_pct:g} % insert ({pk[steps_up]:.0f} mm) does not fit the last leg (limit {p_limit:.0f} mm)")
        Lk = {k: v for k, v in Lk.items() if k <= up}
        pk = {k: v for k, v in pk.items() if k <= up}
    T = wall + pl["n"] * (w + wall)
    if T > max_depth:
        warn.append(f"cell depth {T:.1f} mm > {max_depth:g} mm")
    return dict(L0=L0, f0=f0n, f_mode=f_mode, plan=pl, Lk=Lk, pk=pk, ns_down=steps_down, ns_up=up, warn=warn, w=w, wall=wall,
                width=width, height=height, bend_dL=bend_dL, step_pct=step_pct, detune_pct=detune_pct if L is None else 0.0,
                p_limit=p_limit, depth=T)


def cell_layout(plan, w, TW, Y, H):
    """Rectangles of air / walls / voids of one cell in the base orientation (leg 0 goes +y)."""
    n, y_in = plan["n"], plan["y_in"]
    T = TW + n * (w + TW)
    z0, z1 = TW, H - TW
    xs = [TW + j * (w + TW) for j in range(n)]
    air, walls, voids = [], [], []

    def R(lst, kind, x0, x1, y0, y1, zz=(z0, z1)):
        lst.append(dict(kind=kind, x=[round(x0, 4), round(x1, 4)], y=[round(y0, 4), round(y1, 4)], z=[round(zz[0], 4), round(zz[1], 4)]))

    R(air, "inlet", 0.0, xs[0], y_in, y_in + w)
    for j in range(n):
        ya = y_in if j == 0 else TW
        yb = Y - TW
        R(air, "leg", xs[j], xs[j] + w, ya, yb)
        if j < n - 1:
            if j % 2 == 0:
                R(air, "turn", xs[j], xs[j + 1] + w, Y - TW - w, Y - TW)
            else:
                R(air, "turn", xs[j], xs[j + 1] + w, TW, TW + w)
    if (n - 1) % 2 == 0:
        R(air, "window", xs[-1], xs[-1] + w, Y - TW, Y)
        wall_y = Y
    else:
        R(air, "window", xs[-1], xs[-1] + w, 0.0, TW)
        wall_y = 0.0
    pocket = (y_in - 2 * TW) >= 1.0
    if pocket:
        R(voids, "pocket", xs[0], xs[0] + w, TW, y_in - TW)
    zz = (0.0, H)
    R(walls, "front", 0.0, TW, 0.0, y_in, zz)
    R(walls, "front", 0.0, TW, y_in + w, Y, zz)
    R(walls, "back", T - TW, T, 0.0, Y, zz)
    for side, (ya, yb) in (("lat-y0", (0.0, TW)), ("lat-ymax", (Y - TW, Y))):
        if (side == "lat-ymax" and wall_y == Y) or (side == "lat-y0" and wall_y == 0.0):
            R(walls, side, 0.0, xs[-1], ya, yb, zz)
            R(walls, side, xs[-1] + w, T, ya, yb, zz)
        else:
            R(walls, side, 0.0, T, ya, yb, zz)
    for j in range(n - 1):
        if j % 2 == 0:
            R(walls, "partition", xs[j] + w, xs[j] + w + TW, TW, Y - TW - w, zz)
        else:
            R(walls, "partition", xs[j] + w, xs[j] + w + TW, TW + w, Y - TW, zz)
    if pocket:
        R(walls, "leg0-end", xs[0], xs[0] + w, y_in - TW, y_in, zz)
    elif y_in > TW + 1e-9:
        R(walls, "leg0-fill", xs[0], xs[0] + w, TW, y_in, zz)
    return dict(T=T, xs=xs, air=air, walls=walls, voids=voids, wall_y=wall_y, pocket=pocket)


def axis_polyline(plan, w, TW, Y):
    n, y_in, p = plan["n"], plan["y_in"], plan["p_nom"]
    xs = [TW + j * (w + TW) + w / 2 for j in range(n)]
    yT, yB = Y - TW - w / 2, TW + w / 2
    yc = y_in + w / 2
    pts = [(0.0, yc), (xs[0], yc)]
    bends = [(xs[0], yc)]
    for j in range(n - 1):
        yt = yT if j % 2 == 0 else yB
        pts += [(xs[j], yt), (xs[j + 1], yt)]
        bends += [(xs[j], yt), (xs[j + 1], yt)]
    _, y_end = axis_length(n, y_in, p, w, TW, Y)
    pts.append((xs[-1], y_end))
    return pts, bends, y_end


def union_area(rects):
    """Exact area of the union of axis-aligned rectangles (coordinate compression)."""
    xs = sorted({v for r in rects for v in r["x"]})
    ys = sorted({v for r in rects for v in r["y"]})
    a = 0.0
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            cx, cy = (xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2
            if any(r["x"][0] <= cx <= r["x"][1] and r["y"][0] <= cy <= r["y"][1] for r in rects):
                a += (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j])
    return a


def poly_len(pts):
    return sum(math.hypot(pts[i + 1][0] - pts[i][0], pts[i + 1][1] - pts[i][1]) for i in range(len(pts) - 1))


def tr_rect(r, mirror, Y, oy, oz):
    x = list(r["x"])
    y0, y1 = r["y"]
    if mirror:
        y0, y1 = Y - y1, Y - y0
    return dict(kind=r["kind"], x=x, y=[round(y0 + oy, 4), round(y1 + oy, 4)], z=[round(r["z"][0] + oz, 4), round(r["z"][1] + oz, 4)])


def tr_pt(pt, mirror, Y, oy):
    y = Y - pt[1] if mirror else pt[1]
    return [round(pt[0], 4), round(y + oy, 4)]


def block_cells(ny, nz, k):
    """Cell positions (iy, iz) of a block ny x nz holding k cells, filled in rows from the bottom (2x2 with k 3 = an L, no overhang)."""
    cells = []
    for iz in range(nz):
        for iy in range(ny):
            if len(cells) < k:
                cells.append((iy, iz))
    return cells


def plan_block(d, ny, nz, k):
    """The block of k cells in an ny x nz grid for design d: size, every cell's rectangles and axis, the expected volume.
    No solid is built (the rectangle plan is the contract for the CAD part and for the JS port)."""
    pl, w, wall, Y, H = d["plan"], d["w"], d["wall"], d["width"], d["height"]
    n = pl["n"]
    base_end = "+y" if n % 2 == 1 else "-y"
    lay = cell_layout(pl, w, wall, Y, H)
    T = lay["T"]
    pts, bends, y_end = axis_polyline(pl, w, wall, Y)
    L_axis = poly_len(pts)
    L_ac = L_axis + pl["bends"] * d["bend_dL"]
    hc = H - 2 * wall
    cells = []
    for idx, (iy, iz) in enumerate(block_cells(ny, nz, k)):
        want = "-y" if (ny > 1 and iy < ny / 2) else "+y"
        if ny == 1:
            want = "+y"
        mirror = (want != base_end)
        oy, oz = iy * Y, iz * H
        air = [tr_rect(r, mirror, Y, oy, oz) for r in lay["air"]]
        walls = [tr_rect(r, mirror, Y, oy, oz) for r in lay["walls"]]
        voids = [tr_rect(r, mirror, Y, oy, oz) for r in lay["voids"]]
        axis_pts = [tr_pt(q, mirror, Y, oy) for q in pts]
        wall_y = Y - lay["wall_y"] if mirror else lay["wall_y"]
        wall_y += oy
        ye = axis_pts[-1][1]
        sgn = 1.0 if ye < wall_y else -1.0
        xl = lay["xs"][-1]
        cells.append(dict(index=idx, index_y=iy, index_z=iz, origin=[0.0, oy, oz], mirror_y=mirror, air=air, walls=walls, voids=voids,
                          axis_points_xy=axis_pts, window=dict(outer_wall_plane_y=round(wall_y, 4), x=[xl, xl + w],
                                                             z=[oz + wall, oz + H - wall], insert_axis_into_cell="+y" if sgn > 0 else "-y")))
    air_area = union_area(lay["air"])
    void_area = union_area(lay["voids"]) if lay["voids"] else 0.0
    vol = (T * Y * H - hc * (air_area + void_area)) * len(cells)
    return dict(T=T, size_mm=[T, ny * Y, nz * H], ny=ny, nz=nz, k=len(cells), cells=cells, x_last=lay["xs"][-1], hc=hc,
                L_axis=L_axis, L_acoustic=L_ac, bends=pl["bends"], legs=n, volume_mm3=vol)


def insert_dims(p, d, T, x_last, hc, H):
    """Dimensions of one insert: flange [fx, fz, ft] and body [bw, bh, p] (print orientation: flange on the bed)."""
    ft, ff, w = DEF["flange_t"], DEF["flange_front"], d["w"]
    cx, cz = DEF["ins_cx"], DEF["ins_cz"]
    fx = T - (x_last - ff)
    fz = H - 1.0
    bw, bh = w - 2 * cx, hc - 2 * cz
    return dict(flange=[fx, fz, ft], body=[bw, bh, p], body_origin=[ff + cx, (fz - bh) / 2, ft])


def plate_plan(d, T, x_last, hc, H):
    """Inserts k = -ns_down ... ns_up side by side along x with 3 mm between them: positions and dims."""
    out, x_cursor = [], 0.0
    for kk in range(-d["ns_down"], d["ns_up"] + 1):
        dims = insert_dims(d["pk"][kk], d, T, x_last, hc, H)
        out.append(dict(k=kk, x0=x_cursor, p=d["pk"][kk], dims=dims))
        x_cursor += dims["flange"][0] + 3.0
    return out


def insert_file_suffix(kk):
    return ("m" if kk < 0 else "p" if kk > 0 else "") + str(abs(kk))


# ---------------------------------------------------------------------------------------------------- CAD (build123d only here)
def _bx(x0, x1, y0, y1, z0, z1):
    from build123d import Box, Pos, Align
    return Pos(x0, y0, z0) * Box(x1 - x0, y1 - y0, z1 - z0, align=(Align.MIN, Align.MIN, Align.MIN))


def build_block(cells, T, Y, H):
    solid = None
    for c in cells:
        oy, oz = c["origin"][1], c["origin"][2]
        h = _bx(0, T, oy, oy + Y, oz, oz + H)
        solid = h if solid is None else solid + h
    for c in cells:
        for r in c["air"] + c["voids"]:
            solid = solid - _bx(r["x"][0], r["x"][1], r["y"][0], r["y"][1], r["z"][0], r["z"][1])
    return solid


def build_insert(p, d, T, x_last, hc, H):
    dm = insert_dims(p, d, T, x_last, hc, H)
    fx, fz, ft = dm["flange"]
    bw, bh, _ = dm["body"]
    bx0, by0, _ = dm["body_origin"]
    ins = _bx(0, fx, 0, fz, 0, ft) + _bx(bx0, bx0 + bw, by0, by0 + bh, ft, ft + p)
    try:
        from build123d import Axis
        edges = ins.edges().filter_by_position(Axis.Z, ft + p - 0.01, ft + p + 0.01)
        ins = ins.chamfer(0.6, None, edges)
    except Exception:
        pass
    return ins


def export_stl(part, path):
    from build123d import export_stl as _e
    _e(part, path, tolerance=0.001, angular_tolerance=0.1)


def export_set(out_dir, name, d, bp, inserts=True):
    """Writes <name>-block.stl, <name>-insert-<k>.stl, <name>-inserts-plate.stl (CAD, needs build123d). Returns the file names."""
    from build123d import Pos
    os.makedirs(out_dir, exist_ok=True)
    files = {}
    solid = build_block(bp["cells"], bp["T"], d["width"], d["height"])
    files["block"] = f"{name}-block.stl"
    export_stl(solid, os.path.join(out_dir, files["block"]))
    if inserts:
        plate = None
        files["inserts"] = {}
        for item in plate_plan(d, bp["T"], bp["x_last"], bp["hc"], d["height"]):
            ins = build_insert(item["p"], d, bp["T"], bp["x_last"], bp["hc"], d["height"])
            fn = f"{name}-insert-{insert_file_suffix(item['k'])}.stl"
            export_stl(ins, os.path.join(out_dir, fn))
            files["inserts"][item["k"]] = fn
            piece = Pos(item["x0"], 0, 0) * ins
            plate = piece if plate is None else plate + piece
        files["plate"] = f"{name}-inserts-plate.stl"
        export_stl(plate, os.path.join(out_dir, files["plate"]))
    return files


def plan_json(d, bp):
    """JSON-able summary of a design and a block plan (fixtures)."""
    return dict(
        design=dict(L0=d["L0"], f0=d["f0"], f_mode=d["f_mode"], legs=d["plan"]["n"], bends=d["plan"]["bends"], y_in=d["plan"]["y_in"],
                    p_nom=d["plan"]["p_nom"], L_geom=d["plan"]["L_geom"], depth=d["depth"], p_limit=d["p_limit"],
                    ns_down=d["ns_down"], ns_up=d["ns_up"], warn=d["warn"], bend_dL=d["bend_dL"],
                    Lk={str(k): v for k, v in sorted(d["Lk"].items())}, pk={str(k): v for k, v in sorted(d["pk"].items())}),
        block=dict(size_mm=bp["size_mm"], grid=[bp["ny"], bp["nz"]], k=bp["k"], L_axis=bp["L_axis"], L_acoustic=bp["L_acoustic"],
                   volume_mm3=bp["volume_mm3"], x_last=bp["x_last"], hc=bp["hc"],
                   cells=[dict(index=c["index"], index_y=c["index_y"], index_z=c["index_z"], origin=c["origin"], mirror_y=c["mirror_y"],
                               air=c["air"], walls=c["walls"], voids=c["voids"], axis_points_xy=c["axis_points_xy"], window=c["window"])
                          for c in bp["cells"]]),
    )
