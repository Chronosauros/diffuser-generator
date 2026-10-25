"""Reference side of the bass traps tests: runs the Python twin (reference/bass_room.py, reference/bass_cell.py) and writes
scripts/fixtures/bass/ (JSON of the physics and of the geometry plans, binary STL of blocks and inserts).

    python scripts/bass_reference.py physics     numpy only; interpreter e.g. the research .venv-bem
    python scripts/bass_reference.py geometry    numpy + build123d (STL); JSON plans are written without build123d
    python scripts/bass_reference.py all         both (needs build123d for the STL files)
    python scripts/bass_reference.py geometry --research <research cell/stl-v2 folder>
                                                 also prints the comparison of the regenerated STL with the research STL

Not needed for `npm test`, which uses the committed fixtures. Numbers are written with 10 significant digits. The physics
cases: 5 rooms x alpha 0.05 / 0.10 / 0.20 x every mode up to 120 Hz x 7 anchors x N 1, 2, 4, 8, 16. The geometry cases: modes
35 - 120 Hz, block layouts 1x1, 1x2, 2x1, 2x2 with 3 and 4 cells, inserts and the inserts plate.
"""
from __future__ import annotations

import json
import math
import os
import struct
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "reference"))
OUT = HERE / "fixtures" / "bass"

import numpy as np  # noqa: E402
import bass_room as br  # noqa: E402

ROOMS = [
    ("example", (4.0, 3.2, 2.6)),
    ("small", (3.0, 2.5, 2.4)),
    ("large", (6.0, 4.5, 2.8)),
    ("long", (7.0, 3.0, 2.5)),
    ("near-cube", (3.0, 2.9, 2.7)),
]
ALPHAS = (0.05, 0.10, 0.20)
FMAX = 120.0
GEOM_FREQS = (35.0, 42.9, 53.6, 68.7, 77.3, 85.8, 95.2, 107.2, 120.0)
# (name, f mode, ny, nz, k, write STL)
GEOM_BLOCKS = [
    ("f042.9-1x1-k1", 42.9, 1, 1, 1, True),
    ("f068.7-2x1-k2", 68.7, 2, 1, 2, True),
    ("f095.2-1x2-k2", 95.2, 1, 2, 2, True),
    ("f095.2-2x2-k3", 95.2, 2, 2, 3, False),
    ("f095.2-2x2-k4", 95.2, 2, 2, 4, True),
    ("f107.2-2x2-k3", 107.2, 2, 2, 3, True),
]
STL_INSERTS = {"f095.2": (-2, 0, 3)}      # insert k written as its own STL; the plate is always written for these designs
STL_PLATES = ("f095.2",)


def rnd(x):
    """10 significant digits, recursively; numpy -> python."""
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return x if not math.isfinite(x) or x == 0 else float(f"{x:.10g}")
    if isinstance(x, (int, np.integer)) and not isinstance(x, bool):
        return int(x)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, dict):
        return {str(k): rnd(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [rnd(v) for v in x]
    return x


def dump(path: Path, obj, indent=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rnd(obj), indent=indent, separators=(",", ":") if indent is None else None, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", path.relative_to(HERE.parent), f"({path.stat().st_size / 1024:.0f} KiB)")


# ------------------------------------------------------------------------------------------------------ physics
def pack_mode(ev):
    """Compact record of br.evaluate_mode: arrays [anchor][N] in the order of ANCHORS and N_LIST."""
    an = list(br.ANCHORS)
    ns = [str(n) for n in br.N_LIST]

    def grid(key):
        return [[ev["anchors"][a]["by_n"][n].get(key) for n in ns] for a in an]
    return dict(
        n=ev["mode"]["n"], f=ev["mode"]["f"], kind=ev["mode"]["kind"], Lam=ev["mode"]["Lam"], delta_base=ev["delta_base"],
        t60_before=ev["t60_before"], L0=ev["L0"], f_cell=ev["f_cell"], gamma_cell=ev["gamma_cell"], no_effect=ev["no_effect"],
        overlap=ev["overlap"], uncertain=ev["uncertain"], psi2=[ev["anchors"][a]["psi2"] for a in an],
        nodal=[ev["anchors"][a]["nodal"] for a in an], delta_eff=grid("delta_eff"), peak_drop_db=grid("peak_drop_db"),
        best=[ev["best"][n] for n in ns])


def primitives():
    """Single numbers of the building blocks, for unit-level comparison."""
    L = (4.0, 3.2, 2.6)
    out = {}
    out["dl_slit_array"] = {"6": br.dl_slit_array(6e-3, 0.12), "8": br.dl_slit_array(8e-3, 0.12)}
    out["solve_L"] = {str(f): br.solve_L(f, 6.0) for f in (35.0, 42.9 * 1.005, 53.6 * 1.005, 68.7 * 1.005, 95.2 * 1.005, 120.0 * 1.005)}
    out["first_res"] = {str(l): br.first_res(6.0, l) for l in (500.0, 868.295, 1908.79, 2400.0)}
    zs = []
    for w, Lc, f, im in ((6.0, 1908.79, 42.9, 0.0), (6.0, 1908.79, 43.1, 5.0), (6.0, 868.295, 95.7, 12.0), (6.0, 1200.0, 70.0, -3.0), (8.0, 1000.0, 80.0, 8.0)):
        om = 2 * np.pi * f + 1j * im
        z = complex(br.zs_qwr(om, w, Lc))
        zs.append(dict(w=w, L=Lc, f=f, om_im=im, zs=[z.real, z.imag]))
    out["zs_qwr"] = zs
    cell = br.Cell(43.1145, 6.0, L=1908.7908590295256)
    out["absorption"] = {str(f): float(cell.alpha_f(f)) for f in (40.0, 43.1145, 47.0)}
    out["loading_mass"] = {a: {str(n): br.loading_mass(n, a) for n in (1, 4, 16)} for a in br.ANCHORS}
    out["sabine_decay"] = {str(a): br.sabine_decay(L, a) for a in ALPHAS}
    clus = {}
    for a in br.ANCHORS:
        cl = br.cluster(L, a, 20)
        clus[a] = [dict(wall=c["wall"], center=c["center"]) for c in cl]
    out["cluster_example"] = clus
    modes = br.room_modes(L, FMAX)
    out["cells_psi_example"] = {a: br.cells_psi(modes[:4], L, br.cluster(L, a, 3)) for a in br.ANCHORS}
    out["anchor_text_example"] = {a: br.anchor_text(a, L) for a in br.ANCHORS}
    return out


def physics():
    rooms = []
    for rid, L in ROOMS:
        cases = []
        for al in ALPHAS:
            modes, rows = br.analyze_room(L, al, FMAX)
            recs = []
            for i, row in enumerate(rows):
                ev = br.evaluate_mode(L, modes, i, al)
                rec = pack_mode(ev)
                rec.update(spacing_Hz=row["spacing_Hz"])
                recs.append(rec)
            cases.append(dict(alpha=al, sabine_t60=br.t60(br.sabine_decay(L, al)), modes=recs))
            print(f"physics {rid} alpha {al}: {len(recs)} modes", flush=True)
        rooms.append(dict(id=rid, L=list(L), cases=cases))
    obj = dict(format="bass-physics-fixture/1", note="calculated, not measured; written by scripts/bass_reference.py from reference/bass_room.py",
               constants=dict(C0=br.C0, RHO0=br.RHO0, MU=br.MU, KTH=br.KTH, CP=br.CP, GAM=br.GAM, P0=br.P0),
               fmax=FMAX, anchors=list(br.ANCHORS), ns=list(br.N_LIST), peak=dict(points=br.PEAK_POINTS, span=list(br.PEAK_SPAN)),
               primitives=primitives(), rooms=rooms)
    dump(OUT / "physics.json", obj)


# ------------------------------------------------------------------------------------------------------ geometry
def stl_stats(path: Path):
    b = path.read_bytes()
    n = struct.unpack("<I", b[80:84])[0]
    a = np.frombuffer(b[84:], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]), count=n)
    v = a["v"].astype(float)
    vol = float(np.einsum("ij,ij->i", v[:, 0], np.cross(v[:, 1], v[:, 2])).sum() / 6)
    pts = v.reshape(-1, 3)
    return dict(triangles=int(n), volume_mm3=vol, bbox_min=pts.min(0).tolist(), bbox_max=pts.max(0).tolist())


def strip_walls(pj):
    for c in pj["block"]["cells"]:
        c.pop("walls", None)
    return pj


def geometry(with_stl: bool, research: Path | None):
    import bass_cell as bc
    designs, plans, stl = {}, [], {}
    for f in GEOM_FREQS:
        designs[f"{f:g}"] = bc.design(f)
    configs = [(1, 1, 1), (1, 2, 2), (2, 1, 2), (2, 2, 3), (2, 2, 4)]
    for i, f in enumerate(GEOM_FREQS):
        d = designs[f"{f:g}"]
        use = configs if f in (42.9, 95.2) else [configs[i % 5], configs[(i + 2) % 5]]
        for ny, nz, k in use:
            block = bc.plan_json(d, bc.plan_block(d, ny, nz, k))["block"]
            if (ny, nz, k) != (1, 1, 1):                 # walls only for the single cell (they are not used by the solid)
                for c in block["cells"]:
                    c.pop("walls", None)
            plans.append(dict(f_mode=f, ny=ny, nz=nz, k=k, block=block))
    plates = {}
    for f in GEOM_FREQS:
        d = designs[f"{f:g}"]
        bp = bc.plan_block(d, 1, 1, 1)
        plates[f"{f:g}"] = [dict(k=it["k"], x0=it["x0"], p=it["p"], flange=it["dims"]["flange"], body=it["dims"]["body"], body_origin=it["dims"]["body_origin"])
                            for it in bc.plate_plan(d, bp["T"], bp["x_last"], bp["hc"], d["height"])]
    obj = dict(format="bass-geometry-fixture/1", note="calculated, not measured; written by scripts/bass_reference.py from reference/bass_cell.py",
               designs={k: bc_design_json(v) for k, v in designs.items()}, plans=plans, plates=plates)
    if with_stl:
        OUT.mkdir(parents=True, exist_ok=True)
        for name, f, ny, nz, k, write in GEOM_BLOCKS:
            if not write:
                continue
            d = designs[f"{f:g}"] if f"{f:g}" in designs else bc.design(f)
            bp = bc.plan_block(d, ny, nz, k)
            files = bc.export_set(OUT, f"bass-{name}", d, bp, inserts=False)
            stl[files["block"]] = dict(kind="block", f_mode=f, ny=ny, nz=nz, k=k, plan_volume_mm3=bp["volume_mm3"],
                                       pockets=sum(len(c["voids"]) for c in bp["cells"]), **stl_stats(OUT / files["block"]))
        for tag in STL_PLATES:
            f = float(tag[1:])
            d = designs[f"{f:g}"]
            bp = bc.plan_block(d, 1, 2, 2)
            files = bc.export_set(OUT, f"bass-{tag}", d, bp, inserts=True)
            keep = {files["inserts"][kk] for kk in STL_INSERTS[tag]} | {files["plate"]}
            for kk, fn in files["inserts"].items():
                if fn in keep:
                    stl[fn] = dict(kind="insert", f_mode=f, k=kk, **stl_stats(OUT / fn))
                else:
                    (OUT / fn).unlink()
            os.unlink(OUT / files["block"])          # the block of this design is covered by the 1x2 block above
            stl[files["plate"]] = dict(kind="plate", f_mode=f, **stl_stats(OUT / files["plate"]))
        obj["stl"] = stl
        if research:
            print("comparison with the research STL (volume, bbox):")
            for name, rf in (("bass-f042.9-1x1-k1-block.stl", "qwr2-m100-f042.9-w6-block.stl"), ("bass-f095.2-1x2-k2-block.stl", "qwr2-m111-f095.2-w6-block.stl"),
                             ("bass-f095.2-inserts-plate.stl", "qwr2-m111-f095.2-w6-inserts-plate.stl")):
                a, b = stl_stats(OUT / name), stl_stats(research / rf)
                print(f"  {name}: {a['volume_mm3']:.3f} vs {b['volume_mm3']:.3f} mm3 (rel {abs(a['volume_mm3'] / b['volume_mm3'] - 1):.2e}), "
                      f"bbox max {np.round(a['bbox_max'], 3).tolist()} vs {np.round(b['bbox_max'], 3).tolist()}")
    dump(OUT / "geometry.json", obj)


def bc_design_json(d):
    return dict(L0=d["L0"], f0=d["f0"], f_mode=d["f_mode"], legs=d["plan"]["n"], bends=d["plan"]["bends"], y_in=d["plan"]["y_in"],
                p_nom=d["plan"]["p_nom"], L_geom=d["plan"]["L_geom"], depth=d["depth"], p_limit=d["p_limit"], ns_down=d["ns_down"], ns_up=d["ns_up"],
                warn=d["warn"], bend_dL=d["bend_dL"], Lk={str(k): v for k, v in sorted(d["Lk"].items())}, pk={str(k): v for k, v in sorted(d["pk"].items())})


def main(argv):
    what = argv[0] if argv else "all"
    research = None
    if "--research" in argv:
        research = Path(argv[argv.index("--research") + 1])
    if what in ("physics", "all"):
        physics()
    if what in ("geometry", "all"):
        try:
            import build123d  # noqa: F401
            stl = True
        except ImportError:
            print("build123d not found: the STL files are skipped (use the build123d interpreter for them)")
            stl = False
        geometry(stl, research)


if __name__ == "__main__":
    main(sys.argv[1:])
