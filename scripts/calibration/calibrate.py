"""Calibrates the print time / filament mass estimate of the Soundblooms diffuser tiles.

This script documents how scripts/calibration/estimate.json was fitted. It needs a private design/print project (the
Python canon, the stage 1 slices and a slicer wrapper), which is not part of this repository, so it cannot run from a
clean checkout. `npm test` does not use it; the test only reads data.json and estimate.json.

    python scripts/calibration/calibrate.py <project-folder>             fit from data.json when it exists, else slice first
    python scripts/calibration/calibrate.py <project-folder> --reslice   rebuild every row: stage 1 slices + new slices, then fit
    (or set DIFFUSER_PROJECT=<project-folder>; DIFFUSER_CALIB=<folder> sets where slices and fit-report.json are written)

<project-folder> holds design/ (qrd.py, diffuser.py, out/, out-v3/), print/ (slices of the v3 set), tools/slice.py and
archive/stage1-full-bed/print/.

Data (data.json, one row per tile): params, features (CONTRACT.md: area_cm2, volume_cm3, height_mm, footprint_cm2 of the
tile in print orientation), sliced minutes and grams. Rows:
  - the 29 stage 1 tiles: features of the stage 1 STL files + the stage 1 slice results (read only);
  - NEW tiles at the derive-rule v2 sizes (beds A1 mini 180x180x180, MK4 250x210x220, P1S 256x256x256; default depth
    plus other depths): STL from the Python canon (params from the rule v2 derive, the same way as scripts/reference.py),
    sliced one at a time with the same presets as the slicer wrapper (P1S 0.4, Textured PEI, Bambu PLA Basic,
    "Diffuser deep 0.28" / "Diffuser blocks 0.20", no supports, chamber fan 0);
  - rule v3 deep tiles (walls grow with the depth, CONTRACT.md): the 21 slices of the v3 set (design/out-v3 + print/
    of the project, read only) and 4 new slices at other depths (see PLAN_V3). A data.json from before v3 gets these
    rows added on the next run; nothing else is sliced again.

Deep rows of rule v1 / v2 (stage 1, v2 slices) carry "in_fit": false: kept as a record, not fitted (the generator
makes rule v3 tiles only). The deep fit uses the 25 v3 rows, held-out groups = the three beds.

Fit per series: weighted least squares (weights 1/y, i.e. relative error), non-negative feature coefficients, feature
subset chosen by the error on HELD-OUT groups (leave one group out: stage 1 S, stage 1 L, v2 A1 mini, v2 MK4, v2 P1S, v3 A1 mini, v3 MK4, v3 P1S):
the simplest subset within 0.15 percentage points of the best held-out mean error. Writes estimate.json (ESTIMATE shape of
CONTRACT.md) and the fit report to the calibration output folder.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import math
import os
import re
import subprocess
import sys
import time
import zipfile
from pathlib import Path

sys.dont_write_bytecode = True                       # no __pycache__ in the project folders

HERE = Path(__file__).resolve().parent
_proj = next((a for a in sys.argv[1:] if not a.startswith("--")), os.environ.get("DIFFUSER_PROJECT", ""))
DIFF = Path(_proj).resolve() if _proj else None
if DIFF is None or not (DIFF / "design" / "qrd.py").is_file():
    sys.exit("calibrate.py needs the private design/print project (design/qrd.py, tools/slice.py, "
             "archive/stage1-full-bed/print/), which is not part of this repository.\n"
             "Usage: python scripts/calibration/calibrate.py <project-folder> [--reslice]   or set DIFFUSER_PROJECT. "
             "`npm test` does not need it.")
DESIGN, TOOLS, PRINT = DIFF / "design", DIFF / "tools", DIFF / "archive" / "stage1-full-bed" / "print"
CALIB = Path(os.environ.get("DIFFUSER_CALIB") or (HERE / "out"))   # slices and fit-report.json (out/ is git-ignored)
DATA, ESTIMATE = HERE / "data.json", HERE / "estimate.json"
sys.path.insert(0, str(DESIGN))

import numpy as np  # noqa: E402
import trimesh  # noqa: E402
import diffuser  # noqa: E402
import qrd  # noqa: E402

diffuser.OUT = CALIB / "stl"                         # save() writes to the module global OUT; keep the design folder untouched
LAYOUTS = json.loads((DESIGN / "out" / "layouts.json").read_text(encoding="utf-8"))

M = 25.0                                             # margin on every bed axis (CONTRACT.md, rule v2)
BEDS = {                                             # key: (label, x, y, z)
    "a1": ("Bambu Lab A1 mini", 180, 180, 180),
    "mk4": ("Prusa MK4 / Core One", 250, 210, 220),
    "p1s": ("Bambu Lab P1S / X1", 256, 256, 256),
}
# New slices: (series, bed, depth mm or None = default, tile). 9 deep + 6 blocks.
PLAN = [
    ("deep", "a1", None, 1), ("deep", "a1", 60, 4), ("deep", "a1", 120, 6),
    ("deep", "mk4", None, 1), ("deep", "mk4", 80, 4), ("deep", "mk4", 170, 6),
    ("deep", "p1s", None, 1), ("deep", "p1s", 80, 4), ("deep", "p1s", 200, 6),
    ("blocks", "a1", None, 1), ("blocks", "a1", 30, 2),
    ("blocks", "mk4", None, 1), ("blocks", "mk4", 40, 2),
    ("blocks", "p1s", None, 1), ("blocks", "p1s", 70, 2),
]
FEATS = ("area_cm2", "volume_cm3", "height_mm", "footprint_cm2")
ALL = FEATS + ("const",)


# ---------------------------------------------------------------- derive rule v2 (CONTRACT.md)
def ceil01(v: float) -> float:
    return math.ceil(v * 10 - 1e-9) / 10             # round up to 0.1 mm (CONTRACT.md)


def derive(series: str, bed: tuple, depth=None, rule: int = 3) -> dict:
    """Rule v2 (deep fin 1.2 / edge 0.8) for the v2 rows, rule v3 (fin / edge grow with the depth) by default."""
    x, y, z = bed
    label = f"custom {x} x {y} x {z}"
    if series == "deep":
        pitch = float((x - M) // 7)
        length = min(y - M, z - M)
        back = 2.0
        hi = y - M - back
        d = min(4 * pitch, hi) if depth is None else depth
        fin, edge = (1.2, 0.8) if rule == 2 else (max(1.6, ceil01(d / 82.5)), max(2.0, ceil01(d / 66)))
        return dict(N=7, pitch=pitch, length=float(length), dmax=float(d), fin=fin, edge=edge, floor=1.6, back=back,
                    cap=1.2, chamfer=3.0 if pitch < 30 else 4.0, bed=label)
    W = x - M
    H = min(y - M, W)
    Ny = 13 if H / W >= 12 / 13 else 11
    if H / W < 10 / 13:
        W = round(H * 13 / 10)
    base = 3.0
    hi = min(z - M - base, 4 * min(W / 13, H / Ny))
    h = min(48 if Ny == 13 else 60, int(np.floor(hi))) if depth is None else depth
    return dict(Nx=13, Ny=Ny, W=float(W), H=float(H), hmax=float(h), base=base, bed=label)


def build(series: str, params: dict, tile: int, stem_hint: str):
    """Tile via the Python canon with params injected under a temporary key (scripts/reference.py does the same)."""
    key = "_calib_" + stem_hint
    table = qrd.DEEP if series == "deep" else qrd.BLOCKS
    table[key] = dict(params)
    try:
        if series == "deep":
            solid, info = diffuser.deep_tile(key, tile)
        else:
            k, l = LAYOUTS["blocks"]["S" if params["Ny"] == 13 else "L"]["shifts"][tile - 1]
            solid, info = diffuser.block_tile(key, tile, k, l)
        stem = qrd.name(series, key, tile)
    finally:
        del table[key]
    return solid, info, stem


def features_of_mesh(tm: trimesh.Trimesh) -> dict:
    """CONTRACT.md `features`: surface area of the mesh (inner chamber walls included), volume, bbox Z, bbox X x Y."""
    lo, hi = tm.bounds
    e = hi - lo
    return dict(area_cm2=float(tm.area) / 100, volume_cm3=float(tm.volume) / 1000, height_mm=float(e[2]),
                footprint_cm2=float(e[0] * e[1]) / 100)


def mesh_of(solid) -> trimesh.Trimesh:
    m = solid.to_mesh()
    return trimesh.Trimesh(vertices=np.asarray(m.vert_properties)[:, :3], faces=np.asarray(m.tri_verts), process=False)


# ---------------------------------------------------------------- slicing (same presets and chamber fan fix as the slicer wrapper)
def load_slicer():
    spec = importlib.util.spec_from_file_location("diffuser_slice_tool", TOOLS / "slice.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)                     # patches slice_plates: profiles, filament, chamber fan 0
    mod.sp.TMP = CALIB / "slice"
    return mod.sp


def bambu_pids() -> set[str]:
    out = subprocess.run(["tasklist", "/FI", "IMAGENAME eq bambu-studio.exe", "/FO", "CSV", "/NH"],
                         capture_output=True, text=True).stdout
    return {ln.split('","')[1] for ln in out.splitlines() if ln.startswith('"bambu')}


def chamber_fan_report(sp, rec: dict, stem: str) -> dict:
    """Chamber fan proof: resolved filament preset and the G-code of the sliced 3MF have no M106 P3 S>0."""
    fil = json.loads((sp.TMP / stem / "filament-full.json").read_text(encoding="utf-8"))
    preset = [int(v) for key in ("filament_start_gcode", "filament_end_gcode") for x in fil.get(key, [])
              for v in re.findall(r"M106 P3 S(\d+)", x)]
    z = zipfile.ZipFile(rec["file"])
    g = z.read(next(n for n in z.namelist() if n.endswith(".gcode"))).decode("utf-8", "replace")
    code = [int(v) for v in re.findall(r"^M106 P3 S(\d+)", g, re.M)]
    return dict(preset_nonzero=sum(1 for v in preset if v), gcode_nonzero=sum(1 for v in code if v), gcode_m106_p3=len(code))


def slice_new(sp, stl: Path, series: str, stem: str, timeout_note: str = "") -> dict:
    t0 = time.time()
    rec = sp.slice_plate(dict(id=1, stem=stem, filament=series, parts=[(stem, 1)]), {stem: dict(path=str(stl))})
    if "error" in rec or "file" not in rec:
        raise RuntimeError(f"slice failed for {stem}: {rec}")
    fan = chamber_fan_report(sp, rec, stem)
    if fan["preset_nonzero"] or fan["gcode_nonzero"]:
        raise RuntimeError(f"chamber fan not 0 in {stem}: {fan}")
    rec["wall_seconds"] = round(time.time() - t0, 1)
    rec["fan"] = fan
    return rec


# ---------------------------------------------------------------- rows
def stage1_rows() -> list[dict]:
    rows = []
    for r in json.loads((DESIGN / "out" / "manifest.json").read_text(encoding="utf-8")):
        stl = DESIGN / "out" / r["file"]
        info = json.loads(stl.with_suffix(".json").read_text(encoding="utf-8"))
        sl = json.loads((PRINT / (stl.stem + ".slice.json")).read_text(encoding="utf-8"))
        tm = trimesh.load(stl, process=False)
        rows.append(dict(id=stl.stem, source="stage1", group=f"stage1-{r['size']}", series=r["series"], tile=r["tile"],
                         depth=info["params"]["dmax" if r["series"] == "deep" else "hmax"], params=info["params"],
                         features=features_of_mesh(tm), minutes=round(sl["seconds"] / 60, 3), grams=sl["grams"],
                         slice_seconds=sl.get("slice_seconds")))
    return rows


def new_rows(reslice_log: list) -> list[dict]:
    sp = load_slicer()
    CALIB.mkdir(parents=True, exist_ok=True)
    before = bambu_pids()
    rows = []
    try:
        for series, bed, depth, tile in PLAN:
            label, x, y, z = BEDS[bed]
            params = derive(series, (x, y, z), depth, rule=2)
            d = params["dmax" if series == "deep" else "hmax"]
            solid, info, stem = build(series, params, tile, f"{series}_{bed}_{d:g}")
            saved = diffuser.save(solid, info, stem)
            stl = diffuser.OUT / saved["file"]
            tm = trimesh.load(stl, process=False)
            feats = features_of_mesh(tm)
            rec = slice_new(sp, stl, series, stem)
            print(f"{stem}: {rec['hm']} {rec['grams']} g (slice {rec['wall_seconds']} s)", flush=True)
            reslice_log.append(rec["wall_seconds"])
            rows.append(dict(id=stem, source="new", group=f"v2-{bed}", series=series, tile=tile, depth=d, params=params,
                             features=feats, minutes=round(rec["seconds"] / 60, 3), grams=rec["grams"],
                             slice_seconds=rec["slice_seconds"], chamber_fan=rec["fan"],
                             warnings=rec.get("warnings", [])))
            solid.delete() if hasattr(solid, "delete") else None
    finally:
        new = bambu_pids() - before                  # only processes this run started; the owner's windows stay
        for pid in new:
            subprocess.run(["taskkill", "/PID", pid, "/F", "/T"], capture_output=True)
    return rows


# ---------------------------------------------------------------- rule v3 rows (deep walls grow with the depth)
# 21 existing slices of the v3 set (design/out-v3 STL + print/*.slice.json of the project, read only: nothing is sliced
# for them) and 4 NEW deep slices at other depths (the maximum depth of the P1S, a shallow P1S tile, the maximum depth of
# the A1 mini and a shallow MK4 tile), because the v3 set only has the default depth of each bed and the wall thickness
# changes with the depth.
BED_OF_SIZE = {"A1mini": "a1", "MK4": "mk4", "P1S": "p1s"}
PLAN_V3 = [("p1s", 229, 5), ("p1s", 60, 4), ("a1", 153, 6), ("mk4", 90, 2)]


def v3_existing_rows() -> list[dict]:
    rows = []
    for sl_path in sorted((DIFF / "print").glob("soundblooms-diffuser-qrd-deep-*.slice.json")):
        sl = json.loads(sl_path.read_text(encoding="utf-8"))
        stl = DESIGN / "out-v3" / (sl["stem"] + ".stl")
        info = json.loads(stl.with_suffix(".json").read_text(encoding="utf-8"))
        rows.append(dict(id=sl["stem"], source="v3", group=f"v3-{BED_OF_SIZE[sl['size']]}", series="deep", tile=sl["tile"],
                         depth=info["params"]["dmax"], params=info["params"], features=features_of_mesh(trimesh.load(stl, process=False)),
                         minutes=round(sl["seconds"] / 60, 3), grams=sl["grams"], slice_seconds=sl.get("slice_seconds")))
    return rows


def v3_new_rows(reslice_log: list) -> list[dict]:
    sp = load_slicer()
    CALIB.mkdir(parents=True, exist_ok=True)
    before = bambu_pids()
    rows = []
    try:
        for bed, depth, tile in PLAN_V3:
            label, x, y, z = BEDS[bed]
            params = derive("deep", (x, y, z), float(depth), rule=3)
            solid, info, stem = build("deep", params, tile, f"deep_{bed}_{depth:g}")
            saved = diffuser.save(solid, info, stem)
            stl = diffuser.OUT / saved["file"]
            feats = features_of_mesh(trimesh.load(stl, process=False))
            rec = slice_new(sp, stl, "deep", stem)
            print(f"{stem}: {rec['hm']} {rec['grams']} g (fin {params['fin']}, edge {params['edge']}, slice {rec['wall_seconds']} s)", flush=True)
            reslice_log.append(rec["wall_seconds"])
            rows.append(dict(id=stem, source="v3-new", group=f"v3-{bed}", series="deep", tile=tile, depth=float(depth), params=params,
                             features=feats, minutes=round(rec["seconds"] / 60, 3), grams=rec["grams"],
                             slice_seconds=rec["slice_seconds"], chamber_fan=rec["fan"], warnings=rec.get("warnings", [])))
            solid.delete() if hasattr(solid, "delete") else None
    finally:
        for pid in bambu_pids() - before:
            subprocess.run(["taskkill", "/PID", pid, "/F", "/T"], capture_output=True)
    return rows


# ---------------------------------------------------------------- fit
def design_matrix(rows, names):
    return np.array([[1.0 if n == "const" else r["features"][n] for n in names] for r in rows])


def wls(rows, names, target):
    X, y = design_matrix(rows, names), np.array([r[target] for r in rows])
    w = 1.0 / y
    coef, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    return coef


def pct_err(rows, names, coef, target):
    X, y = design_matrix(rows, names), np.array([r[target] for r in rows])
    return np.abs(X @ coef - y) / y * 100


def logo(rows, names, target):
    """Leave-one-group-out: fit without a group, error on that group. Returns (all held-out errors %, per group)."""
    errs, per = [], {}
    for g in sorted({r["group"] for r in rows}):
        tr = [r for r in rows if r["group"] != g]
        te = [r for r in rows if r["group"] == g]
        e = pct_err(te, names, wls(tr, names, target), target)
        errs.extend(e)
        per[g] = dict(n=len(te), mean=round(float(e.mean()), 2), max=round(float(e.max()), 2))
    return np.array(errs), per


def choose(rows, target, gain=0.15):
    """Every subset of area_cm2, volume_cm3, height_mm, footprint_cm2, const with non-negative feature coefficients
    (const is free, a negative physical coefficient would break extrapolation). Chosen: the simplest model whose
    held-out (leave-one-group-out) mean error is within `gain` percentage points of the best one."""
    results = []
    for k in range(1, len(ALL) + 1):
        for names in itertools.combinations(ALL, k):
            names = list(names)
            coef = wls(rows, names, target)
            if any(c < 0 for n, c in zip(names, coef) if n != "const"):
                continue
            if (design_matrix(rows, names) @ coef).min() <= 0:
                continue
            e, per = logo(rows, names, target)
            results.append((float(e.mean()), float(e.max()), names, coef, per))
    results.sort(key=lambda t: t[0])
    best = results[0]
    pick = min((r for r in results if r[0] <= best[0] + gain), key=lambda r: (len(r[2]), r[0]))
    return pick, results


def mark_legacy(rows):
    """Deep tiles of rule v1 / v2 (fin 1.2, edge 0.8: stage 1 and the v2 slices) stay in data.json as a record but are not
    fitted: the generator makes rule v3 tiles only, and thin walls print with a different mass per volume."""
    for r in rows:
        if r["series"] == "deep" and not r["source"].startswith("v3"):
            r["in_fit"] = False
    return rows


def fit(rows):
    out, report = {}, {}
    for series in ("deep", "blocks"):
        rs = [r for r in rows if r["series"] == series and r.get("in_fit", True)]
        out[series], report[series] = {}, dict(n=len(rs))
        for target in ("grams", "minutes"):
            (mean, mx, names, coef, per), allm = choose(rs, target)
            train = pct_err(rs, names, coef, target)
            out[series][target] = {n: float(c) for n, c in zip(names, coef)}
            report[series][target] = dict(
                formula=" + ".join(f"{c:.6g}*{n}" if n != "const" else f"{c:.6g}" for n, c in zip(names, coef)),
                heldout_mean_pct=round(mean, 2), heldout_max_pct=round(mx, 2), per_group=per,
                train_mean_pct=round(float(train.mean()), 2), train_max_pct=round(float(train.max()), 2),
                best_candidates=[dict(features=a[2], heldout_mean=round(a[0], 2), heldout_max=round(a[1], 2))
                                 for a in allm[:8]])
    return out, report


def estimate(est, series, feats):
    return {t: sum(c * (1.0 if n == "const" else feats[n]) for n, c in est[series][t].items()) for t in ("grams", "minutes")}


def sanity(est, rows):
    """A tile at 1.5x the P1S size (bed 384 x 384 x 384) against the P1S default tile: positive, larger values."""
    res = {}
    for series in ("deep", "blocks"):
        feats = {}
        for lab, bed in (("p1s", (256, 256, 256)), ("p1s_x1.5", (384, 384, 384))):
            params = derive(series, bed)
            solid, info, stem = build(series, params, 1, f"{series}_sanity_{lab}")
            feats[lab] = features_of_mesh(mesh_of(solid))
            feats[lab]["dims"] = info["dims_whd"]
        ev = {lab: estimate(est, series, f) for lab, f in feats.items()}
        res[series] = dict(features=feats, estimate=ev,
                           ok=bool(all(v > 0 for e in ev.values() for v in e.values()) and
                                   all(ev["p1s_x1.5"][t] > ev["p1s"][t] for t in ("grams", "minutes"))))
    return res


def main(argv):
    reslice = "--reslice" in argv
    log = []
    if DATA.exists() and not reslice:
        rows = json.loads(DATA.read_text(encoding="utf-8"))
        print(f"{DATA.name} exists: {len(rows)} rows, no slicing (--reslice rebuilds it)")
        if not any(r["source"].startswith("v3") for r in rows):          # data.json from before rule v3: add the v3 rows
            t0 = time.time()
            rows = rows + v3_existing_rows() + v3_new_rows(log)
            DATA.write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
            print(f"added the v3 rows ({len(log)} sliced in {time.time() - t0:.0f} s wall), wrote {DATA}")
    else:
        t0 = time.time()
        rows = stage1_rows() + new_rows(log) + v3_existing_rows() + v3_new_rows(log)
        DATA.write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"sliced {len(log)} new tiles in {time.time() - t0:.0f} s wall, wrote {DATA}")
    mark_legacy(rows)
    DATA.write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
    est, report = fit(rows)
    ESTIMATE.write_text(json.dumps(est, indent=1), encoding="utf-8")
    report["sanity_1_5x"] = sanity(est, rows)
    CALIB.mkdir(parents=True, exist_ok=True)
    (CALIB / "fit-report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
    for series in ("deep", "blocks"):
        for t in ("grams", "minutes"):
            r = report[series][t]
            print(f"{series:6s} {t:7s} held-out max {r['heldout_max_pct']:5.2f} % mean {r['heldout_mean_pct']:5.2f} %  "
                  f"(train max {r['train_max_pct']} mean {r['train_mean_pct']})  {r['formula']}")
            print("        per group:", r["per_group"])
    print(json.dumps(report["sanity_1_5x"], indent=1))
    print("wrote", ESTIMATE)


if __name__ == "__main__":
    main(sys.argv[1:])
