"""Fits the print time / filament estimate of the bass trap blocks and inserts plates (core `bass-geometry`).

    python scripts/calibration/bass-calibrate.py <research-folder>

<research-folder> is the private research package `.../research/2026-10-25-mode-traps/cell/` (read only), which holds
`stl/` and `stl-v2/` (the STL files) and `slice_results.json` / `slice_results_v2.json` (Bambu Studio CLI slices: P1S 0.4
nozzle, Textured PEI, Soundblooms PLA Basic, process "Diffuser blocks 0.20", 3 walls, 15 % infill, no supports, chamber fan 0).
`npm test` does not run this script; the tests read bass-data.json and bass-estimate.json only. Needs numpy only.

Rows: the 10 blocks and 10 inserts plates of the two research sets (v1: 5 modes, no bend correction, 7 inserts per plate;
v2: 5 modes, the current defaults, 11 inserts per plate). Features of each STL as in CONTRACT.md: area_cm2 (mesh surface),
volume_cm3, height_mm (bbox Z), footprint_cm2 (bbox X x Y). Fit per kind (block, plate): weighted least squares (weights 1/y,
relative error), non-negative feature coefficients, the subset chosen by the error on HELD-OUT groups (leave one mode out: both
versions of a mode leave together), the simplest subset within 0.15 percentage points of the best held-out mean error - the
same rule as calibrate.py of the diffuser tiles. Writes bass-data.json, bass-estimate.json and bass-fit-report.json.
Ten rows per kind: this is a rough estimate (see the held-out errors in the report), shown as approximate.
"""
from __future__ import annotations

import itertools
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FEATS = ("area_cm2", "volume_cm3", "height_mm", "footprint_cm2")
ALL = FEATS + ("const",)
# blocks use only features that follow from the rectangle plan (volume, height, footprint), so the packer needs no solid;
# the plate is built as a solid anyway (boxes), so it may use the surface area too
ALLOWED = {"block": ("volume_cm3", "height_mm", "footprint_cm2", "const"), "plate": ALL}


def stl_features(path: Path) -> dict:
    b = path.read_bytes()
    n = struct.unpack("<I", b[80:84])[0]
    a = np.frombuffer(b[84:], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]), count=n)
    v = a["v"].astype(float)
    vol = float(np.einsum("ij,ij->i", v[:, 0], np.cross(v[:, 1], v[:, 2])).sum() / 6)
    area = float(0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1).sum())
    pts = v.reshape(-1, 3)
    e = pts.max(0) - pts.min(0)
    return dict(area_cm2=area / 100, volume_cm3=vol / 1000, height_mm=float(e[2]), footprint_cm2=float(e[0] * e[1]) / 100)


def rows_of(cell: Path) -> list[dict]:
    rows = []
    for ver, folder, res in (("v1", "stl", "slice_results.json"), ("v2", "stl-v2", "slice_results_v2.json")):
        slices = json.loads((cell / res).read_text(encoding="utf-8"))
        for label, sl in slices.items():
            kind = "block" if label.endswith("-block") else "plate"
            f = float(re.search(r"f(\d+\.\d)", label).group(1))
            rows.append(dict(id=label, source=ver, group=f"f{f:.1f}", kind=kind, mode_hz=f, features=stl_features(cell / folder / (label + ".stl")),
                             minutes=round(sl["seconds"] / 60, 3), grams=sl["grams"]))
    return rows


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
    errs, per = [], {}
    for g in sorted({r["group"] for r in rows}):
        tr = [r for r in rows if r["group"] != g]
        te = [r for r in rows if r["group"] == g]
        e = pct_err(te, names, wls(tr, names, target), target)
        errs.extend(e)
        per[g] = dict(n=len(te), mean=round(float(e.mean()), 2), max=round(float(e.max()), 2))
    return np.array(errs), per


def choose(rows, target, allowed, gain=0.15):
    results = []
    for k in range(1, 4):                           # at most 3 terms for 10 rows
        for names in itertools.combinations(allowed, k):
            names = list(names)
            coef = wls(rows, names, target)
            if any(c < 0 for c in coef):            # every coefficient non-negative, the constant too (no negative mass for a small plate)
                continue
            if (design_matrix(rows, names) @ coef).min() <= 0:
                continue
            e, per = logo(rows, names, target)
            results.append((float(e.mean()), float(e.max()), names, coef, per))
    results.sort(key=lambda t: t[0])
    best = results[0]
    pick = min((r for r in results if r[0] <= best[0] + gain), key=lambda r: (len(r[2]), r[0]))
    return pick, results


def main(argv):
    if not argv or not (Path(argv[0]) / "slice_results_v2.json").is_file():
        sys.exit("usage: python scripts/calibration/bass-calibrate.py <research cell folder with stl/, stl-v2/, slice_results*.json>")
    rows = rows_of(Path(argv[0]).resolve())
    (HERE / "bass-data.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
    est, report = {}, {}
    for kind in ("block", "plate"):
        rs = [r for r in rows if r["kind"] == kind]
        est[kind], report[kind] = {}, dict(n=len(rs))
        for target in ("grams", "minutes"):
            (mean, mx, names, coef, per), allm = choose(rs, target, ALLOWED[kind])
            train = pct_err(rs, names, coef, target)
            est[kind][target] = {n: float(c) for n, c in zip(names, coef)}
            report[kind][target] = dict(
                formula=" + ".join(f"{c:.6g}*{n}" if n != "const" else f"{c:.6g}" for n, c in zip(names, coef)),
                heldout_mean_pct=round(mean, 2), heldout_max_pct=round(mx, 2), per_group=per,
                train_mean_pct=round(float(train.mean()), 2), train_max_pct=round(float(train.max()), 2),
                best_candidates=[dict(features=a[2], heldout_mean=round(a[0], 2), heldout_max=round(a[1], 2)) for a in allm[:8]])
            r = report[kind][target]
            print(f"{kind:5s} {target:7s} held-out mean {r['heldout_mean_pct']:5.2f} % max {r['heldout_max_pct']:5.2f} %  "
                  f"(train mean {r['train_mean_pct']} max {r['train_max_pct']})  {r['formula']}")
    (HERE / "bass-estimate.json").write_text(json.dumps(est, indent=1), encoding="utf-8")
    (HERE / "bass-fit-report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
    print("wrote bass-data.json, bass-estimate.json, bass-fit-report.json in", HERE)


if __name__ == "__main__":
    main(sys.argv[1:])
