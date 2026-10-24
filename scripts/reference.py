"""Reference side of the equality test: runs the Python canon of the tile geometry on the fixture cases.

    python scripts/reference.py [<design-folder>]
    DIFFUSER_DESIGN=<design-folder> python scripts/reference.py

The Python canon (qrd.py, diffuser.py and out/layouts.json) is in reference/ in this repository, which is the default
<design-folder>. This script only regenerates scripts/fixtures/ (binary STL, JSON sidecars and expected.json) from it;
it is not needed for `npm test`, which uses the committed fixtures. <design-folder> is the folder that holds qrd.py and
diffuser.py, with out/layouts.json inside. It is imported, never edited, and nothing is written there.

Reads scripts/fixtures/cases.json, injects each params dict into qrd.DEEP / qrd.BLOCKS under a temporary key, runs
diffuser.deep_tile / diffuser.block_tile and writes binary STL + JSON sidecar (diffuser.save, its output folder is
pointed at fixtures/ for this run) for every case. Also writes fixtures/expected.json with the Python values the JS
test compares against: presets, bands, dims, file stems, depths, block heights, sequences. Needs Python with
manifold3d and trimesh.
"""
from __future__ import annotations

import json
import sys
import os
from pathlib import Path

sys.dont_write_bytecode = True                       # keep the design folder free of __pycache__

HERE = Path(__file__).resolve().parent
FIX = HERE / "fixtures"
_design_arg = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("DIFFUSER_DESIGN", "")
DESIGN = Path(_design_arg).resolve() if _design_arg else HERE.parent / "reference"
if (not (DESIGN / "qrd.py").is_file() or not (DESIGN / "diffuser.py").is_file()
        or not (DESIGN / "out" / "layouts.json").is_file()):
    sys.exit("reference.py needs the Python canon of the tile geometry (a design folder with qrd.py, diffuser.py and "
             "out/layouts.json, by default reference/ in this repository).\n"
             "Usage: python scripts/reference.py [<design-folder>]   or set DIFFUSER_DESIGN. "
             "`npm test` does not need it.")
sys.path.insert(0, str(DESIGN))

import diffuser  # noqa: E402
import qrd  # noqa: E402

diffuser.OUT = FIX                                   # save() writes to the module global OUT; keep the design folder untouched
LAYOUTS = json.loads((DESIGN / "out" / "layouts.json").read_text(encoding="utf-8"))


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    return x


def band_pack(series: str, size: str) -> dict:
    if series == "deep":
        return dict(band=qrd.deep_band(size), dims=list(qrd.deep_dims(size)))
    return dict(band=qrd.block_band(size), dims=list(qrd.block_dims(size)))


def run_case(case: dict) -> dict:
    cid, series, tile = case["id"], case["series"], case["tile"]
    key = "_case_" + cid.replace("-", "_")
    table = qrd.DEEP if series == "deep" else qrd.BLOCKS
    table[key] = dict(case["params"])                # temporary key, removed in finally
    try:
        out = band_pack(series, key)
        out["name"] = qrd.name(series, key, tile)
        if series == "deep":
            solid, info = diffuser.deep_tile(key, tile)
            out["depths"] = qrd.deep_depths(key, tile)
        else:
            k, l = case["shift"]
            solid, info = diffuser.block_tile(key, tile, k, l)
            out["heights"] = qrd.block_heights(key, k, l)
            out["shift"] = [k, l]
        saved = diffuser.save(solid, info, cid)
        out.update(file=saved["file"], volume_cm3=saved["volume_cm3"], triangles=saved["triangles"])
        solid_status = str(solid.status())
        out["manifold_status"] = solid_status
    finally:
        del table[key]
    return out


def main():
    cases = json.loads((FIX / "cases.json").read_text(encoding="utf-8"))
    expected = {
        "presets": {"deep": qrd.DEEP, "blocks": qrd.BLOCKS},
        "presets_band": {
            "deep": {sz: band_pack("deep", sz) for sz in qrd.DEEP},
            "blocks": {sz: band_pack("blocks", sz) for sz in qrd.BLOCKS},
        },
        # file stems of the 29 stage 1 tiles, in tile order
        "names": {
            "deep": {sz: [qrd.name("deep", sz, t) for t in range(1, qrd.DEEP_TILES + 1)] for sz in qrd.DEEP},
            "blocks": {sz: [qrd.name("blocks", sz, t) for t in range(1, len(LAYOUTS["blocks"][sz]["shifts"]) + 1)]
                       for sz in qrd.BLOCKS},
        },
        "deep_shift": {str(t): qrd.deep_shift(t) for t in range(1, qrd.DEEP_TILES + 1)},
        "deep_depths": {sz: {str(t): qrd.deep_depths(sz, t) for t in range(1, qrd.DEEP_TILES + 1)} for sz in qrd.DEEP},
        "residues": {str(n): qrd.residues(n) for n in (5, 7, 11, 13)},
        "shifted": {f"{n}:{k}": qrd.shifted(n, k) for n in (7, 11, 13) for k in range(n)},
        "phase_max": {f"{nx}x{ny}": qrd.phase_max(nx, ny) for nx, ny in ((13, 13), (13, 11), (7, 7))},
        "phase_fraction": [dict(Nx=nx, Ny=ny, n=n, m=m, k=k, l=l, v=qrd.phase_fraction(nx, ny, n, m, k, l))
                           for nx, ny in ((13, 13), (13, 11)) for (n, m, k, l) in
                           ((0, 0, 0, 0), (3, 5, 0, 0), (12, 10, 4, 9), (7, 2, 11, 1), (1, 1, 5, 5))],
        "block_heights": {sz: [qrd.block_heights(sz, k, l) for k, l in LAYOUTS["blocks"][sz]["shifts"]]
                          for sz in qrd.BLOCKS},
        "cases": {},
    }
    for case in cases:
        expected["cases"][case["id"]] = run_case(case)
        print(case["id"], expected["cases"][case["id"]]["file"], expected["cases"][case["id"]]["volume_cm3"], "cm3")
    (FIX / "expected.json").write_text(json.dumps(jsonable(expected), indent=1), encoding="utf-8")
    print("wrote", FIX / "expected.json")


if __name__ == "__main__":
    main()
