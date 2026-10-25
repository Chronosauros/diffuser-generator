# Python reference implementation

The source of truth for the maths and the tile geometry. The generator in `index.html` is a port of it, and
`npm test` checks the port against STL files made here (`scripts/fixtures/`).

- `qrd.py` - sequences, well depths, design bands and the sizing rules (pure functions, no geometry).
- `diffuser.py` - tile geometry with manifold-3d: wells, fins, closed air chambers, chamfers, the sunk tile number.
- `layout.py` - the far-field model and the optimiser that picked the wall orders (see `docs/how-it-works.md`).
- `bass_room.py`, `bass_cell.py` - the Python twins of the bass trap cores of `index.html` (`bass-physics`, `bass-geometry`): room
  modes, the QWR cell impedance, the single-mode trap model, flags and peak drop in dB; the folded cell, blocks, inserts, the plate
  and the STL export (build123d, imported only inside the export function). Trimmed copies of the research of 25.10.2026, each with
  a header that names its source files. Calculated, not measured: `docs/bass-traps.md`. Numpy only, except the STL export.
- `out/layouts.json` - the frozen wall orders and block shifts of stage 1, which the current sizes copy. The model scores
  are left out on purpose: they are predictions for comparing layouts, not measurements.

The files were moved here from a private design project, so comments that say `design/` mean this folder.

## Running it

Python 3.11 or newer:

```
pip install -r reference/requirements.txt
python reference/qrd.py                    # print sizes, bands and depths
python reference/diffuser.py deep P1S 3    # build one tile, writes STL and JSON to reference/out-v3/
python reference/diffuser.py               # build the whole current set into reference/out-v3/
python reference/layout.py                 # evaluate the current wall orders with the model, plots in out-v3/model/
python reference/layout.py --stage1 DIR    # search the orders again (slow), into DIR
python scripts/reference.py                # regenerate the test fixtures from this folder
python scripts/bass_reference.py physics   # bass traps: scripts/fixtures/bass/physics.json (numpy only)
python scripts/bass_reference.py geometry  # bass traps: geometry.json and the STL files (needs build123d)
```

The bass trap twins need numpy; `bass_cell.export_set` (STL files of a block, the inserts and the plate) needs build123d and
imports it only when called. There is no manifold3d in the environments the fixtures were made in; the expected volumes come
from the rectangle plan and from build123d.
