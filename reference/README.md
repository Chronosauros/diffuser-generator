# Python reference implementation

The source of truth for the maths and the tile geometry. The generator in `index.html` is a port of it, and
`npm test` checks the port against STL files made here (`scripts/fixtures/`).

- `qrd.py` - sequences, well depths, design bands and the sizing rules (pure functions, no geometry).
- `diffuser.py` - tile geometry with manifold-3d: wells, fins, closed air chambers, chamfers, the sunk tile number.
- `layout.py` - the far-field model and the optimiser that picked the wall orders (see `docs/how-it-works.md`).
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
```
