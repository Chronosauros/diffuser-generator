# Generator contract

The browser generator is `index.html`. The geometry it ports is the Python canon of stage 1 (`qrd.py`, `diffuser.py`),
kept in a private design project; its output for the 29 stage 1 tiles is committed as test fixtures in
`scripts/fixtures/stage1/` (STL, JSON sidecars and a trimmed `layouts.json` with the wall layouts).

## Files

- `index.html` - the whole app, one file, opens from disk (`file://`) and from any static host.
- `package.json` - `npm test` = `node --test scripts/*.test.mjs`; devDependency `manifold-3d` pinned to **3.5.3**
  (installed by `npm install`). The browser loads the same version: `https://cdn.jsdelivr.net/npm/manifold-3d@3.5.3/+esm`.
  Python canon runs manifold3d 3.5.2.
- `scripts/` - the equality test, `reference.py` (regenerates `fixtures/` from the private Python canon, not needed for
  `npm test`), `fixtures/` and `calibration/` (data and fit of the print time / filament estimate).

## Core section of index.html

- A classic script element `<script id="diffuser-core">...</script>` (not `type="module"`, no `import`, no DOM, no
  `window`/`document`, no fetch). It defines one global function `createDiffuser(M)`, where `M` is the initialised
  manifold module (the object returned by `await Module()` from manifold-3d, with `M.Manifold`, `M.CrossSection`);
  the caller runs `M.setup()` before `createDiffuser(M)` (idempotent).
- The rest of the app (UI, three.js preview, STL/ZIP export) lives in other script elements and calls
  `createDiffuser(await Module())`.
- The test loads it like this: read index.html, take the text of `<script id="diffuser-core">`, run it in `node:vm`
  with `Math`, `console` in the context, then call `createDiffuser(M)` with a module from `node_modules/manifold-3d`.

`createDiffuser(M)` returns:

| member | meaning |
|---|---|
| `PRESETS` | `{deep: {S, L}, blocks: {S, L}}` - field-for-field copies of `qrd.DEEP` and `qrd.BLOCKS` (same keys, same values, `bed` string included) |
| `LAYOUTS` | `{deep: {S: {room_left_to_right}, L: {...}}, blocks: {S: {shifts, grids}, L: {...}}}` - copied from `layouts.json` (only these keys; no `model`, `score`, `evaluation`) |
| `residues(N)`, `shifted(N,k)`, `deepShift(tile)`, `deepDepths(params, tile)` | ports of qrd.py |
| `phaseFraction(Nx,Ny,n,m,k,l)`, `phaseMax(Nx,Ny)`, `blockHeights(params,k,l)` | ports of qrd.py (`[m][n]`) |
| `deepBand(params)`, `blockBand(params)` | ports of `deep_band` / `block_band`, same keys and units |
| `deepDims(params)`, `blockDims(params)` | ports, arrays `[a, b, c]` in mm |
| `deepTile(params, tile)` | returns a `M.Manifold` equal to `diffuser.deep_tile` for those params (label included) |
| `blockTile(params, tile, k, l)` | returns a `M.Manifold` equal to `diffuser.block_tile` (label included) |
| `name(series, params, tile)` | file stem exactly like `qrd.name` (`soundblooms-diffuser-qrd-deep-175x175x102-tile1`) |
| `derive(series, bed, depth)` | params object for a bed `{x, y, z}` in mm and a depth in mm (rule below) |
| `limits(series, bed)` | `{depthMin, depthMax, depthDefault}` for the slider |

`params` always has the field names of the Python dicts: deep `N, pitch, length, dmax, fin, edge, floor, back, cap,
chamfer, bed`; blocks `Nx, Ny, W, H, hmax, base, bed`. The caller frees Manifold objects with `.delete()`.

## derive - rule v2, N is fixed

Tiles at the full bed size are unrealistic and tend to lift off the bed, so the margin is **M = 25 mm on every bed
axis**: every tile is at least 20 mm smaller than its stage 1 size. `PRESETS` stay the stage 1 dicts (the geometry canon
for the equality test); they are no longer reachable from the UI and no round trip to them is required. Beds in the UI:
Bambu Lab A1 mini 180 x 180 x 180, Prusa MK4 / Core One 250 x 210 x 220, Bambu Lab P1S / X1
256 x 256 x 256, Custom (x, y, z).

Deep (printed standing: width W on bed X, depth T on bed Y, length along the wells on bed Z):
`N = 7`; `pitch = floor((x - M) / 7)` whole mm; `fin 1.2, edge 0.8, floor 1.6, back 2.0, cap 1.2` fixed;
`chamfer = pitch < 30 ? 3 : 4`; `length = min(y - M, z - M)`; `dmax = depth`.
Limits: `dmax` from 40 to `y - M - back`; default = `min(4 * pitch, max)` (the "matched" tile:
unit depth = pitch, one tile = half the design wavelength).
Results: A1 mini pitch 22, tile 154 x 155, max 153, default 88; MK4 pitch 32, 224 x 185, max 183, default 128;
P1S pitch 33, 231 x 231, max 229, default 132.

Blocks (printed lying on the back: W on bed X, H on bed Y, height on bed Z):
`Nx = 13`; `W = x - M`; `H = min(y - M, W)`; `Ny = H / W >= 12/13 ? 13 : 11`; when `H / W < 10/13`, reduce `W` to
`Math.round(H * 13 / 10)` (JS rounding, half up); `base = 3.0`; `hmax = depth`. In wallBand, `wallHeightMm` = rows x
tile length (deep) or rows x H (blocks).
Limits: `hmax` from 20 to `min(z - M - base, 4 * min(W / Nx, H / Ny))`; default =
`min(Ny == 13 ? 48 : 60, floor(max))`.
Results: A1 mini 155 x 155 (13 x 13), max 47.69, default 47; MK4 225 x 185 (13 x 11), max 67.27, default 60;
P1S 231 x 231 (13 x 13), max 71.08, default 48.

## Set semantics (v2)

- A set is the selected wall from `wallOptions`. Set grams and minutes = the sum over every cell of `grid` (deep tiles
  differ in mass; repeated rows count every printed tile).
- ZIP = one STL per distinct tile number in the wall + `wall-plan.svg` of the selected wall + `README.txt` with the
  quantity to print per tile number, params, band and estimates.
- The wall plan section shows the selected wall only. The Tile selector (preview and single STL) offers only tile
  numbers present in the selected wall.

Tile sets: deep wall order = `LAYOUTS.deep.S` when `pitch < 30`, else `LAYOUTS.deep.L`. Blocks: shifts and grids from
`LAYOUTS.blocks.S` when `Ny = 13`, `LAYOUTS.blocks.L` when `Ny = 11`.

## Wall choice, wall band, estimates - new members of createDiffuser(M) (v2)

| member | meaning |
|---|---|
| `wallOptions(series, params)` | list of walls the user can pick, smallest first. Each `{id, cols, rows, tiles, grid}` where `grid` is rows (top to bottom) of tile numbers seen from the room, left to right, and `tiles` = cols x rows. Always first `{cols 1, rows 1}` = tile 1 alone ("one tile to try"). Deep: every wall width in the chosen `room_left_to_right` (S: 4, 7; L: 3, 5, 7), each with rows 1, 2 and 3 (more rows repeat the same row). Blocks: every grid in the chosen `grids` (S: 2x2, 2x3, 3x3; L: 2x2, 2x3). |
| `wallBand(series, params, cols, rows)` | `{low, high, wallWidthMm, wallHeightMm}`: `high = band.f_width`, `wallWidthMm = cols * tile width` (deep: W; blocks: W), `low = max(band.f0, 343 / (wallWidthMm / 1000))` in Hz. Below c / wall width a wall scatters little, so the low end is tied to wall width, never to f0 alone. |
| `features(manifold)` | `{area_cm2, volume_cm3, height_mm, footprint_cm2}` of a tile in print orientation (STL as built: Z up = print height; footprint = bbox X x Y). |
| `ESTIMATE` | coefficients from `scripts/calibration/estimate.json` (written by `scripts/calibration/calibrate.py`): `{deep: {grams: {...}, minutes: {...}}, blocks: {...}}`, each a map feature name (`area_cm2`, `volume_cm3`, `height_mm`, `footprint_cm2`, `const`) -> number. |
| `estimate(series, feats)` | `{grams, minutes}` = sum of coefficient x feature (+ const) per series, over the keys present in `ESTIMATE` only (a missing key counts as 0; the fit picked different features per series). Bambu P1S, 0.4 nozzle, Bambu PLA Basic, the stage 1 process profiles; shown in the UI as approximate. |

## Equality ("the same STL")

For each of the 29 stage 1 tiles (params = preset, shifts from layouts.json) and for every case in
`scripts/fixtures/cases.json` (random params, Python reference STL made by `scripts/reference.py`): the JS mesh and
the reference STL, both measured by the same JS code - heights by ray cast at the measurement points of the Python reference tests
(+-0.01 mm), bounding box (+-0.01 mm), volume (+-0.1 %), shells: one closed outer solid, deep also the same number of
inner chamber shells as the reference.
