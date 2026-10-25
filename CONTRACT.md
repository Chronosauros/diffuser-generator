# Generator contract

The browser generator is `index.html`. The geometry it ports is the Python canon of stage 1 (`qrd.py`, `diffuser.py`),
in `reference/`; its output for the 29 stage 1 tiles is committed as test fixtures in
`scripts/fixtures/stage1/` (STL, JSON sidecars and a trimmed `layouts.json` with the wall layouts).

## Files

- `index.html` - the whole app, one file, opens from disk (`file://`) and from any static host.
- `package.json` - `npm test` = `node --test scripts/*.test.mjs`; devDependency `manifold-3d` pinned to **3.5.3**
  (installed by `npm install`). The browser loads the same version: `https://cdn.jsdelivr.net/npm/manifold-3d@3.5.3/+esm`.
  Python canon runs manifold3d 3.5.2.
- `scripts/` - the equality test, `reference.py` (regenerates `fixtures/` from the Python canon in `reference/`, not needed for
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
| `derive(series, bed, depth)` | params object for a bed `{x, y, z}` in mm and a depth in mm (rule v3 below) |
| `limits(series, bed)` | `{depthMin, depthMax, depthDefault}` for the slider |

`params` always has the field names of the Python dicts: deep `N, pitch, length, dmax, fin, edge, floor, back, cap,
chamfer, bed`; blocks `Nx, Ny, W, H, hmax, base, bed`. The caller frees Manifold objects with `.delete()`.

## derive - rule v3 (deep wall thickness grows with depth), N is fixed

Tiles at the full bed size are unrealistic and tend to lift off the bed, so the margin is **M = 25 mm on every bed
axis**: every tile is at least 20 mm smaller than its stage 1 size. `PRESETS` stay the stage 1 dicts (the geometry canon
for the equality test); they are no longer reachable from the UI and no round trip to them is required. Beds in the UI:
Bambu Lab A1 mini 180 x 180 x 180, Prusa MK4 / Core One 250 x 210 x 220, Bambu Lab P1S / X1
256 x 256 x 256, Custom (x, y, z).

Deep (printed standing: width W on bed X, depth T on bed Y, length along the wells on bed Z):
`N = 7`; `pitch = floor((x - M) / 7)` whole mm; `floor 1.6, back 2.0, cap 1.2` fixed;
`fin = max(1.6, ceil01(dmax / 82.5))`, `edge = max(2.0, ceil01(dmax / 66))` where `ceil01(x)` rounds up to 0.1 mm
(`Math.ceil(x * 10 - 1e-9) / 10`, identical to the Python canon, the epsilon keeps an exact multiple from jumping up);
`chamfer = pitch < 30 ? 3 : 4`; `length = min(y - M, z - M)`; `dmax = depth`.
Why: a tile printed standing up to 231 mm tall swayed with 1.2 mm fins and 0.8 mm edge walls, so the walls get thicker
with the depth and stay as stiff on a deeper tile. Rule v2 (fin 1.2, edge 0.8 fixed) was the first public version; rule
v3 changes only these two numbers of the deep series (blocks and everything else are as in v2).
Limits: `dmax` from 40 to `y - M - back`; default = `min(4 * pitch, max)` (the "matched" tile:
unit depth = pitch, one tile = half the design wavelength).
Results: A1 mini pitch 22, tile 154 x 155, max 153, default 88; MK4 pitch 32, 224 x 185, max 183, default 128;
P1S pitch 33, 231 x 231, max 229, default 132. At the three default depths (and at every depth up to 132 mm)
fin 1.6 and edge 2.0; at the maximum depth of each bed fin / edge are A1 mini (153 mm) 1.9 / 2.4, MK4 (183 mm)
2.3 / 2.8, P1S (229 mm) 2.8 / 3.5.

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
| `wallBand(series, params, cols, rows)` | `{low, high, wallWidthMm, wallHeightMm}`: `high = band.f_high` (deep: min(f_width, N · f0); blocks: f_width), `wallWidthMm = cols * tile width` (deep: W; blocks: W), `low = max(start, 343 / (wallWidthMm / 1000))` in Hz, `start = band.f0`, except a deep wall narrower than `deepFullWidth(params) = 3.5 * band.lambda0` mm, where `start = sqrt(2) * band.f0` (2D BEM check of the default walls in the tile order until 25.10.2026: 3-5 tiles started about half an octave above f0, 7 tiles at f0; kept as a cautious start for the current order, where BEM starts at about f0). Below c / wall width a wall scatters little, so the low end is tied to wall width, never to f0 alone. |
| `features(manifold)` | `{area_cm2, volume_cm3, height_mm, footprint_cm2}` of a tile in print orientation (STL as built: Z up = print height; footprint = bbox X x Y). |
| `ESTIMATE` | coefficients from `scripts/calibration/estimate.json` (written by `scripts/calibration/calibrate.py`): `{deep: {grams: {...}, minutes: {...}}, blocks: {...}}`, each a map feature name (`area_cm2`, `volume_cm3`, `height_mm`, `footprint_cm2`, `const`) -> number. |
| `estimate(series, feats)` | `{grams, minutes}` = sum of coefficient x feature (+ const) per series, over the keys present in `ESTIMATE` only (a missing key counts as 0; the fit picked different features per series). Bambu P1S, 0.4 nozzle, Bambu PLA Basic, the stage 1 process profiles; shown in the UI as approximate. The deep coefficients are fitted on the 25 rule v3 slices (default depth of each bed plus four other depths up to the maximum), because with thicker walls the mass per volume of a deep tile changed (rule v2 coefficients were 9 % low); the older thin-wall deep rows stay in `data.json` with `in_fit: false`. |

## Equality ("the same STL")

For each of the 29 stage 1 tiles (params = preset, shifts from layouts.json) and for every case in
`scripts/fixtures/cases.json` (random params, Python reference STL made by `scripts/reference.py`): the JS mesh and
the reference STL, both measured by the same JS code - heights by ray cast at the measurement points of the Python reference tests
(+-0.01 mm), bounding box (+-0.01 mm), volume (+-0.1 %), shells: one closed outer solid, deep also the same number of
inner chamber shells as the reference.

## Bass traps core (two more script elements of index.html)

The bass trap mode is two classic scripts next to `diffuser-core`: `<script id="bass-physics">` defines the global function
`createBassPhysics()` (no manifold, no DOM, no `window`, no fetch, no import) and `<script id="bass-geometry">` defines
`createBassGeometry(M, physics)` (`M` as for `createDiffuser(M)`, `physics` = the object returned by `createBassPhysics()`; the
geometry uses only its `solveL` / `firstRes` to cut the channel, so there is one source of the channel length; the factory throws
when `physics` is missing). Both load in `node:vm` like `diffuser-core` (`scripts/mesh-tools.mjs` has the loader). They are
1:1 ports of the research of 25.10.2026 and of their Python twins `reference/bass_room.py` and `reference/bass_cell.py`
(headers name the source files and the date). Everything is **calculated, not measured**: `docs/bass-traps.md` has the maths, the
assumptions, the limits and the honest effect. Lengths in mm unless a name ends in `_m`; frequencies in Hz; decay rates `delta` in
1/s (pressure amplitude), `T60 = ln(1000) / delta`; `L = [Lx, Ly, Lz]` in m. The slit is 6 mm only, modes 35 - 120 Hz.
Nothing throws for ordinary input: a state that cannot be built is `{ok: false, reasons / error}`.

### createBassPhysics()

| member | meaning |
|---|---|
| `C0, RHO0, Z0, CELL, DETUNE_PCT (0.5), SLOT_W (6), N_LIST ([1,2,4,8,16]), ANCHORS, PEAK_POINTS (3001), PEAK_SPAN ([0.85,1.15])` | constants of the research; `ANCHORS` = `[{id, p, walls}]` with ids `corner, edge-z, edge-y, edge-x, wall-x, wall-y, wall-z` |
| `analyzeRoom(L, alpha, fmax = 120)` | `{L, alpha, fmax, modes, rows}`; `alpha` = one Sabine coefficient or six (walls x0, xL, y0, yL, z0, zL); `rows[i] = {index, n, f, kind ('axial'/'tangential'/'oblique'), Lam, deltaBase, t60Before, spacingHz, overlap}` |
| `evaluateMode(room, i, opts?)` | `opts = {ns, anchors (ids), L0 (mm), w}`; returns `{mode, deltaBase, t60Before, gammaCell, fCell, L0, noEffect, overlap, uncertain: [reasons], anchors: {id: {text, psi2, nodal, cellsAvailable, byN: {N: {available, solved, deltaEff, fEff, t60After, decayFactor, peakDropDb, deltaWallsEff, G, approximate, approximateReasons}}}}, best: {N: anchorId}}`. `best[N]` = the anchor with the largest `deltaEff` (ties: first of `ANCHORS`); `decayFactor = deltaEff / deltaBase = T60 before / T60 after`; `peakDropDb` < 0 = lower peak; a trap with more cells than the walls have gives `available: false`, an unsolved equation `solved: false` |
| `roomModes(L, fmax, c?)`, `baseDecay(mode, L, alpha, c?)`, `sabineDecay(L, alpha, c?)`, `t60(delta)` | modes, base decay of a mode, Sabine decay, T60 |
| `cluster(L, anchorId, nmax)`, `cellsPsi(modes, L, cells)`, `loadingMass(N, anchorId)` | the nmax cells nearest an anchor (0.12 m grid), `<psi>` averaged on the cell faces `[cell][mode]`, radiation mass kg/m2 (corner 2, edge sqrt 2, wall 1 times the piston mass) |
| `zsQwr(omr, omi, w, L, Y, rough)`, `makeCell(fDesign, w, L, extraMass?, rough?)`, `cellForMode(fMode, L0?, w?)` | cell impedance `[re, im]` at complex angular frequency; a cell object `{zs(omr, omi), dzs(om), alphaAt(f), gamma(f?)}`; the cell of the design of a mode (`f cell = f mode x 1.005`, `L0 = solveL(f cell, w)`) |
| `firstRes(w, L, rough?)`, `solveL(f0, w, rough?)` | `design_Q.first_res` / `solve_L` on the grid 30, 30.25, ..., 519.75 Hz (NaN when no resonance on the grid); 40 bisection steps |
| `solveSingle(wn, deltaW, G, cell, c?)`, `singleModeEffect(mode, deltaWEff, G, cell, c?)` | roots of `(w^2 - wn^2 - 2 j w delta') Zs(w) = j w c Z0 G` (Newton polish on the exact Zs); `deltaEff` = Im of the slower of the two roots nearest the mode (`null` when none) |
| `peakDropDb(mode, deltaBase, deltaWEff, G, cell, c?, grid?)`, `peakGrid(mode, cell)` | the peak drop defined in `docs/bass-traps.md` (grid optional: the unloaded Zs on the frequency grid, to reuse per mode) |
| `overlapFlag(modes, i, L, alpha)`, `orient(kind, overlap, N)`, `approximateReasons(...)`, `uncertainReasons(kind, alphaMean, deltaBase, gamma)`, `anchorText(id, L)` | flags (`approximate` = overlap or N > 8 or non-axial with N >= 4 or a corner with a neighbouring mode within 10 %), English description of a place |
| `dlSlitArray(h_m, d_m, nterms?)`, `absorption([re, im])` | tmm.py helpers |
| `verdict(ev, opts?)`, `VERDICT_WORTH (1.25)`, `VERDICT_SMALL (1.08)`, `VERDICT_MAX_N (8)` | the verdict of the mode table of the interface (policy, not physics, no Python twin; see "Mode verdict" below) |

### Mode verdict (interface policy)

`verdict(ev, {maxN = 8})` takes the result of `evaluateMode` (at the alpha of the room, any `ns`) and sorts the mode by the best decay factor
`x` of every place of every N <= `maxN` (cases with `available: false` or `solved: false` do not count; ties go to the smaller N, then to the
earlier place of `ANCHORS`):

| level (`label`) | rule |
|---|---|
| `worth` ("Worth trying") | `x >= 1.25` and no `noEffect` flag |
| `small` ("Small effect") | `1.08 <= x < 1.25` and no `noEffect` flag |
| `none` ("No effect") | `x < 1.08`, or no solved case (`x = null`), or the `noEffect` flag of `evaluateMode` (delta of the walls above the width of the cell) |

It returns `{level, label, x, N, anchor, noEffect, approximate, approximateReasons, uncertain, uncertainReasons}`: `x`, `N`, `anchor` are the best
case (`null` when none could be solved); `approximate` and its reasons are the flags of that best case (overlapping modes, tangential or oblique
with 4 or more cells, a corner with a neighbouring mode), `uncertain` and its reasons are `ev.uncertain` (delta above 0.8 of the cell width, oblique
at alpha >= 0.10, alpha >= 0.15). The thresholds are the owner's decision of 25.10.2026 and a policy of the interface: they are not derived from the
physics and are not compared with a Python twin. The test is the last block of `scripts/bass-physics.test.mjs`.

### createBassGeometry(M, physics)

Data keys follow the Python twin (`plan_json`): `design` returns `{ok, L0, f0, f_mode, plan: {n, y_in, p_nom, bends, L_geom}, Lk,
pk, ns_down, ns_up, warn, w, wall, width, height, bend_dL, depth, p_limit, ...}`; rectangles are `{kind, x: [x0, x1], y: [y0, y1], z:
[z0, z1]}` in block coordinates (x depth from the face, y face width, z height = print axis).

| member | meaning |
|---|---|
| `DEF, BEND_K (-0.245), BEDS, BED_MARGIN (5), MAX_HOURS (12), PLATE_GAP (3), BLOCK_ORDER, FMODE_MIN, FMODE_MAX, ESTIMATE` | defaults of the research (w 6, wall 1.2, face 120 x 120, inserts -2 ... +8 steps of 0.5 %), the beds of `createDiffuser(M).BEDS`, layouts tried in order (2x2 k4, 2x2 k3, 1x2, 2x1, 1x1), the fitted coefficients |
| `design(fMode, opts?)` | the cell for a MODE frequency: `L0`, legs, inlet, the inserts `k = -ns_down ... ns_up` (`Lk`, `pk` in mm), warnings; `{ok: false, error}` outside 35 - 120 Hz |
| `planBlock(d, ny, nz, k)` | block of `k` cells in an `ny x nz` grid (filled in rows from the bottom): `size_mm` [T, ny 120, nz 120], per cell `air`, `walls`, `voids`, `axis_points_xy`, `window`, plus `volume_mm3`, `x_last`, `hc`, `L_axis`, `L_acoustic` |
| `blockSolid(d, bp)`, `insertSolid(d, bp, k)`, `platesSolid(d, bp, plate)` | `M.Manifold` (the caller frees it with `.delete()`): the block (hull minus channel air minus closed pockets), an insert in print orientation (flange on the bed, top edges chamfered 0.6), a plate |
| `platePlan(d, bp)`, `insertDims(...)`, `insertSuffix(k)` | the research plate: inserts in a row along x, 3 mm apart, with positions and dimensions; the insert dimensions |
| `packBlocks(d, N, bed, opts?)` | the blocks that carry N cells on `bed = {x, y, z}` (5 mm margin per axis, `z` is the print axis, `x` and `y` may swap): `{ok, reasons, blocks: [{ny, nz, k, count, size_mm, minutes, grams, hours, overLimit, name}], prints, totalMinutes, totalGrams, notes}`; the first layout of `BLOCK_ORDER` that fits and takes at most 12 h, a remainder block for the rest, `overLimit` only when nothing quicker fits; `ok: false` with reasons when one cell is deeper than 150 mm or than the bed, or higher than the bed |
| `packPlates(d, bp, bed, spec?)` | `spec = {caps (copies of the nominal insert, one end cap per cell; default 1), tuning (one of every other insert; default true)}`; `{ok, reasons, plates: [{items: [{k, p, x0, y0}], rotated, size_mm, minutes, grams, name}]}`; plates split and turn 90 degrees where that holds more |
| `features(m)`, `planFeatures(bp)`, `estimate(kind, feats)` | quantities of a solid (`area_cm2, volume_cm3, height_mm, footprint_cm2`); the same without a solid (volume, height, footprint); `{grams, minutes}` for `kind` `'block'` or `'plate'` |
| `stem(d)`, `blockName(d, bp)`, `insertName(d, k)`, `platesName(d, index, count)` | file stems: `soundblooms-bass-trap-qwr-f095.2-w6-block-1x2-k2-58.8x120x240`, `...-insert-m2` / `-0` / `-p3`, `...-inserts-plate[-2of3]` |
| `axisLength, planCell, cellLayout, axisPolyline, unionArea, blockCells, bedLimits, fitsBed` | ports of the helpers of `qwr_cell.py`; bed helpers |

The nominal insert (`k = 0`) closes the channel: it is part of every cell. The inserts `k != 0` are optional tuning parts, to be used
only with a measurement of the decay before and after.

### Print estimate

`ESTIMATE` (`scripts/calibration/bass-estimate.json`, written by `scripts/calibration/bass-calibrate.py` from `bass-data.json`):
per kind (`block`, `plate`) and quantity (`grams`, `minutes`) a map feature -> coefficient, fitted to the 10 blocks and 10 inserts
plates of the research set (Bambu Studio, P1S, 0.4 nozzle, PLA, process "Diffuser blocks 0.20"). The block uses only features that
follow from the rectangle plan (volume, height, footprint), so packing needs no solid. Shown as approximate.

### Equality and fixtures

`scripts/bass_reference.py physics|geometry` runs the Python twins and writes `scripts/fixtures/bass/`: `physics.json` (5 rooms
x alpha 0.05 / 0.10 / 0.20 x every mode up to 120 Hz x 7 anchors x N 1, 2, 4, 8, 16, and building blocks), `geometry.json` (designs
for 9 modes 35 - 120 Hz, rectangle plans of blocks 1x1, 1x2, 2x1, 2x2 with 3 and 4 cells, the plate plans, STL statistics) and binary
STL of 5 blocks, 3 inserts and the plate. Tolerances: physics numbers 1e-6 relative, roots of the Newton polish 1e-9, rectangle
plans 1e-6 mm; JS meshes against the reference STL measured by the same code as in `equality.test.mjs` (copied into
`scripts/mesh-tools.mjs`): volume 0.1 %, bounding box 0.01 mm, heights / sections by ray cast 0.01 mm, one closed outer solid with
the pockets as the only cavities. The tests are `scripts/bass-physics.test.mjs` and `scripts/bass-geometry.test.mjs` (`npm test`).
