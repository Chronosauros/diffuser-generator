# Soundblooms Diffuser Generator

Browser generator of the Soundblooms acoustic diffuser tiles (deep 1D QRD and block 2D QRD). Use it online at
https://chronosauros.github.io/diffuser-generator/

One HTML file, `index.html`, opens from disk or any static host. Pick a series, how many tiles (the wall; a column
chart with faces compares how much sense each wall option makes), a printer bed (tiles stay 25 mm inside the bed on
every axis) and a depth. The page shows each tile of the wall in 3D, a summary of the chosen wall (size, how much
sense it makes, approximate print time and filament weight, calibrated on Bambu P1S slices), the calculated frequency
band, what you are about to print (the whole wall and the tile on screen, each with its time and filament) and the
wall plan, and downloads print-ready binary STL files (the chosen wall as a ZIP with one STL per tile number, the wall
plan SVG and a README.txt with the quantity of each tile, or the tile on screen alone).

The geometry core is the `<script id="diffuser-core">` section: a 1:1 port of the Python reference implementation of
the tile geometry, plus the bed/depth rules in `CONTRACT.md`. It has no DOM and runs in `node:vm` as well as in the
browser. The page loads manifold-3d 3.5.3 and three.js 0.160.0 from the jsDelivr CDN, so it needs a connection the
first time.

## Development

With Node.js 18 or newer, `npm install` and `npm test` run the equality tests (`scripts/`): the JS core against the
Python reference STL files for all 29 reference tiles (`scripts/fixtures/stage1/`) and for random parameter sets
(`scripts/fixtures/`). They do not cover real browser rendering.

The Python reference implementation itself is not part of this repository. `scripts/reference.py` regenerates the
fixtures from it and takes the folder that holds it as an argument (`python scripts/reference.py <design-folder>`, or
the `DIFFUSER_DESIGN` environment variable); it is not needed for `npm test`. `scripts/calibration/` holds the data and
the fit behind the print time and filament estimate; `calibrate.py` there documents how it was made and needs the same
private project.

## Licence

The generator and the diffuser tiles it makes are licensed under Creative Commons Attribution-NonCommercial-ShareAlike
4.0 International (CC BY-NC-SA 4.0): https://creativecommons.org/licenses/by-nc-sa/4.0/

## Credits

Set in Geom, copyright 2022 The Geom Project Authors (https://github.com/ThanosPoulakidas/Geom), SIL Open Font License
1.1 (https://openfontlicense.org); the font is embedded in `index.html`. Geometry with
[manifold-3d](https://github.com/elalish/manifold), preview with [three.js](https://threejs.org). Formulas after
T. J. Cox and P. D'Antonio, *Acoustic Absorbers and Diffusers*, 3rd ed., CRC Press, 2016. The frequency band is a
calculated design band, not a measurement.

Made by [Soundblooms](https://soundblooms.com).
