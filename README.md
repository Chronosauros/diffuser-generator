# Diffuser Generator

A free browser tool that makes printable acoustic diffuser tiles (quadratic residue diffusers) sized to your 3D
printer. Pick a series, a wall and your printer bed, check each tile in 3D and download print-ready STL files.

Use it online: https://chronosauros.github.io/diffuser-generator/

The generator was built with AI assistance and the tiles have not been measured. The maths follows Cox and D'Antonio,
but I would like people who know acoustics to check it. If you see a mistake or a better way, open an issue or start a
discussion.

## What it does

- Two series: deep tiles (1D QRD, N = 7, wells with fins) and block tiles (2D QRD, 13 x 13 or 13 x 11, shallow).
- Sizes every tile from the printer bed, with 25 mm to spare on every axis. Presets for the Bambu Lab A1 mini, Prusa
  MK4 / Core One and Bambu Lab P1S / X1, or any custom bed.
- Numbered tiles in a set order, so a wall grows without moving the tiles that are already up.
- Shows the design band (design intent, calculated, not measured), a rough print time and filament estimate, and the
  wall plan.
- Downloads the whole wall as a ZIP (one STL per tile number, the wall plan as SVG and a README with quantities) or a
  single tile.

One HTML file, `index.html`, that opens from disk or any static host. It loads manifold-3d 3.5.3 and three.js 0.160.0
from the jsDelivr CDN, so it needs a connection the first time.

## How it works

The maths, the sizing rules, the wall orders and the open questions are written up in
[docs/how-it-works.md](docs/how-it-works.md). In short:

- deep tiles: well depths d_n = s_n · λ0 / (2N) with s_n = n² mod 7, design frequency f0 = c · s_max / (2N · d_max),
  upper limit c / (2w); one tile is one period and the seven tiles are the seven cyclic shifts of the sequence;
- block tiles: phase frac(s_n / Nx + s_m / Ny) turned into block heights;
- wall orders: picked with the Fourier far-field model and the AES-4id diffusion coefficient, used only to compare
  layouts with each other.

## Repository layout

- `index.html` - the whole app. Its geometry core, `<script id="diffuser-core">`, has no DOM and also runs in Node.
- `reference/` - the Python reference implementation: formulas, tile geometry, the far-field model and the layout
  optimiser. The generator is a port of it.
- `scripts/` - the equality tests, the fixtures made by the Python code, and the data behind the print estimate.
- `docs/how-it-works.md` - the maths and the design rules.
- `CONTRACT.md` - what the JavaScript core must match.

## Development

With Node.js 18 or newer, `npm install` and `npm test` run the equality tests: the JavaScript tiles against the Python
reference STL files for all 29 reference tiles and for random parameter sets (heights within 0.01 mm, bounding box
within 0.01 mm, volume within 0.1 %, one closed shell). They do not cover real browser rendering.

To change the geometry, change `reference/` and the core in `index.html` together, regenerate the fixtures with
`python scripts/reference.py` and run `npm test`. Python setup: [reference/README.md](reference/README.md). The print
time and filament fit in `scripts/calibration/` was made from Bambu Studio slices; `calibrate.py` documents how, but
needs the slicer setup it was made with.

Contributions: [CONTRIBUTING.md](CONTRIBUTING.md).

## Licence

The generator and the Python reference are licensed under the GNU General Public License v3.0 ([LICENSE](LICENSE)).
You may use, change and share them, also commercially, as long as anything you distribute that is based on them stays
under GPL-3.0 with its source. The ready-made tile files published on Printables have their own licence, CC BY-NC-SA
4.0, stated on the model page. The embedded Geom font is under the SIL Open Font License 1.1.

## Credits

Set in Geom, copyright 2022 The Geom Project Authors (https://github.com/ThanosPoulakidas/Geom), SIL Open Font License
1.1 (https://openfontlicense.org); the font is embedded in `index.html`. Geometry with
[manifold-3d](https://github.com/elalish/manifold), preview with [three.js](https://threejs.org). Formulas after
T. J. Cox and P. D'Antonio, *Acoustic Absorbers and Diffusers*, 3rd ed., CRC Press, 2016. The frequency band is
design intent, calculated from the shape of the tiles, not a measurement.

Made by [Soundblooms](https://soundblooms.com).
