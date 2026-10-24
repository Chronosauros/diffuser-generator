# How the Diffuser Generator works

This page describes the maths and the design rules behind the generator, so that anyone can check them. Every number
here is calculated from formulas. Nothing has been measured: the tiles have not been tested in a room or to ISO 17497-2.

The generator, the reference implementation and this page were written with AI assistance (Claude). The formulas come
from T. J. Cox and P. D'Antonio, *Acoustic Absorbers and Diffusers: Theory, Design and Application*, 3rd ed., CRC Press,
2016. If something here is wrong, please open an issue.

## What the tiles are meant to do

A flat wall sends a reflection back as one strong copy of the sound. A Schroeder (phase grating) diffuser is a row of
wells of different depths. Sound reflected from a deep well comes back later than sound from a shallow one, so the
parts of one reflection leave the surface with different phases. The depths follow a number sequence whose Fourier
transform is flat, which spreads the reflected energy evenly over the diffraction directions instead of piling it
into the specular one.

It is a diffuser, not an absorber. It is not meant to do anything for bass, room modes or the reverberation time.

## Deep series: 1D quadratic residue diffuser, N = 7

Notation: N is a prime, n = 0..N-1 is the well index, c = 343 m/s, w is the well width (pitch minus fin).

| Quantity | Formula |
|---|---|
| Sequence | s_n = n² mod N |
| Well depth | d_n = s_n · λ0 / (2N), with λ0 = c / f0 |
| Design frequency | f0 = c · s_max / (2N · d_max) |
| Upper limit (cross modes in a well) | f_w = c / (2w) |
| Flat-plate frequency | f_N = N · f0 and its multiples |

For N = 7 the sequence is 0 1 4 2 2 4 1, so s_max = 4 and the deepest well is 4/14 of the design wavelength. N = 7 was
chosen because it gives the best ratio s_max / N among the primes 5 to 13, which means the lowest f0 for a given depth.
That matters because the depth is limited by the printer bed.

### Sizing rule

A tile is exactly one period: seven wells side by side. The generator sizes a tile from the printer bed (x, y, z in mm)
with a margin of 25 mm on every axis:

- pitch = floor((x - 25) / 7) mm (well plus fin)
- length along the wells = min(y - 25, z - 25) mm
- default depth d_max = 4 · pitch, so one depth unit equals the pitch and one tile is exactly λ0 / 2 wide (a "matched"
  tile). The depth slider can make the tile deeper, up to the bed.

Because of that rule, f0 = c / (14 · pitch) and the flat-plate frequency N · f0 = c / (2 · pitch). The well width w is
the pitch minus a fin, so N · f0 always falls just below f_w = c / (2w). At the default depth the flat-plate frequency
is the top end of the band.

| Printer | Tile (mm) | Pitch | Interior well | d_max | f0 | N · f0 | f_w |
|---|---|---|---|---|---|---|---|
| Bambu A1 mini (180 x 180 x 180) | 154 x 155 x 90 | 22 mm | 20.4 mm | 88 mm | 1114 Hz | 7.80 kHz | 8.4 kHz |
| Prusa MK4 / Core One (250 x 210 x 220) | 224 x 185 x 130 | 32 mm | 30.4 mm | 128 mm | 766 Hz | 5.36 kHz | 5.6 kHz |
| Bambu P1S / X1 (256 x 256 x 256) | 231 x 231 x 134 | 33 mm | 31.4 mm | 132 mm | 742 Hz | 5.20 kHz | 5.5 kHz |

### Band shown by the generator

The generator shows a design band for one tile from f0, set by the deepest well, to the lower of f_w and the
flat-plate frequency N · f0 (at the default depth that is N · f0; a shallow tile from the depth slider can push N · f0
above f_w). It says which well sets each end and calls the band design intent, calculated, not measured. For a wall it
ties the lower end to the wall width as well: low = max(f0, c / wall width). For walls of three tiles or more c / wall
width is already below f0, so the lower end is f0. f0 is not a hard edge: a diffuser keeps some effect below it, and
above N · f0 it diffuses again outside the design theory.

Until 24.10.2026 the top of the band was f_w, which put it above the flat-plate frequency (P1S / X1: 5.46 kHz instead of
5.20 kHz). Thanks to the review on r/acoustics for pushing on "design intent, not a spec".

## Seven tiles: cyclic shifts of one sequence

Seven identical tiles in a row would repeat the same period, and a periodic surface concentrates energy in grating
lobes at sin θ + sin ψ = m · λ / W, where W is the tile width. To break that period, the seven tiles are the seven
cyclic shifts of the sequence. Tile t uses shift k = (t + 3) mod 7, so tile 1 is the symmetric period centred on the
zero-depth well:

| Tile | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| Shift k | 4 | 5 | 6 | 0 | 1 | 2 | 3 |
| Depth units (STL X order) | 2 4 1 0 1 4 2 | 4 1 0 1 4 2 2 | 1 0 1 4 2 2 4 | 0 1 4 2 2 4 1 | 1 4 2 2 4 1 0 | 4 2 2 4 1 0 1 | 2 2 4 1 0 1 4 |

Turning tile t upside down reverses its wells and gives tile 9 - t, so 2 and 7, 3 and 6, 4 and 5 are the same part.
That means four distinct shapes.

Why the shifts should help (my own derivation in the Fourier model, not taken from the literature): a cyclic shift by
k wells moves the content of one period by k · pitch, which is k / N of a period. In a grating lobe direction of order
m, the far field of the shifted tile is the far field of the unshifted tile times exp(-i 2π m k / N). A row that
contains each of the seven shifts once sums these factors to zero for m = 1..6, so in this model the lobes of the tile
period cancel in those exact directions and their energy goes elsewhere. With only three, four or five tiles the sum
is not zero, and the order of the tiles also changes the response in every other direction. The order was therefore
picked by an optimiser (next section), not by the derivation alone.

Wall orders seen from the room, left to right. A wider wall adds tiles on the left, so tiles already on the wall never
move:

| Printer | 3 tiles | 4 tiles | 5 tiles | 7 tiles |
|---|---|---|---|---|
| A1 mini | - | 4 5 3 6 | - | 1 7 2 4 5 3 6 |
| MK4 / Core One, P1S / X1 | 4 3 6 | - | 5 2 4 3 6 | 1 7 5 2 4 3 6 |

More rows repeat the same row.

## Block series: 2D quadratic residue diffuser

The block series is a grid of square blocks of different heights, without fins (skyline style). The phase of cell
(n, m) is the sum of two 1D residue phases, because the reflection factors of two 1D sequences multiply:

- φ(n, m) = frac(s_n / Nx + s_m / Ny), which for Nx = Ny = N is ((n² + m²) mod N) / N
- block height above the 3 mm base: h = h_max · (1 - φ / φ_max), so the tallest block is the zero-phase cell
- design frequency f0 = c · φ_max / (2 · h_max)
- rough upper limit c / (2 · block width). Without fins this is not a well cross mode, only the point where a block is
  half a wavelength wide.

Sizes: 13 x 13 cells on the A1 mini and P1S / X1, 13 x 11 on the MK4 / Core One (the bed is not square there).
Each tile is the same pattern shifted by (k, l) cells, a different shift per tile number, with the shifts and the
grid order picked by the same optimiser.

| Printer | Tile (mm) | Cells | Block | h_max | f0 | c / 2w |
|---|---|---|---|---|---|---|
| A1 mini | 155 x 155 x 50 | 13 x 13 | 11.9 mm | 47 mm | 3.4 kHz | 14.4 kHz |
| MK4 / Core One | 225 x 185 x 63 | 13 x 11 | 17.3 x 16.8 mm | 60 mm | 2.8 kHz | 10.2 kHz |
| P1S / X1 | 231 x 231 x 51 | 13 x 13 | 17.8 mm | 48 mm | 3.3 kHz | 9.7 kHz |

The block series is the shallow, decorative one: upper mids and treble only.

## How the wall orders were chosen

The orders come from a far-field prediction model from Cox and D'Antonio, chapter on prediction models: the Fourier
(Fraunhofer) model of a plane surface with a position-dependent reflection factor R(x).

- Deep tiles: R = exp(-2ikd) for a well of depth d (a plane wave down the well and back), R = 1 on the fins. The far
  field is p(θ) = Σ_e w_e · sinc(k w_e s / 2) · exp(i k x_e s) · R_e with s = sin θ + sin ψ.
- Block tiles: each block top is a flat facet at height h with a Kirchhoff phase exp(i k h (cos θ + cos ψ)).
- Diffusion coefficient as in AES-4id-2001 / ISO 17497-2, d = [(Σ E)² - Σ E²] / [(n - 1) Σ E²], over equally weighted
  receivers (1D: a semicircle every 5 degrees; 2D: 400 equal-area points on the hemisphere), with E averaged over the
  frequencies of each third-octave band.
- Normalised to a flat plate of the same size, d_n = (d - d_flat) / (1 - d_flat), per incidence, then averaged over
  incidences of 0, ±30 and ±60 degrees.
- The optimiser tries the orders and keeps the one with the best mean d_n over all wall sizes of a series together, so
  a wall can grow without moving tiles. For the deep series the third-octave band around N · f0 is left out of the score.

This model is only used to compare layouts with each other. It is known to be weakest at low frequencies, at oblique
angles and near the edges of small walls, so I do not quote any diffusion coefficient from it.

## Construction details that may matter acoustically

- Deep tiles print standing on their end, wells upward, so every layer has the same cross-section and nothing needs
  supports.
- Fins are 1.6 mm and the outer walls 2.0 mm of PLA at the default depths, thicker for a deeper tile (fin =
  max(1.6, d_max / 82.5), edge = max(2.0, d_max / 66), rounded up to 0.1 mm). Two tiles side by side make a 4.0 mm fin,
  so the two wells at a joint sit 1.2 mm further apart than the pitch.
- The solid under the shallower wells is hollowed into closed air chambers (1.6 mm floor under each well, 1.2 mm cap on
  the bed, ridge roof of two 45 degree slopes so nothing bridges over air). This saves a lot of filament.
- 45 degree chamfers of 3 or 4 mm join each well floor to its fins. They change the effective depth of a well by about
  0.5 mm at most.
- Back plate 2 mm, flat, for adhesive strips.

## What the generator does with all this

`index.html` holds the whole app. Its geometry core (`<script id="diffuser-core">`) is a port of the Python reference
implementation in `reference/` (sequences, depths, bands, sizing rules and tile geometry, built with manifold-3d). `npm
test` checks the JavaScript tiles against STL files made by the Python code: heights at the measurement points within
0.01 mm, bounding box within 0.01 mm, volume within 0.1 %, one closed shell.

The "how much sense" faces are a plain rule of thumb by wall width (one tile, under 60 cm, 60 to 120 cm, from
120 cm), not a calculation. Print time and filament come from a linear fit to Bambu Studio slices on a P1S profile.

## Open questions

These are the parts I am least sure about:

1. Is N = 7 a sensible choice under a print-depth limit, or would a larger prime (11, 13) with shallower steps work
   better on these widths?
2. Are cyclic shifts of one sequence a valid way to break periodicity, compared with the usual modulation by two
   different base shapes?
3. The flat-plate frequency N · f0 sits just below f_w because of the matched sizing rule. Should the rule move it,
   or is that the normal compromise?
4. Is c / wall width a fair lower limit for a small wall, or too optimistic?
5. Can 1.6 mm PLA fins, about 20 to 31 mm wide wells and closed air chambers under 1.6 mm floors add absorption or
   resonances that the model does not see?
6. Is the Fourier model good enough to choose the orders, or should that be done with BEM?

## References

- T. J. Cox, P. D'Antonio, *Acoustic Absorbers and Diffusers: Theory, Design and Application*, 3rd ed., CRC Press, 2016.
- AES-4id-2001, AES information document for room acoustics and sound reinforcement systems: characterisation and
  measurement of surface scattering uniformity.
- ISO 17497-2:2012, Acoustics: sound-scattering properties of surfaces, part 2: measurement of the directional
  diffusion coefficient in a free field.
- T. J. Cox, "A brief history of room acoustic diffusers", Acoustics Today.
