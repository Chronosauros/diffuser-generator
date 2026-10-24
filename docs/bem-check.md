# 2D BEM check of the deep tiles

After the review on r/acoustics ("a pretty plot at 2 kHz says nothing about 200 Hz"; "without a measurement the band is
design intent") the deep tiles and walls were run through a 2D boundary element method (BEM) simulation and compared
with the Fourier model that picked the wall orders. BEM is what Cox and D'Antonio use to check diffusers before a
measurement. Everything here is a simulation, done on 24-25.10.2026 with AI assistance (Claude). Nothing has been
measured, and the numbers below are model predictions, not performance claims. The code, the raw far fields and every
figure are in [`checks/bem/`](../checks/bem/), so anyone can rerun it or point out a mistake.

## Short answer

- Below f0 no wall diffuses. Between 200 Hz and f0 the normalised diffusion coefficient d_n is at most 0.10: the wall
  reflects like a flat panel of the same width. The reviewer was right about 200 Hz.
- The cyclic shifts work, but not from c / wall width. A wall in the published order starts to diffuse (d_n >= 0.3 in
  two third-octave bands in a row) at about f0 with 7 tiles (P1S / X1: the 800 Hz band, f0 = 742 Hz; A1 mini: 1000 Hz,
  f0 = 1114 Hz), a third of an octave higher with 3 and 5 P1S tiles (1000 Hz) and half an octave higher with 4 A1 mini
  tiles (1600 Hz). Between c / wall width and f0, d_n is 0.07-0.10.
- Between f0 and c / period (= 2 f0, a tile is λ0 / 2 wide) the order of the tiles makes a large difference: walls in
  the published order reach d_n 0.43-0.48, walls made of copies of one tile 0.13-0.26, and at normal incidence the
  copies hardly diffuse at all (0.02-0.14). It is the same physics as in the standard QRD below: a periodic surface has
  only one propagating reflection direction below c / period, the specular one.
- A single tile only diffuses from 2 f0 = c / tile width. The Fourier model says it diffuses from f0 (d_n about 0.8),
  BEM says it does not (0.07-0.21 on average, less than a flat panel at 45°).
- The Fourier model ranks things the same way as BEM: the published order beats copies of one tile in every range
  from f0 up, and both models show the flat-plate frequency 7 f0. In numbers it is optimistic: it agrees with BEM for
  7 tiles between f0 and 2 f0, overestimates d_n of 3-5 tiles there (P1S, 3 tiles: 0.78 against 0.48), overestimates
  by 0.1-0.2 between 2 f0 and 7 f0 (by 0.22-0.36 at normal incidence) and underestimates above 7 f0 at oblique
  incidence.
- Ranking all 210 orders of a 3-tile P1S / X1 wall: the published order `6 3 4` is first in the Fourier model and 22nd
  of 108 mirror pairs in BEM (d_n 0.41, the BEM best `5 6 2` reaches 0.45). The two rankings correlate moderately
  (Spearman 0.68). The Fourier model is good at ruling out bad orders, not at picking the best one.

What changed in the generator because of this: the lower end of a deep wall's band. A deep wall narrower than 3.5
design wavelengths (fewer than 7 tiles at the default depth) now starts at √2 · f0 instead of f0, never below
c / wall width (commit `ac8c5f2`, see "Band shown by the generator" in [how-it-works.md](how-it-works.md)). A single
tile starts at c / tile width = 2 f0, as before. The top of the band at the flat-plate frequency 7 f0 (commit
`2278725`) is what BEM shows as well.

![P1S / X1 walls: d_n against the band, BEM and Fourier](../checks/bem/out/fig-P1S-summary.png)

![A1 mini walls: d_n against the band, BEM and Fourier](../checks/bem/out/fig-A1mini-summary.png)

Solid lines BEM, dashed lines the Fourier model; blue the published order, orange copies of tile 1. Shaded bands: the
flat reference itself already scatters almost evenly there (d_flat > 0.9), so d_n divides by almost zero.

## Method

Solver: [`checks/bem/bem2d.py`](../checks/bem/bem2d.py), about 200 lines of numpy and scipy. No ready tool fitted:
bempp-cl is 3D only, AcousticBEM is unmaintained and has no near-singular integration (a 1.6 mm fin has its other face
1.6 mm away), NiHu and OpenBEM need C++ or MATLAB.

- Exterior Helmholtz equation around a rigid (sound-hard) body, plane-wave incidence, far field.
- Burton-Miller formulation (CBIE plus i/k times the hypersingular equation), so the interior resonances of the closed
  body do not spoil the result. Without it the P1S reference box (231 x 134 mm) is 3.4 % off at its first interior
  resonance, 1480 Hz, and a 7-tile wall has such resonances all over the band.
- Constant elements on straight segments, collocation at the midpoints. The singular and near-singular parts of the
  Laplace kernels are integrated exactly over every segment, the smooth Helmholtz-minus-Laplace remainder by Gauss
  quadrature; the hypersingular operator in Maue's form.
- Element length at most min(λ / 10, 2 mm). The elements must be shorter than the fins are thick: a first run with
  λ / 8 capped at 8 mm was 13 % off in the diffusion coefficient of a P1S tile at f0.
- Geometry from `reference/diffuser.py`, rule v3: fins 1.6 mm, outer walls 2.0 mm (a 4.0 mm fin where two tiles meet),
  45° chamfers in the bottom corners of every well (3 mm A1 mini, 4 mm P1S), flat back. The closed air chambers under
  the shallow wells are invisible from outside.
- Wall orders: the `walls` of `reference/out/layouts.json` in model (STL X) coordinates, as in `reference/layout.py`:
  A1 mini class S (4 tiles `6 3 5 4`, 7 tiles `6 3 5 4 2 7 1`), P1S / X1 class L (3 `6 3 4`, 5 `6 3 4 2 5`, 7
  `6 3 4 2 5 7 1`). Copies: tile 1 repeated, as in `layout.evaluate_deep`.
- Reference: a rigid box with the same outline (wall width x tile depth, wells filled), free-standing.
- Diffusion coefficient AES-4id / ISO 17497-2 over 37 receivers every 5°, energy averaged over the 6 frequencies of
  each third-octave band 200 Hz - 10 kHz (the same frequencies as the evaluation in `layout.py`), normalised per
  incidence: d_n = (d - d_flat) / (1 - d_flat). Incidences 0, ±30, ±45 and ±60°. `diffusion` and `band_energy` are
  imported from `layout.py`, and the Fourier model is called through its own code (`layout.Deep._tile`, `wall`,
  `flat`) on the same frequencies and receivers.
- Polar plots: dB re the peak of the flat reference of the same model at that frequency and incidence. Positive angles
  are on the +X side of the STL (seen from the room, STL X points left).

### What a rigid 2D BEM does not see

- 2D means infinitely long wells. A real tile is 155-231 mm along its wells, so in the plane along the wells it acts
  like a panel of that length. That is the usual assumption for 1D diffusers (single-plane measurement).
- Rigid walls: no viscous losses and no fin vibration. A real diffuser absorbs a little, mostly near the well
  resonances (deepest well a quarter wave at about 650 Hz on the P1S / X1, 975 Hz on the A1 mini) and in the narrow A1
  mini wells (20 mm, below the 25 mm minimum QRDude recommends). Where energy circulates between the wells, a real wall
  may reflect less than this BEM, not diffuse more.
- A free-standing sample, not one set into an infinite wall: that is how the diffusion coefficient against a flat panel
  of the same size is defined.

## Validation of the solver

[`checks/bem/validate.py`](../checks/bem/validate.py), numbers in `checks/bem/out/validation.json`.

| Test | Result |
|---|---|
| Rigid cylinder against the exact Bessel series (far field), ka = 1 / 10 / 30 | error 0.008 % / 0.11 % / 0.36 % |
| Cylinder at an interior resonance (ka = 2.405, zero of J0) | Burton-Miller 0.02 %, plain CBIE 0.12 % |
| P1S box 231 x 134 mm, 1440-1520 Hz (interior resonance 1479.6 Hz) | Burton-Miller at most 0.05 %, CBIE 3.4 % at 1480 Hz |
| Flat box 1617 x 134 mm against Kirchhoff, kW > 20 | specular peak within 0.2 dB, half main lobe within 0.1 dB (0°) and 1.1 dB (45°) |
| Flat box against Kirchhoff at small kW (one 231 mm tile, 200 Hz) | 0.4-5 dB apart: Kirchhoff does not hold there, as expected |
| Mesh: 2 mm elements and 4 Gauss points against 1 mm and 8 (tile, box, 7-tile P1S wall, f0 to 10 kHz) | diffusion coefficient within 0.004, levels above -20 dB within 0.4 dB |

A finite flat panel has no exact closed-form solution, so the cylinder is the strict test and Kirchhoff only an
asymptotic one (large kW, near the specular direction).

![Solver checks: cylinder and flat box](../checks/bem/out/fig-validation-solver.png)

### Standard N = 7 QRD (Cox, Angus, D'Antonio, Proc. IOA 25(5), 2003)

f0 490 Hz, depth unit 50 mm (wells 0 50 200 100 100 200 50 mm), wells 75 mm, period 525 mm, fins 1.6 mm, 6 periods =
3.15 m, normal incidence. The energy of lobe q is taken between the half-way points to its neighbours and turned into
|R_q|² (BEM: divided by cos θ_q, because a lobe at an angle carries power ~ |R_q|² cos θ_q; the Fourier model has no
obliquity factor). Lobes at grazing (|sin θ| > 0.95) are left out.

| Frequency | Lobes | Spread of \|R_q\|², BEM | Spread, Fourier | Specular share, BEM / Fourier |
|---|---|---|---|---|
| f0 = 490 Hz | specular only (λ0 = 700 mm > period) | - | - | 100 % / 100 % |
| 2 f0 = 980 Hz | 3 | 3.7 dB | 0.4 dB | 18 % / 34 % |
| 3 f0 = 1470 Hz | 5 | 3.1 dB | 1.3 dB | 24 % / 23 % |
| 4 f0 = 1960 Hz | 5 | 4.8 dB | 1.3 dB | 22 % / 23 % |
| 5 f0 = 2450 Hz | 7 | 12.9 dB (lobes ±3 at 53° 10 dB weaker) | 2.7 dB | 25 % / 18 % |
| 6 f0 = 2940 Hz | 9 | 7.4 dB | 5.1 dB | 17 % / 16 % |
| 7 f0 = 3430 Hz (flat plate) | 9 | side lobes 19-29 dB below specular | 19-33 dB | 97 % / 97 % |

- The known properties hold: at 2-4 f0 the lobes carry similar energy (BEM between -3 and +1.7 dB around the mean, the
  Fourier model more even), and at 7 f0 the QRD reflects like a flat panel (97 % of the energy specular in both).
- From 5 f0 BEM moves away from equal lobes: a 73.4 mm well is wider than half a wavelength above
  c / (2 · 73.4 mm) = 2337 Hz, so cross modes propagate in the wells, which the plane-wave-in-a-well Fourier model does
  not know.
- The largest difference is at f0. The periodic QRD reflects exactly like a flat panel in BEM (peak -0.04 dB re the
  flat panel), while the Fourier model shows -8.6 dB. The Fourier model loses 7.9 dB of the total energy there: at f0
  the only real direction is the specular one (period 525 mm < λ0 700 mm), and the rest goes into evanescent waves that
  this model simply drops. In a lossless rigid body that energy does not vanish, it comes back in the specular
  direction (with real losses part of it becomes absorption). At 1.5 f0 (735 Hz): BEM -3.5 dB, Fourier -12.8 dB.

![Standard QRD of Cox, Angus and D'Antonio, BEM and Fourier](../checks/bem/out/fig-validation-cox-qrd.png)

## Results for the tiles and walls

19 configurations (flat boxes of 1-7 tiles, tile 1, walls in the published order and of copies of tile 1), 113
frequencies each, 244-4947 elements, 11 minutes in 8 processes. Raw far fields: `checks/bem/out/raw/*.npz` (4 MB),
numbers: `checks/bem/out/results.json`, tables below from [`checks/bem/tables.py`](../checks/bem/tables.py).

Ranges: c/W = c / wall width, c / period = 2 f0. Means over the third-octave bands whose centre lies in the range.
Bands where the flat reference itself scatters almost evenly (d_flat > 0.9: a small tile at 0.8-1.25 kHz) are left out.

#### A1 mini: f0 1114 Hz, 2 f0 = c / period 2227 Hz, 7 f0 7795 Hz

d_n, mean of incidences 0, ±30, ±60° (the optimiser's score), BEM / Fourier:

| Wall | c/W | below c/W | c/W - f0 | f0 - 2 f0 | 2 f0 - 7 f0 | 7 f0 band | above 7 f0 |
|---|---|---|---|---|---|---|---|
| 1 tile | 2227 Hz | -0.23 / 0.00 | - | 0.21 / 0.83 | 0.52 / 0.75 | 0.25 / 0.09 | 0.64 / 0.60 |
| 4 tiles, wall order | 557 Hz | 0.04 / 0.01 | 0.09 / 0.19 | 0.44 / 0.62 | 0.37 / 0.54 | 0.05 / 0.02 | 0.48 / 0.43 |
| 4 tiles, copies of tile 1 | 557 Hz | 0.03 / 0.00 | 0.01 / 0.01 | 0.21 / 0.25 | 0.25 / 0.33 | 0.06 / 0.02 | 0.36 / 0.30 |
| 7 tiles, wall order | 318 Hz | 0.04 / 0.00 | 0.10 / 0.05 | 0.47 / 0.48 | 0.28 / 0.39 | 0.03 / 0.01 | 0.43 / 0.34 |
| 7 tiles, copies of tile 1 | 318 Hz | 0.04 / 0.00 | 0.04 / 0.00 | 0.14 / 0.17 | 0.16 / 0.23 | 0.03 / 0.01 | 0.22 / 0.20 |

Lowest third-octave band with d_n >= 0.3 in two bands in a row, BEM / Fourier:

| Wall | c/W | mean of 0, ±30, ±60° | 0° | 45° | -45° |
|---|---|---|---|---|---|
| 1 tile | 2227 Hz | 2000 / 1250 | 1600 / 1600 | 1000 / 1000 | 1000 / 1000 |
| 4 tiles, wall order | 557 Hz | 1600 / 1000 | 1600 / 1000 | 1600 / 1250 | 1600 / 1000 |
| 4 tiles, copies of tile 1 | 557 Hz | never / 2500 | never / 2000 | never / 1250 | never / 1250 |
| 7 tiles, wall order | 318 Hz | 1000 / 1250 | 1250 / 1250 | 1250 / 1250 | 1600 / 1250 |
| 7 tiles, copies of tile 1 | 318 Hz | never / never | never / never | never / never | never / never |

One tile between f0 and 2 f0 (BEM / Fourier): 0° 0.68 / 0.96, 45° -0.19 / 0.79.

#### P1S / X1: f0 742 Hz, 2 f0 = c / period 1485 Hz, 7 f0 5197 Hz

d_n, mean of incidences 0, ±30, ±60° (the optimiser's score), BEM / Fourier:

| Wall | c/W | below c/W | c/W - f0 | f0 - 2 f0 | 2 f0 - 7 f0 | 7 f0 band | above 7 f0 |
|---|---|---|---|---|---|---|---|
| 1 tile | 1485 Hz | -0.12 / -0.05 | - | 0.07 / 0.78 | 0.49 / 0.75 | 0.13 / 0.09 | 0.59 / 0.41 |
| 3 tiles, wall order | 495 Hz | 0.02 / 0.03 | 0.09 / 0.24 | 0.48 / 0.78 | 0.42 / 0.60 | 0.03 / 0.03 | 0.43 / 0.28 |
| 3 tiles, copies of tile 1 | 495 Hz | 0.05 / 0.00 | 0.00 / 0.01 | 0.26 / 0.30 | 0.31 / 0.40 | 0.03 / 0.02 | 0.38 / 0.22 |
| 5 tiles, wall order | 297 Hz | 0.08 / 0.01 | 0.08 / 0.07 | 0.43 / 0.58 | 0.37 / 0.51 | 0.02 / 0.01 | 0.31 / 0.19 |
| 5 tiles, copies of tile 1 | 297 Hz | 0.06 / 0.00 | 0.03 / 0.00 | 0.17 / 0.19 | 0.22 / 0.29 | 0.02 / 0.01 | 0.25 / 0.14 |
| 7 tiles, wall order | 212 Hz | 0.03 / 0.00 | 0.07 / 0.03 | 0.45 / 0.44 | 0.30 / 0.40 | 0.01 / 0.01 | 0.24 / 0.13 |
| 7 tiles, copies of tile 1 | 212 Hz | 0.03 / 0.00 | 0.04 / 0.00 | 0.13 / 0.16 | 0.17 / 0.24 | 0.01 / 0.01 | 0.18 / 0.10 |

Lowest third-octave band with d_n >= 0.3 in two bands in a row, BEM / Fourier:

| Wall | c/W | mean of 0, ±30, ±60° | 0° | 45° | -45° |
|---|---|---|---|---|---|
| 1 tile | 1485 Hz | 1600 / 800 | 1250 / 1000 | 630 / 800 | 630 / 800 |
| 3 tiles, wall order | 495 Hz | 1000 / 630 | 1000 / 630 | 1000 / 800 | 1000 / 630 |
| 3 tiles, copies of tile 1 | 495 Hz | 800 / 1250 | never / 1250 | 800 / 800 | 800 / 800 |
| 5 tiles, wall order | 297 Hz | 1000 / 800 | 1000 / 800 | 1000 / 800 | 1000 / 800 |
| 5 tiles, copies of tile 1 | 297 Hz | never / 1600 | never / never | 800 / never | 800 / never |
| 7 tiles, wall order | 212 Hz | 800 / 800 | 1000 / 800 | 800 / 800 | 1000 / 800 |
| 7 tiles, copies of tile 1 | 212 Hz | never / never | never / never | never / never | never / never |

One tile between f0 and 2 f0 (BEM / Fourier): 0° 0.50 / 0.89, 45° -0.23 / 0.77. At 45° the P1S tile has d_n >= 0.3
in the 630 and 800 Hz bands and then drops below zero at 1-1.25 kHz, so "630" in that row is not the start of steady
diffusion.

Every wall and the single tile have their own figure with d_n at 0° and 45° and polar plots at 200 Hz, 500 Hz, f0,
2 f0 and 7 f0: `checks/bem/out/fig-A1mini-tile.png`, `fig-A1mini-wall-4.png`, `fig-A1mini-wall-7.png`,
`fig-P1S-tile.png`, `fig-P1S-wall-3.png`, `fig-P1S-wall-5.png`, `fig-P1S-wall-7.png`. For example the widest wall:

![P1S / X1, 7 tiles: d_n at 0° and 45° and polar plots](../checks/bem/out/fig-P1S-wall-7.png)

## Where the Fourier model agrees with BEM, and where it does not

It agrees on:
- The ranking. The order picked by the optimiser beats copies of one tile in BEM between f0 and 2 f0, between 2 f0 and
  7 f0 and above 7 f0, for every wall and both printers (for example P1S, 7 tiles, f0 - 2 f0: 0.45 against 0.13); in
  the 7 f0 band both are a flat panel. So the Fourier model picked an order that is also better than copies in the
  physical model. For the 3-tile P1S wall it is in the top fifth of all orders in BEM, but not the best (next
  section).
- The 7 f0 band: d_n 0.01-0.06 for every wall in both models. This is the flat-plate frequency, and the band ends there.
- Below c / wall width, and between c / wall width and f0 for walls of 5-7 tiles: almost zero in both.
- 7-tile walls between f0 and 2 f0: nearly the same means (P1S 0.45 / 0.44, A1 mini 0.47 / 0.48), band by band within
  0.1 (0.2 at normal incidence).
- Polar plots at 7 f0: the same single specular lobe.

It does not agree on:
- Anything with the period of one tile below c / period (2 f0): the single tile and the copies. The Fourier model loses
  the energy of the evanescent waves there and shows a weaker specular reflection than a flat panel. Lossless BEM sends
  that energy back in the specular direction. P1S tile, f0 - 2 f0: BEM 0.07, Fourier 0.78. It is the same difference
  that gives -0.04 dB (BEM) against -8.6 dB (Fourier) for the standard QRD at f0.
- Small walls (3 P1S tiles, 4 A1 mini tiles) between f0 and 2 f0 and below f0: the Fourier model is too high (P1S,
  3 tiles: 0.78 against 0.48 and 0.24 against 0.09). At 69 cm a few tiles are still almost one period.
- Between 2 f0 and 7 f0 the Fourier model is too high: walls in the published order by 0.10-0.18 on average and by
  0.22-0.36 at normal incidence (P1S, 7 tiles, 0°: 0.24 against 0.46), copies by 0.07-0.09. BEM has deep dips just
  above c / period at normal incidence (P1S 7 tiles, 1.6-2 kHz: 0.07-0.09 against 0.34-0.59; A1 mini 7 tiles,
  2.5-3.15 kHz: 0.05-0.13 against 0.49-0.50).
- Above 7 f0 BEM gives more than the Fourier model at 45° (P1S, 7 tiles: 0.34 against 0.13) and less at 0°. The wells
  are wider than half a wavelength there, and a model with a plane wave in each well does not see the cross modes.
- The Fourier model has no obliquity factor (cos θ), which puts too much energy at wide receiver angles.

## Do the cyclic shifts make a wall diffuse from c / wall width?

Not from c / wall width, but from about f0, and only for wider walls.

- Between c / wall width and f0 a wall in the published order has d_n 0.07-0.10 in BEM, copies 0.00-0.04. There is a
  difference, but the wall still reflects like a panel. The phase differences between the wells are too small to
  spread the reflection (both models give almost zero for walls of 5-7 tiles), so the shifts have nothing to rearrange.
- Between f0 and 2 f0, where one tile and copies of one tile reflect specularly, a wall of cyclic shifts really
  diffuses (BEM 0.43-0.48). That is what the shifts buy: without them this octave would be lost, because the tile
  period (λ0 / 2) is too short for grating lobes at f0.
- Start of diffusion: P1S, 7 tiles (162 cm), from the 800 Hz band; 5 and 3 tiles from 1000 Hz; A1 mini, 7 tiles
  (108 cm), from 1000 Hz; 4 tiles (62 cm) from 1600 Hz. The Fourier model gave the same or one band higher for 7 tiles
  (A1 mini 1250 Hz) and one to two bands lower for smaller walls (P1S, 3 tiles, 630 Hz; A1 mini, 4 tiles, 1000 Hz).

That is where the generator's √2 · f0 start for deep walls narrower than 3.5 design wavelengths comes from: √2 · f0 is
1050 Hz on the P1S / X1 and 1575 Hz on the A1 mini, next to the 1000 Hz and 1600 Hz bands found here.

## Is the published order the best one? Every order of a 3-tile P1S / X1 wall

[`checks/bem/order_ranking.py`](../checks/bem/order_ranking.py) scores all 7 · 6 · 5 = 210 orders of three P1S / X1
tiles with the optimiser's own score: mean d_n over the incidences 0, ±30 and ±60° and the third-octave bands from f0
to below 7 f0 (800 Hz - 4 kHz), 3 frequencies per band, against the flat box of the same width. A tile turned upside
down is tile 9 - t (tile 1 stays), so the order a b c in a mirror is m(c) m(b) m(a), the same wall mirrored, with the
same score; BEM ran one order of every mirror pair (108 pairs) plus both sides of three pairs, which agree to five
decimals. The Fourier model ran all 210. About 8 minutes in 8 processes. Summary: `checks/bem/out/ranking-P1S-3.json`,
band energies of every BEM run: `checks/bem/out/ranking-P1S-3-bem.jsonl`.

| | BEM | Fourier |
|---|---|---|
| Score of all orders (min / median / max) | 0.26 / 0.39 / 0.45 | 0.42 / 0.54 / 0.66 |
| Published order `6 3 4` (mirror `5 6 3`) | 0.41, rank 22 of 108 | 0.66, rank 1 |
| Best order | `5 6 2` (mirror `7 3 4`), 0.45 | `6 3 4`, 0.66 |
| BEM best in the Fourier ranking | rank 6 (0.61) | |
| Between f0 and 2 f0 only (800-1250 Hz): published order | 0.43, rank 28 | 0.77, rank 1 (next best 0.66) |
| Between f0 and 2 f0 only: best order | `2 1 4` (mirror `5 1 7`), 0.54 | `6 3 4`, 0.77 |

- The rankings agree moderately: Spearman rank correlation 0.68 (Kendall 0.49) over the whole score, 0.52 between f0
  and 2 f0. Only 2 of the BEM top 10 are in the Fourier top 10. Both models put the same orders at the bottom
  (`1 2 3`, `3 4 5` and `2 1 7` are the last three in BEM and between 101st and 107th in the Fourier model).
- The published order is a good one: in the top fifth in BEM, 0.15 above the worst order and 0.04 below the best.
  It is not the best, and the reason is visible in the right half of the figure: between f0 and 2 f0 the Fourier model
  rates it far above every other order (0.77, the next one 0.66), while BEM puts it 28th. The optimiser picked the order
  where the Fourier model is most wrong in the range where it is least reliable.
- The difference to the BEM best (0.04 in d_n) is ten times the mesh error of the solver (0.004) but small against what
  a rigid 2D model leaves out (losses, the finite tile length).
- What it means for the generator: the Fourier model is fine for ruling out bad orders and a poor judge among the good
  ones. A better order would come from BEM, for example on a Fourier shortlist. The 5- and 7-tile walls were not ranked
  (2520 and 5040 orders at about 1.5 and 2.5 minutes each in one process, so about 4 and 14 hours in 8 processes even
  with mirror pairs merged), and changing the order would renumber every published wall plan, so nothing in the
  generator changes because of this.

![All orders of a 3-tile P1S wall, BEM against Fourier](../checks/bem/out/fig-P1S-3-ranking.png)

Each dot is one order and its mirror image; orange the published order, the black ring the best order in BEM for the
score of that panel.

## Still open

- The best order of the 5- and 7-tile walls in BEM (a shortlist from the Fourier model, then BEM), and whether the
  difference is worth renumbering the walls.
- Block tiles were not simulated (they need 3D BEM).
- A measurement is still the only proof. What to look for: the specular reflection of a wall in the published order
  dropping from about f0, no such drop for copies at normal incidence, the specular reflection coming back at 7 f0
  (P1S 5.2 kHz, A1 mini 7.8 kHz).

## Running it

Python 3.11 or newer with `pip install -r checks/bem/requirements.txt` (numpy, scipy, matplotlib), from the repository
root:

```
python checks/bem/validate.py     # solver checks, about 5 minutes -> checks/bem/out/validation.json, fig-validation-*.png
python checks/bem/sweep.py        # the BEM runs, about 11 minutes with 8 processes -> checks/bem/out/raw/*.npz
python checks/bem/analysis.py     # diffusion coefficients and figures -> checks/bem/out/results.json, fig-*.png
python checks/bem/tables.py       # the tables of this page
python checks/bem/order_ranking.py   # all orders of a 3-tile P1S wall, about 8 minutes -> checks/bem/out/ranking-P1S-3.*
```

`sweep.py` skips configurations that already have a raw file, so delete `checks/bem/out/raw/` to run everything again.

## References

- T. J. Cox, P. D'Antonio, *Acoustic Absorbers and Diffusers: Theory, Design and Application*, 3rd ed., CRC Press, 2016
  (chapters on Schroeder diffusers and on prediction models).
- T. J. Cox, J. A. S. Angus, P. D'Antonio, "The performance of Schroeder diffusers: a wide band BEM investigation",
  Proc. IOA 25(5), 2003 (the standard N = 7 QRD).
- A. J. Burton, G. F. Miller, "The application of integral equation methods to the numerical solution of some exterior
  boundary-value problems", Proc. R. Soc. Lond. A 323, 1971, 201-210.
- AES-4id-2001 and ISO 17497-2:2012 for the diffusion coefficient.
