# Bass traps: the maths, the assumptions and the honest effect

This page describes the model behind the bass trap mode of the generator (the cores `bass-physics` and `bass-geometry`
in `index.html`, the Python twin in `reference/bass_room.py` and `reference/bass_cell.py`). Every number it produces is
**calculated, not measured**. No trap was printed, placed in a room or measured; the model has not been checked
against a measurement of a real room. The research it ports was done on 25.10.2026 with a transfer-matrix model of the cell
and a 3D finite element model of a room with a patch of cells, both written for this project.

The generator, the reference implementation and this page were written with AI assistance (Claude). If something here is
wrong, please open an issue.

## What the traps are meant to do

A room mode is a standing wave between the walls. Below roughly 120 Hz a small room has a handful of well separated modes,
and each one rings for a time that depends on how much the walls absorb at that frequency. A bass trap here is a **narrow
resonant absorber tuned to one mode**: a folded quarter-wave channel (QWR) in a cell 120 x 120 mm, a slit of 6 mm wound
into a depth of 50 - 150 mm. N cells together sit in the place where the mode has its pressure maximum (a corner, an edge or
the middle of a wall) and take energy out of that mode.

It is not a broadband absorber. It does nothing for the other modes, for the treble or for flutter echo. For broadband bass, porous or membrane
traps in the corners do more per kilogram of material; this is a narrow tool for one ringing mode. The honest numbers
are small (see "What it does and does not do"); the page is written so that you can judge them yourself.

## The model chain

Notation: c = 343.2 m/s, Z0 = rho0 c, exp(+j w t) so a decaying mode has Im(w) > 0. Lengths in mm unless a name ends in _m.
Decay rates are for the **pressure amplitude** (1/s); the decay of the energy is twice that. T60 = ln(1000) / delta.

### 1. Room modes

Rigid walls: psi_n(r) = cos(nx pi x / Lx) cos(ny pi y / Ly) cos(nz pi z / Lz), f_n = (c / 2) sqrt((nx / Lx)^2 + (ny / Ly)^2 +
(nz / Lz)^2), norm Lambda_n = V / 2^k with k the number of nonzero indices. All modes up to 120 Hz are listed; axial modes have one
nonzero index, tangential two, oblique three.

### 2. Base decay from the walls (an assumption: the absorption of the walls)

The walls are given one Sabine coefficient alpha (or one per wall), taken as a local-reaction wall with Re beta = alpha / 8
(small beta, diffuse field). Green's identity for the pressure and the rigid mode gives, to first order in beta,

delta_n = (c / (2 Lambda_n)) * integral of Re(beta) psi_n^2 dS.

For the axial mode (1,0,0) of the example room 4.0 x 3.2 x 2.6 m at alpha = 0.10 this is c beta (2 / Lx + 1 / Ly + 1 / Lz) =
5.136 1/s, T60 = 1.345 s. **alpha is a guess for your room** (no furniture, doors, windows, no frequency dependence); the
factor of about 2 between definitions of alpha alone moves every result, see "Limits".

### 3. The cell: impedance at complex frequency

The channel is a slit of width w = 6 mm. Its effective density and bulk modulus include the viscous and thermal boundary
layers (parallel plates): rho_eff = rho0 / Xi, K_eff = gamma P0 / (gamma - (gamma - 1) Psi) with Xi = 1 - tanh(G) / G, G = (w / 2)
sqrt(j w rho0 / mu), and the same with the Prandtl number for Psi. The input impedance of the closed channel of acoustic length
L is Zc / (j tan(k L)) (with Zc = sqrt(rho_eff K_eff) / w, k = w sqrt(rho_eff / K_eff)), plus the mass of the end correction of the
slit in an array with the cell width as period (4000 terms), scaled to the face: Zs = Zin * 120 mm. It is evaluated at **complex**
w, which is the point of the next step.

The channel length for a resonance f is found as in the research (`design_Q.first_res` / `solve_L`): the first upward zero of
Im Zs on the grid 30 - 520 Hz in steps of 0.25 Hz with linear interpolation, and 40 bisection steps in L. The JS and the Python
code use the same grid, so they give the same L (the tests compare them to 1e-9).

### 4. The trap changes the characteristic equation of the mode

N cells on the walls, each with the surface impedance Zs(w), couple to the mode with the strength G = sum over cells of
A <psi>_cell^2 / Lambda (A = 0.12 m x 0.12 m, <psi> averaged over the cell face). The cells replace the wall under them, so the
wall part of the decay becomes delta' = delta_n - (c / 2) sum A <psi>^2 / Lambda * Re beta_wall. The mode is then the root of

(w^2 - wn^2 - 2 j w delta') Zs(w) = j w c Z0 G

near wn. Solved exactly (the cubic with Zs linearised around wn as the start, polished by Newton's method on the exact Zs),
this gives two roots: the mode and the cell, pushed apart. The **slower** of the two roots near the mode is the decay that
remains: delta_eff. The classic perturbation formula (Re beta on the real axis) underestimates the effect by about
gamma / (gamma - delta), because the width of the cell resonance gamma is only 2 - 3 times the width of the mode delta; the
model evaluates the cell where the mode really is, in the complex plane.

T60 after = ln(1000) / delta_eff. The **decay factor** is delta_eff / delta_n, which is also T60 before / T60 after (a factor of
1.5 means a T60 one and a half times shorter).

### 5. The radiation mass of the patch

A patch of N cells in a wall is a piston of radius a = sqrt(N A / pi); its reactance loads the cell like a mass of
rho0 * 8 a / (3 pi) per m^2 (ka << 1), which pulls the resonance down. Images in the walls make the piston larger: the area grows by
the number of images, the radius by its square root, so the mass is multiplied by 2 in a corner, sqrt(2) on an edge and 1 on a wall.
This mass is part of Zs in step 4 (the research package: without it the effect is overstated by 2 - 4 times). The cells are therefore
designed for **f cell = f mode x 1.005**.

### 6. The seven anchors and the best place

For a mode the model is run for seven places: the floor corner, the vertical corner edge at half height, the two floor edges, and
the middle of each of the three wall pairs. For a given N the cells are the N cells closest to the anchor on the 0.12 m grid of
the adjoining walls. The **best place for N** is the anchor with the largest delta_eff, described in words (for example "middle of
the 3.2 x 2.6 m wall", "vertical corner edge, half height (2.6 m)", "floor corner"). A place at a node of the mode (psi^2 of the
first cell below 0.05) is marked nodal and gives nothing.

### 7. The drop of the peak in dB (new in the generator)

For a force on the mode the response is H(w) = 1 / (w^2 - wn^2 - 2 j w delta_n) without the trap and
H(w) = Zs / ((w^2 - wn^2 - 2 j w delta') Zs - j w c Z0 G) with it. The **peak drop** is 20 log10 (max |H_trap| / max |H_base|) on
w = wn x [0.85 ... 1.15], 3001 points; negative means a lower resonance peak. It is a check on the same model, not a tuning. For 8
cells at 42.9 Hz and alpha 0.10 (wall beta 0.0125) the finite element model gave -1.3 dB (the conservative figure). The -1.7 dB sometimes quoted is
for wall beta 0.01, about alpha 0.08 (N = 1 / 8 / 16: -0.23 / -1.7 / -3.0 dB there); at alpha 0.10 only -1.3 dB holds.

## The geometry

Cell face 120 x 120 mm, walls and floors 1.2 mm (three perimeters of a 0.4 mm nozzle), slit **6 mm** (the only width: 8 mm
works in some cases and fails in others depending on a resistance nobody has measured, see `slot-width-decision` in the research),
depth 50 - 150 mm from the number of legs. The channel is extruded along z, which is the print axis, so there are no supports.
Every 90 degree bend shortens the channel acoustically by 0.245 w = 1.47 mm (2D finite element model of a bend, 19 geometries),
so the axis is made longer by that amount per bend. The inlet is a slit in the front wall, the last leg ends in a window in a side
wall, closed by the **nominal insert (k = 0), which is the end cap and is part of every cell**. The other inserts, -1 % to +4 %
in steps of 0.5 % of the frequency (the body is longer for a higher frequency), are optional tuning parts for use **only with a
decay measurement before and after**: without a measurement there is no reason to pick any other than the nominal one.

Cells are packed in blocks of NY x NZ cells (NY across the face width, NZ stacked), filled in rows from the bottom; the layouts tried,
in this order, are 2x2 with 4 cells, 2x2 with 3 (an L, no overhang), 1x2, 2x1 and 1x1. A block is a cuboid and prints without
supports; the generator leaves 5 mm free on every bed axis and takes the first layout that fits the bed and prints in at most 12 h
(estimate). If one cell does not fit (deeper than the bed or than 150 mm, or higher than the bed) the state is "does not fit",
with the reasons, and nothing else is offered.

### Print time and filament

The estimate is a least squares fit to the 10 blocks and 10 inserts plates that the research sliced (Bambu Studio, P1S, 0.4 mm
nozzle, PLA, process "Diffuser blocks 0.20", 3 walls, 15 % infill): block minutes = 1.393 volume + 0.0648 footprint, block grams =
1.175 volume + 0.0493 footprint (volume in cm3, footprint in cm2); the plate uses the surface area, volume and height. Held-out
errors (leave one mode out): 0.15 % - 0.55 % (blocks), 0.65 % - 1.4 % (plates). With ten rows per kind it is a rough estimate for
shapes like the research ones and **only for that printer and profile**; the numbers are shown as approximate.

## Flags

- **No effect**: delta_n > gamma of the cell (Re Zs / |dZs/dw| at the cell frequency). The mode decays faster than the cell can
  answer; the trap then does not shorten the decay (the fixtures: the decay is a little *slower*) although the peak still drops a
  little. This replaces the research's first T60 < 1 s test, which was wrong for oblique modes.
- **Approximate**: the single-mode model leaves the multi-mode results by 25 % or more (up to about 40 % against the 3D simulation) when modes overlap (the nearest mode closer
  than the sum of the two mode widths), for more than 8 cells, for tangential and oblique modes with 4 or more cells (the research
  rule `orient`), and for a corner with a neighbouring mode within 10 % of the cell frequency. The reasons are listed.
- **Nodal**: the first cell sits where psi^2 < 0.05.
- **Uncertain** (reasons listed): delta between 0.8 gamma and gamma; an oblique mode at alpha >= 0.10; any alpha >= 0.15, where
  the sign of the effect depends on the unmeasured resistance of the channel.

## What it does and does not do

The figures below are this model's (the cells tuned 0.5 % above the mode, as the generator makes them) next to the finite element (FEM,
3D simulation) figures where the research has them, for alpha = 0.10, the example room and the mode (1,0,0) at 42.9 Hz. They are
**calculated, not measured**.

| Cells N | Best place (this model) | T60 before -> after (this model) | Decay factor, this model | Decay factor, FEM | Peak, this model | Peak, FEM |
|---|---|---|---|---|---|---|
| 1 | floor corner | 1.345 -> 1.294 s | x1.04 | x1.04 | -0.18 dB | -0.2 dB |
| 2 | floor corner | 1.345 -> 1.242 s | x1.08 | - | -0.36 dB | - |
| 4 | vertical corner edge | 1.345 -> 1.139 s | x1.18 | - | -0.71 dB | - |
| 8 | middle of the 3.2 x 2.6 m wall | 1.345 -> 0.868 s | x1.55 | x1.37 tuned on the mode, x1.50 tuned +0.5 % | -1.36 dB | -1.3 dB |
| 16 | middle of the 3.2 x 2.6 m wall | 1.345 -> 0.842 s | x1.60 | - | -2.55 dB | - |

The 8-cell FEM figure depends on the tuning: x1.37 is for a cell tuned exactly on the mode, x1.50 for +0.5 % (what the generator prints).

- **One cell is practically inaudible**: T60 1.35 -> 1.30 s, the peak 0.2 dB lower. It is below the threshold of audibility and
  nothing a listener could hear.
- **Eight cells** on an axial mode at 42.9 Hz shorten the decay of that mode by x1.37 (tuned on the mode) to x1.50 (tuned +0.5 %) in the finite element model (this
  model: x1.38 in a corner to x1.55 on the wall) and lower the peak by 1.3 dB (1.7 dB only for wall beta 0.01, about alpha 0.08). That is 8 prints of about 8 h and 3.4 kg of PLA
  for one mode.
- **There is a ceiling**: about x1.5 - 1.75 however many cells (finite element model, beta = 0.01, up to 64 cells; the single-mode
  model gives x1.6 at 16). More cells make the patch larger, which pulls its resonance away from the mode.
- **When the mode already decays faster than the cell resonance is wide (delta > gamma), the trap does nothing** to the decay.
  At alpha = 0.20 the (1,0,0) mode of the example room has delta = 10.3 1/s against gamma = 11.7: eight cells give x1.05 and
  the result is flagged uncertain.
- Everything scales with the unknown absorption of your walls: at alpha 0.05 the same cells give x1.06 (one) and x1.50 (eight), and at
  beta >= 0.03 (alpha about 0.24) nothing (finite element model).
- It is a narrow band absorber for one mode. The neighbouring modes get a little if they are within a few Hz.

## Limits of the model (read before building anything)

1. **The resistance of the channel is the largest unknown, plus or minus a factor of 2.** FDM walls are rough and not perfectly
   airtight; the model uses a smooth channel. The size of the effect changes by 0.3 - 3 times with it, and for wider slits its sign.
2. **The tuning is narrow.** A detuning of +-1 % costs 10 - 40 % of the effect (FEM, 8 cells: +1 % gives x1.33 instead of x1.50), and the
   frequency of a printed cell is uncertain by 1 - 6 % (below).
   **The resonance of a printed cell is lower than calculated**: roughness lowers it by about 2.5 % (a viscosity 4 times higher) to
   6 % (25 times); a leak at the flange of the insert raises it. The tuning range of the inserts (-1 % ... +4 %) does not cover all of
   that. The inserts help **only with a measurement of the decay before and after** (sweep and decay plot with any measurement
   microphone); without it use the nominal insert.
3. **The seal of the insert flange is not tested.** Foam tape and aluminium tape are suggested; the end of the channel is a
   pressure maximum, so a leak shortens the resonator.
4. **One mode at a time.** The single-mode model is exact for isolated modes (it agrees with the multi-mode model to 0.0 % on the
   example room) and differs by 25 % or more (up to about 40 % against the 3D simulation) when modes overlap or lie within 10 % of the cell; those results are flagged.
   Cells cut into the next modes only through the flags, not through numbers.
5. **The walls are uniform** (one alpha, local reaction, no frequency dependence). A real room has windows, doors, furniture and
   openings that give each mode a different base decay, and the result depends on delta / gamma.
6. **The patch is a circular piston** for the radiation mass; this is weak for long thin patches (a column along a corner) and for
   a pair of modes. Against the finite element model the research model agreed in the decay of the mode in 12 of 12 FEM comparison cases and in the
   gain within 20 % in 9 of 12 (the other three: +26 %, -41 % and +34 %; the review of 25.10.2026).
7. **Calculated, not measured**: neither the cell (the transfer-matrix model was not validated on a printed part) nor the room (the
   finite element model was not compared to a measurement) has a measurement behind it.
8. The cells are a research result. The printed cell, the room and the effect should be checked with a measurement before any
   claim is made.

## Code checks (not acoustic tests)

These are checks that the code does what the research code does. None of them measures a room or a printed cell.

- `scripts/bass-physics.test.mjs` compares the JS core with the Python twin on 5 rooms (4.0 x 3.2 x 2.6, 3.0 x 2.5 x 2.4,
  6.0 x 4.5 x 2.8, 7.0 x 3.0 x 2.5, 3.0 x 2.9 x 2.7 m), alpha 0.05 / 0.10 / 0.20, every mode up to 120 Hz, seven anchors and
  N = 1, 2, 4, 8, 16 (relative 1e-6, the roots 1e-9), and checks the flags on cases known from the research.
- `scripts/bass-geometry.test.mjs` compares the rectangle plans (1e-6 mm), the solids and binary STL files (volume 0.1 %, bounding
  box 0.01 mm, sections by ray cast, one closed solid with the pockets as the only cavities), the packing on the beds and the
  estimate. The STL of the research for 42.9 Hz and 95.2 Hz equal the regenerated reference files.
- `python scripts/bass_reference.py physics|geometry` regenerates `scripts/fixtures/bass/`; the fit of the estimate is
  `scripts/calibration/bass-calibrate.py`.
