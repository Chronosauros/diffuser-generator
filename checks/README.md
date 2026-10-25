# Checks

Independent checks of the maths and of the model behind the generator, made after the review on r/acoustics
(24-25.10.2026, with AI assistance). They are not part of `npm test`: they are reports to read, and each one prints
or writes what it found. Everything is a calculation or a simulation; nothing here is a measurement.

## Published recipes through the generator core: `recipes.mjs`

```
node checks/recipes.mjs [path/to/index.html]
```

Loads `<script id="diffuser-core">` from `index.html` (as `npm test` does) and rebuilds diffusers published elsewhere:
the standard N = 7 QRD of Cox, Angus and D'Antonio (Proc. IOA 25(5), 2003: f0 490 Hz, depth unit 50 mm, wells 75 mm,
period 525 mm), the residue tables and build depths of the QRDude calculator (subwoofer-builder.com/qrd.htm, N = 7,
11, 13, and the 2D (n² + m²) mod N pattern), a worked example from the Gearspace forum (28 in period, depths
0 2 8 4 4 8 2 in) and the soundassured.com calculator. Result: all twelve checks pass. The soundassured calculator
uses another convention (deepest well = λ/2), so its "1000 Hz" panel is a 571 Hz panel in Cox and D'Antonio's terms.

It also prints the default tiles against the QRDude rules of thumb. Two of them are not met, and both are design facts,
not errors: a tile is one period of λ0 / 2 (QRDude wants a period of at least λ0), which is why a single tile only
diffuses from c / tile width = 2 f0; and the A1 mini wells are 20.4 mm (19.2 mm at the edges), below the 25 mm minimum
QRDude gives for viscous losses. The live page cannot rebuild the standard QRD exactly: N is fixed at 7 and the well
pitch comes from the printer bed.

## The Fourier model against Schroeder theory: `fourier_lobes.py`

```
python checks/fourier_lobes.py [path/to/reference]
```

Needs numpy. Evaluates the far-field model of `reference/layout.py` (the one that picked the wall orders) at the
grating-lobe directions of a periodic QRD. Theory: at m · f0 every lobe carries the same energy apart from the
well-width envelope, and at 7 f0 only the specular lobe is left. Result: for the standard QRD with zero-thickness fins
the lobe energies after the envelope correction are equal to four decimals at 2-6 f0 and the side lobes at 7 f0 are
1e-30 of the specular one, so the model is implemented as the theory says. For our tiles the 1.6 mm fins and the
narrower edge wells spread the lobe energies by 4-18 % (up to 0.7 dB). This checks the implementation, not whether the
model is right about real diffusers; that is the next check.

## 2D BEM against the Fourier model: `bem/`

A 2D boundary element simulation of the deep tiles and walls, with its own validation, compared band by band with the
Fourier model. Write-up, tables and figures: [docs/bem-check.md](../docs/bem-check.md). Code and data:

- `bem/bem2d.py` - the solver (Burton-Miller, constant elements, exact near-singular integrals);
- `bem/geometry.py` - cross-sections of tiles, walls and flat references from `reference/qrd.py`, and the Fourier model
  of `reference/layout.py` on any frequencies;
- `bem/validate.py` - solver checks (cylinder, interior resonance, Kirchhoff, standard QRD, mesh convergence);
- `bem/sweep.py`, `bem/analysis.py`, `bem/tables.py` - the runs, the diffusion coefficients and figures, the tables;
- `bem/order_ranking.py` - every order of a 3-tile P1S / X1 wall in BEM and in the Fourier model, ranked;
- `bem/order_search.py` - a BEM search for a better 7-tile P1S / X1 order over the 3-, 5- and 7-tile walls together;
- `bem/order_walls.py` - per-wall BEM rankings: P1S / X1 5 and 7 tiles, A1 mini 4 and 7 tiles;
- `bem/order_nested.py` - the best chains of nested walls from every cached BEM run, and a check of the orders in
  `reference/out/layouts.json` (cached numbers and a fresh BEM run of the smallest wall of each class). These four
  scripts were made against the wall order published until 25.10.2026 and keep it as a named constant, so their numbers
  do not change when `layouts.json` does; each can rebuild its summary from the cached band energies without any BEM
  (`--cache-only`, `--no-bem` for order_nested.py);
- `bem/out/` - `validation.json`, `results.json`, `raw/*.npz` (far fields of every run), the band energies of the order
  runs (`*-bem.jsonl`, 9 MB), their summaries (`ranking-P1S-3.json`, `search-P1S.json`, `walls.json`, `nested.json`) and
  the figures.

`pip install -r checks/bem/requirements.txt` (numpy, scipy, matplotlib), then the commands in docs/bem-check.md.
