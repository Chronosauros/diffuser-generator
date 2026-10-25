# Current state

Live at https://chronosauros.github.io/diffuser-generator/ (GitHub Pages, main branch, repository root). Open source
since 24.10.2026: GPL-3.0, the Python reference in `reference/`, the maths in `docs/how-it-works.md`, Discussions on.
`npm test` passes 170 of 170.

24.10.2026, after the review on r/acoustics: the design band now ends at the lower of c / 2w and the flat-plate
frequency N · f0, says which well sets each end and is called design intent, calculated, not measured. 25.10.2026, after
a 2D BEM check: a deep wall narrower than 3.5 design wavelengths (fewer than 7 tiles at the default depth) starts half
an octave above f0. The checks are public in `checks/` (published recipes through the core, the Fourier model against
theory, the 2D BEM simulation with code, raw far fields and figures), written up in `docs/bem-check.md`, and the page
links the BEM check next to the deep design band.

25.10.2026, owner's decision after the night research programme: the deep walls use the best nested orders of the BEM
ranking (walls grow by prefix, tiles already on the wall never move). Series L (P1S / X1, MK4 / Core One)
`5 6 2 7 1 4 3`, seen from the room 3 `2 6 5`, 5 `1 7 2 6 5`, 7 `3 4 1 7 2 6 5`; series S (A1 mini) `1 2 7 3 6 5 4`,
room 4 `3 7 2 1`, 7 `4 5 6 3 7 2 1`. Until then L `6 3 4 2 5 7 1`, S `6 3 5 4 2 7 1`. Changed together:
`reference/out/layouts.json`, the test fixture, LAYOUTS in `index.html`, `docs/how-it-works.md`, `docs/bem-check.md`
and `checks/bem/` (new `order_walls.py` and `order_nested.py`, the five wall runs of the sweep redone for the new
order, the old order frozen as a constant in the order scripts). The Printables set (model 1887298) and its wall
plans still show the old order; that is the owner's call.

25.10.2026, owner's decision: the √2 · f0 band start of deep walls under 7 tiles stays as a cautious start for the new
order, although BEM has those walls diffusing from about f0 (P1S 3 and 5 tiles from the 800 Hz band, A1 mini 4 tiles
from 1000 Hz). Only the docs, the contract and the code comments say so now; the rule and the page did not change.

Next result: owner's decision on the Printables set with the old order; block tiles need 3D BEM; maybe an
advanced mode with N, design frequency and well width as inputs; a measurement is still the only proof. Move any fixes
into `reference/` and the core of `index.html` together.
