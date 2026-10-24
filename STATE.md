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

Next result: open points of the BEM check (is the published wall order the best one or only better than copies;
block tiles need 3D BEM), maybe an advanced mode with N, design frequency and well width as inputs; a measurement is
still the only proof. Move any fixes into `reference/` and the core of `index.html` together.
