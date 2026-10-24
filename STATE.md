# Current state

Live at https://chronosauros.github.io/diffuser-generator/ (GitHub Pages, main branch, repository root). Open source
since 24.10.2026: GPL-3.0, the Python reference in `reference/`, the maths in `docs/how-it-works.md`, Discussions on.
`npm test` passes 170 of 170.

24.10.2026, after the review on r/acoustics: the design band now ends at the lower of c / 2w and the flat-plate
frequency N · f0, says which well sets each end and is called design intent, calculated, not measured. 25.10.2026, after a 2D BEM check: a deep wall narrower than 3.5 design wavelengths (fewer than 7 tiles at the default depth) starts half an octave above f0.

Next result: more answers to the review (a BEM check of the Fourier model, maybe an advanced mode with N, design
frequency and well width as inputs); move any fixes into `reference/` and the core of `index.html` together.
