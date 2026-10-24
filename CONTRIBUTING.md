# Contributing

Thanks for looking. The most useful help right now is a second pair of eyes on the acoustics.

## Reviewing the maths

Read [docs/how-it-works.md](docs/how-it-works.md). It ends with the questions I am least sure about. If you think
something is wrong, open an issue or a discussion and, where you can, point to a source (a book chapter, a paper, a
standard) or show the calculation. Measurements of printed tiles are very welcome, with the setup you used.

## Changing the code

- The Python code in `reference/` is the source of truth; `index.html` is a port of it. A change to the geometry or to
  the formulas goes into both, then `python scripts/reference.py` regenerates the fixtures and `npm test` has to pass.
- Keep `index.html` a single file that opens from disk.
- Say in the pull request what you checked and how. AI-assisted changes are fine; say so.
- Public texts (page, README, docs) must not claim measured performance, because nothing has been measured to
  ISO 17497-2. Numbers from a model or a simulation (Fourier, BEM) are welcome when they say so, with the method, the
  code and the limits, as in [docs/bem-check.md](docs/bem-check.md).

By contributing you agree that your contribution is licensed under the same licence as this repository.
