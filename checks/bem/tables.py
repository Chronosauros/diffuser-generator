"""Markdown tables of docs/bem-check.md from out/results.json (printed to stdout).

    python checks/bem/tables.py
"""
import json
from pathlib import Path

R = json.loads((Path(__file__).resolve().parent / "out" / "results.json").read_text(encoding="utf-8"))
REG = ["below c/W", "c/W .. f0", "f0 .. 2 f0", "2 f0 .. 7 f0", "7 f0 band", "above 7 f0"]
REG_EN = {"below c/W": "below c/W", "c/W .. f0": "c/W - f0", "f0 .. 2 f0": "f0 - 2 f0", "2 f0 .. 7 f0": "2 f0 - 7 f0",
          "7 f0 band": "7 f0 band", "above 7 f0": "above 7 f0"}
LABEL = {"A1mini": "A1 mini", "P1S": "P1S / X1"}


def f2(v):
    return "-" if v is None else f"{v:.2f}".replace("-0.00", "0.00")


def hz(v):
    return "never" if v is None else f"{v:g}"


def kinds(n):
    if n == "1":
        return [("tile", "1 tile")]
    return [("layout", f"{n} tiles, wall order"), ("periodic", f"{n} tiles, copies of tile 1")]


for size, d in R["sizes"].items():
    s = d["info"]
    print(f"\n#### {LABEL[size]}: f0 {s['f0']:.0f} Hz, 2 f0 = c / period {s['f2']:.0f} Hz, 7 f0 {s['f7']:.0f} Hz\n")
    print("d_n, mean of incidences 0, ±30, ±60° (the optimiser's score), BEM / Fourier:\n")
    print("| Wall | c/W | " + " | ".join(REG_EN[r] for r in REG) + " |")
    print("|---|---|" + "---|" * len(REG))
    for n, c in sorted(d["cases"].items(), key=lambda x: int(x[0])):
        for k, lab in kinds(n):
            m = c[k]["region_means"]
            cells = []
            for r in REG:
                b, f = m[r]["dn_bem_opt"], m[r]["dn_fourier_opt"]
                cells.append("-" if (b is None and f is None) else f"{f2(b)} / {f2(f)}")
            print(f"| {lab} | {c['c_over_W']:.0f} Hz | " + " | ".join(cells) + " |")
    print("\nLowest third-octave band with d_n >= 0.3 in two bands in a row, BEM / Fourier:\n")
    print("| Wall | c/W | mean of 0, ±30, ±60° | 0° | 45° | -45° |")
    print("|---|---|---|---|---|---|")
    for n, c in sorted(d["cases"].items(), key=lambda x: int(x[0])):
        for k, lab in kinds(n):
            fb = c[k]["first_band_dn_0.3"]
            cells = [f"{hz(fb['dn_bem_' + a])} / {hz(fb['dn_fourier_' + a])}" for a in ("opt", "0", "45", "m45")]
            print(f"| {lab} | {c['c_over_W']:.0f} Hz | " + " | ".join(cells) + " |")
    t = d["cases"]["1"]["tile"]["region_means"]["f0 .. 2 f0"]
    print(f"\nOne tile between f0 and 2 f0 (BEM / Fourier): 0° {f2(t['dn_bem_0'])} / {f2(t['dn_fourier_0'])}, "
          f"45° {f2(t['dn_bem_45'])} / {f2(t['dn_fourier_45'])}, -45° {f2(t['dn_bem_m45'])} / {f2(t['dn_fourier_m45'])}")
