"""Known-answer check of the Fourier model in reference/layout.py (the model that picked the wall orders).

    python checks/fourier_lobes.py [path/to/reference]

Theory (Schroeder; Cox & D'Antonio): a periodic N = 7 QRD at f = m * f0 (m not a multiple of 7) sends equal energy into
every grating lobe sin(theta) = j * lambda / period, apart from the well-width envelope sinc^2; at f = 7 * f0 every well is
a whole number of half wavelengths and only the specular lobe is left (flat plate). The script builds the wells with the
model's own deep_elements() and evaluates the model's own far-field formula exactly at the lobe directions.
"""
import sys
from pathlib import Path

import numpy as np

REF = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "reference"
sys.path.insert(0, str(REF))
import qrd  # noqa: E402
import layout  # noqa: E402

C = qrd.C
# COX standard N = 7 (Cox, Angus, D'Antonio 2003): unit depth 50 mm, well 75 mm, ideal zero-thickness fins
qrd.DEEP["COX"] = dict(N=7, pitch=75.0, length=500.0, dmax=200.0, fin=0.0, edge=0.0, floor=1.6, back=2.0, cap=1.2,
                       chamfer=0.0, bed="recipe")
TILE = 4  # deep_shift(4) = 0: the unshifted sequence


def period_field(size, f, s):
    """One period of the model's p(theta) at s = sin(theta) + sin(psi), the same sum as layout.Deep._tile."""
    k = 2 * np.pi * f / C
    p = 0
    for x, w, d in layout.deep_elements(size, TILE):
        p = p + w * np.sinc(k * w * s / (2 * np.pi)) * np.exp(1j * k * x * s) * np.exp(-2j * k * d)
    return p


def lobes(size, f, well):
    period = qrd.deep_dims(size)[0] / 1000
    lam = C / f
    js = [j for j in range(-20, 21) if abs(j * lam / period) <= 1]
    out = []
    for j in js:
        s = j * lam / period
        E = abs(period_field(size, f, s)) ** 2
        env = (well * np.sinc(2 * np.pi * f / C * well * s / (2 * np.pi))) ** 2   # well-width envelope
        out.append((j, E, E / env if env > 0 else np.nan))
    return out


def report(size, well, label):
    b = qrd.deep_band(size)
    f0 = b["f0"]
    print(f"\n{label}: f0 {f0:.1f} Hz, period {b['period']:.0f} mm, c/period {C / (b['period'] / 1000):.0f} Hz")
    for m in [1, 2, 3, 4, 5, 6, 7]:
        L = lobes(size, m * f0, well)
        E = np.array([e for _, e, _ in L])
        R = np.array([r for _, _, r in L])
        if m == 7:
            side = E[[i for i, (j, _, _) in enumerate(L) if j != 0]]
            spec = E[[i for i, (j, _, _) in enumerate(L) if j == 0]][0]
            print(f"  {m} x f0 = {m * f0:7.0f} Hz: lobes {len(L):2d}, side-lobe energy / specular = "
                  f"{(side.max() / spec if len(side) else 0):.2e}  (flat plate expected: 0)")
        elif len(L) == 1:
            print(f"  {m} x f0 = {m * f0:7.0f} Hz: only the specular lobe propagates (period shorter than the wavelength)")
        else:
            print(f"  {m} x f0 = {m * f0:7.0f} Hz: lobes {len(L):2d}, energy max/min raw {E.max() / E.min():.3f}, "
                  f"after the well-width envelope {np.nanmax(R) / np.nanmin(R):.4f}  (QRD expected: 1)")


report("COX", 0.075, "COX standard N=7, zero-thickness fins")
for size in ("P1S", "A1mini"):
    p = qrd.deep_params(size)
    report(size, (p["pitch"] - p["fin"]) / 1000, f"Our {size} deep tile (fins {p['fin']} mm, edge {p['edge']} mm)")
