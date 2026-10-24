"""BEM sweep of the deep tiles and walls (rule v3), raw far fields into checks/bem/out/raw/.

    python checks/bem/sweep.py [--workers 8] [--only NAME ...]

Frequencies: the third-octave bands 200 Hz - 10 kHz of layout.py with 6 frequencies per band (as its evaluation),
plus the polar frequencies 200 Hz, 500 Hz, f0, 2 f0, 7 f0 of each size. Incidences 0, +-30, +-45, +-60 degrees (one LU,
seven right-hand sides). Receivers: 37 every 5 degrees (diffusion coefficient), and every 0.5 degree at the polar
frequencies. A configuration that already has its raw file is skipped, so a stopped run resumes.
"""
from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bem2d  # noqa: E402
import geometry as geo  # noqa: E402

RAW = HERE / "out" / "raw"                               # far fields per configuration, 4 MB
SIZES = ("A1mini", "P1S")
INC = [0.0, 30.0, -30.0, 45.0, -45.0, 60.0, -60.0]
TH37 = np.arange(-90, 90.01, 5.0)
THP = np.arange(-90, 90.01, 0.5)
PER_BAND = 6
DEADLINE_S = 3 * 3600


def walls_of(size):
    """Wall orders of a size, copied from the frozen stage 1 layouts.json by qrd.layout_class (like layout.v2_layouts)."""
    s1 = json.loads((geo.REF / "out" / "layouts.json").read_text(encoding="utf-8"))
    cls = geo.qrd.layout_class("deep", geo.qrd.deep_params(size))
    return {int(n): w for n, w in s1["deep"][cls]["walls"].items()}


def polar_freqs(size):
    b = geo.qrd.deep_band(size)
    return [200.0, 500.0, b["f0"], 2 * b["f0"], 7 * b["f0"]]


def configs():
    out = {}
    for size in SIZES:
        walls = walls_of(size)
        out[f"{size}-flat-1"] = dict(size=size, kind="flat", ncol=1)
        out[f"{size}-tile-1"] = dict(size=size, kind="wall", order=[1])
        for n, order in walls.items():
            out[f"{size}-flat-{n}"] = dict(size=size, kind="flat", ncol=n)
            out[f"{size}-wall-{n}"] = dict(size=size, kind="wall", order=order)
            out[f"{size}-periodic-{n}"] = dict(size=size, kind="wall", order=[1] * n)
    return out


def frequencies(size):
    import layout  # noqa: F401 (through geometry)
    fb = geo.layout.band_freqs(geo.layout.BANDS, PER_BAND)
    return np.concatenate([fb, polar_freqs(size)])


def polygon(cfg):
    if cfg["kind"] == "flat":
        return geo.centre_x(geo.flat(cfg["size"], cfg["ncol"]))
    return geo.centre_x(geo.wall(cfg["size"], cfg["order"]))


_CACHE = {}


def task(args):
    name, cfg, f, polar = args
    if name not in _CACHE:
        _CACHE.clear()
        _CACHE[name] = polygon(cfg)
    k = 2 * np.pi * f / bem2d.C
    ms = bem2d.mesh(_CACHE[name], bem2d.element_size(f))
    p = bem2d.solve(ms, k, INC)
    F37 = bem2d.farfield(ms, k, p, TH37).astype(np.complex64)
    FP = bem2d.farfield(ms, k, p, THP).astype(np.complex64) if polar else None
    return name, f, F37, FP, len(ms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    cfgs = configs()
    names = [n for n in cfgs if (not a.only or n in a.only) and not (RAW / f"{n}.npz").exists()]
    # big meshes first, so the pool does not end on one long configuration
    names.sort(key=lambda n: -bem2d.mesh(polygon(cfgs[n]), 0.002).L.size)
    jobs, want = [], {}
    for n in names:
        fs = frequencies(cfgs[n]["size"])
        npol = len(polar_freqs(cfgs[n]["size"]))
        want[n] = len(fs)
        jobs += [(n, cfgs[n], float(f), i >= len(fs) - npol) for i, f in enumerate(fs)]
    print(f"{len(names)} configurations, {len(jobs)} solves, {a.workers} workers", flush=True)
    got = {n: [] for n in names}
    t0 = time.time()
    pool = get_context("spawn").Pool(a.workers)
    try:
        for name, f, F37, FP, nel in pool.imap(task, jobs, chunksize=1):
            got[name].append((f, F37, FP, nel))
            if time.time() - t0 > DEADLINE_S:
                raise TimeoutError("sweep deadline")
            if len(got[name]) == want[name]:
                rows = got.pop(name)
                np.savez_compressed(
                    RAW / f"{name}.npz", cfg=json.dumps(cfgs[name]), f=np.array([r[0] for r in rows]),
                    F37=np.stack([r[1] for r in rows]), elements=np.array([r[3] for r in rows]),
                    fpolar=np.array([r[0] for r in rows if r[2] is not None]),
                    FP=np.stack([r[2] for r in rows if r[2] is not None]), inc=np.array(INC), th37=TH37, thp=THP)
                print(f"{time.time() - t0:7.0f}s  {name}  elements {rows[0][3]}", flush=True)
    finally:
        pool.terminate()
        pool.join()


if __name__ == "__main__":
    main()
