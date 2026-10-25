"""Per-wall BEM ranking of tile orders: P1S / X1 walls of 5 and 7 tiles, A1 mini walls of 4 and 7 tiles.

    python checks/bem/order_walls.py [--workers 8] [--only P1S-5 P1S-7 A1mini-7 A1mini-4] [--cache-only]

Made on the night of 25.10.2026 for the order that was in the generator until then (OLD below); the order in the
generator now (NEW) is marked next to it. Uses the solver, geometry and score of order_ranking.py.

Score per wall, as in layout.deep_layouts and order_ranking.py: mean d_n over the incidences 0, +-30, +-60 degrees and
the bands from f0 to below N f0 (P1S 800 Hz - 4 kHz, A1 mini 1 - 6.3 kHz), 3 frequencies per band, 37 receivers, against
the flat box of the same width. A wall and its mirror image (tile t upside down is 9 - t) score the same, so BEM runs one
order of every mirror class.

Which orders BEM runs:
- A1 mini, 4 tiles: all 432 mirror classes of the 840 orders (exhaustive, like the 3-tile P1S ranking), the old order,
  the Fourier top 40 and a random sample first, so a cut run still answers.
- 5 and 7 tiles: the Fourier top 40 classes, the old and the new order and 16 random classes (the unbiased sample for
  the rank correlation), plus every class already run in the P1S order search (P1S at the same frequencies, cached in
  out/search-P1S-bem.jsonl). The rank of an order is then among the evaluated orders only, 2-10 % of the classes for
  5 and 7 tiles.
Band energies of every run: out/walls-<size>-bem.jsonl (they are in the repository, so the summary can be rebuilt without
BEM: --cache-only stops if a run is missing); summary out/walls.json; figures out/fig-walls-<size>-<n>.png.
"""
from __future__ import annotations

import os

for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import base64  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import bem2d  # noqa: E402
import order_ranking as R  # noqa: E402

geo, L = R.geo, R.L
INC, TH, PER_BAND = R.INC, R.TH, R.PER_BAND
# Wall orders in model (STL X) order: OLD published until 25.10.2026 (P1S / X1 class L 6 3 4 2 5 7 1, A1 mini class S
# 6 3 5 4 2 7 1), NEW in reference/out/layouts.json since then. Frozen on purpose: this ranking was made against the
# old order, and reading layouts.json would change its numbers. Walls are prefixes of the 7-tile order.
OLD = {"P1S": (6, 3, 4, 2, 5, 7, 1), "A1mini": (6, 3, 5, 4, 2, 7, 1)}
NEW = {"P1S": (5, 6, 2, 7, 1, 4, 3), "A1mini": (1, 2, 7, 3, 6, 5, 4)}
WALLS = {"P1S-5": ("P1S", 5), "P1S-7": ("P1S", 7), "A1mini-7": ("A1mini", 7), "A1mini-4": ("A1mini", 4)}
DEADLINE_S = 3 * 3600
SEED = 20261025


def bands(size):
    b = geo.qrd.deep_band(size)
    return [fc for fc in L.BANDS if b["f0"] / 2 ** (1 / 6) <= fc <= b["f_flat"] / 2 ** (1 / 6)]


def freqs(size):
    return L.band_freqs(bands(size), PER_BAND)


def old(size, n):
    return OLD[size][:n]


def new(size, n):
    return NEW[size][:n]


def key(arg):
    return f"flat-{arg}" if isinstance(arg, int) else " ".join(map(str, R.canonical(arg)))


def task(item):
    size, arg = item
    poly = geo.centre_x(geo.flat(size, arg) if isinstance(arg, int) else geo.wall(size, list(arg)))
    rows = []
    for f in freqs(size):
        k = 2 * np.pi * f / bem2d.C
        ms = bem2d.mesh(poly, bem2d.element_size(f))
        rows.append(bem2d.farfield(ms, k, bem2d.solve(ms, k, INC), TH))
    return item, R.band_energy(np.stack(rows)).astype(np.float32)


class Cache:
    def __init__(self, size, cache_only=False):
        self.size, self.E, self.cache_only = size, {}, cache_only
        self.path = OUT / f"walls-{size}-bem.jsonl"
        sources = [self.path] + ([OUT / "search-P1S-bem.jsonl"] if size == "P1S" else [])
        for p in sources:
            if p.exists():
                for line in p.read_text(encoding="utf-8").splitlines():
                    r = json.loads(line)
                    self.E.setdefault(r["key"], np.frombuffer(base64.b64decode(r["E"]), np.float32).reshape(r["shape"]).astype(float))

    def run(self, args, workers, t0, label):
        todo, seen = [], set()
        for a in args:
            k = key(a)
            if k not in self.E and k not in seen:
                seen.add(k)
                todo.append((self.size, a))
        print(f"{label}: {len(todo)} BEM runs ({len(args) - len(todo)} cached or repeated)", flush=True)
        if not todo:
            return True
        if self.cache_only:
            raise SystemExit(f"--cache-only: {len(todo)} runs are not in {self.path.name}, first {todo[0][1]}")
        pool = get_context("spawn").Pool(workers)
        try:
            with self.path.open("a", encoding="utf-8") as fh:
                for i, ((_, a), E) in enumerate(pool.imap(task, todo)):
                    self.E[key(a)] = E.astype(float)
                    fh.write(json.dumps(dict(key=key(a), shape=list(E.shape), E=base64.b64encode(E.tobytes()).decode())) + "\n")
                    fh.flush()
                    if (i + 1) % 10 == 0 or i + 1 == len(todo):
                        print(f"{time.time() - t0:7.0f}s  {label} {i + 1}/{len(todo)}", flush=True)
                    if time.time() - t0 > DEADLINE_S:
                        print(f"deadline reached in {label}: partial result", flush=True)
                        return False
        finally:
            pool.terminate()
            pool.join()
        return True

    def dn(self, order):
        return R.score(self.E[key(tuple(order))], self.E[key(len(order))])


def fourier(size, n):
    m = geo.fourier(size, freqs(size), TH, INC)
    E0 = R.band_energy(m.flat(n).transpose(0, 2, 1))
    pos = [m.shift(j) for j in range(n)]
    out = {}
    for o in itertools.permutations(range(1, 8), n):
        c = R.canonical(o)
        if c not in out:
            p = sum(m.T[t] * pos[j] for j, t in enumerate(c))
            out[c] = R.score(R.band_energy(p.transpose(0, 2, 1)), E0)
    return out                                           # mirror class -> (B, I) d_n


def spearman(a, b):
    ra = np.argsort(np.argsort(-np.asarray(a))); rb = np.argsort(np.argsort(-np.asarray(b)))
    return float(np.corrcoef(ra, rb)[0, 1]) if len(a) > 2 else None


def random_classes(n, k):
    """For 5 and 7 P1S tiles: the 16 random 7-orders of order_search.py (prefixes of uniform random permutations are
    uniform), so their runs are already cached; otherwise a fresh seeded sample of classes."""
    orders = list(itertools.permutations(range(1, 8)))
    rng = np.random.default_rng(SEED)
    rand7 = [orders[i] for i in rng.choice(len(orders), 16, replace=False)]
    return list(dict.fromkeys(R.canonical(o[:n]) for o in rand7))[:k]


def rank_wall(name, workers, t0, results, cache_only=False):
    size, n = WALLS[name]
    cache = Cache(size, cache_only)
    if not cache.run([n], workers, t0, f"{name} flat"):
        return False
    F = fourier(size, n)
    classes = sorted(F, key=lambda c: -F[c].mean())
    cur, nxt = R.canonical(old(size, n)), R.canonical(new(size, n))
    rand = random_classes(n, 16)
    first = list(dict.fromkeys([cur, nxt] + classes[:40] + rand))
    rest = [c for c in classes if c not in first] if name == "A1mini-4" else []
    ok = cache.run(first + rest, workers, t0, name)
    evaluated = [c for c in classes if key(c) in cache.E]
    b = geo.qrd.deep_band(size)
    bnd = bands(size)
    low = [i for i, fc in enumerate(bnd) if fc < b["f_period"]]
    high = [i for i, fc in enumerate(bnd) if fc >= b["f_period"]]
    sb = {c: cache.dn(c) for c in evaluated}
    order_b = sorted(evaluated, key=lambda c: -sb[c].mean())
    frank = {c: i + 1 for i, c in enumerate(classes)}
    best = order_b[0]

    def row(c):
        x, f = sb[c], F[c]
        return dict(order=list(c), mirror=list(R.mirror(c)), bem=round(float(x.mean()), 4),
                    bem_f0_2f0=round(float(x[low].mean()), 4), bem_2f0_7f0=round(float(x[high].mean()), 4),
                    fourier=round(float(f.mean()), 4), fourier_f0_2f0=round(float(f[low].mean()), 4),
                    fourier_2f0_7f0=round(float(f[high].mean()), 4), rank_bem=order_b.index(c) + 1, rank_fourier=frank[c])

    rs = [c for c in rand if c in sb]

    def diff(a, b):
        return dict(all=round(float(sb[a].mean() - sb[b].mean()), 4), f0_2f0=round(float(sb[a][low].mean() - sb[b][low].mean()), 4),
                    f2_7f0=round(float(sb[a][high].mean() - sb[b][high].mean()), 4))

    res = dict(
        size=size, tiles=n, complete=ok, classes_total=len(classes), classes_bem=len(evaluated),
        exhaustive=len(evaluated) == len(classes), bands_hz=bnd, low_bands_hz=[bnd[i] for i in low],
        high_bands_hz=[bnd[i] for i in high], old=row(cur), new=row(nxt), old_order=list(old(size, n)), new_order=list(new(size, n)), best=row(best),
        best_minus_old=diff(best, cur), new_minus_old=diff(nxt, cur), best_minus_new=diff(best, nxt),
        best_f0_2f0=row(max(evaluated, key=lambda c: sb[c][low].mean())),
        best_2f0_7f0=row(max(evaluated, key=lambda c: sb[c][high].mean())),
        random_sample=len(rs),
        spearman_random=spearman([F[c].mean() for c in rs], [sb[c].mean() for c in rs]),
        spearman_evaluated=spearman([F[c].mean() for c in evaluated], [sb[c].mean() for c in evaluated]),
        spearman_evaluated_f0_2f0=spearman([F[c][low].mean() for c in evaluated], [sb[c][low].mean() for c in evaluated]),
        bem_range=[round(float(sb[order_b[-1]].mean()), 4), round(float(np.median([sb[c].mean() for c in evaluated])), 4),
                   round(float(sb[best].mean()), 4)],
        random_range=[round(float(min(sb[c].mean() for c in rs)), 4), round(float(max(sb[c].mean() for c in rs)), 4)] if rs else None,
        top10=[row(c) for c in order_b[:10]],
        evaluated=[row(c) for c in order_b],
    )
    results[name] = res
    (OUT / "walls.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
    figure(name, res, rs)
    print(name, json.dumps({k: res[k] for k in ("classes_bem", "classes_total", "best_minus_old", "new_minus_old",
                                               "spearman_random", "spearman_evaluated")}), "rank old", res["old"]["rank_bem"],
          "rank new", res["new"]["rank_bem"], "best", res["best"]["order"], flush=True)
    return ok


def figure(name, res, rs):
    import style
    from style import plt
    ev = res["evaluated"]
    rset = {tuple(c) for c in rs}
    fig, ax = plt.subplots(figsize=(7.5, 6))
    rest = [e for e in ev if tuple(e["order"]) not in rset]
    ax.scatter([e["fourier"] for e in rest], [e["bem"] for e in rest], s=20, color=style.BLUE, alpha=0.6, lw=0,
               label="Fourier top and cached orders" if not res["exhaustive"] else "every order (mirror pairs merged)")
    rr = [e for e in ev if tuple(e["order"]) in rset]
    if rr:
        ax.scatter([e["fourier"] for e in rr], [e["bem"] for e in rr], s=24, color=style.FLAT, lw=0, label="random orders")
    c, n, b = res["old"], res["new"], res["best"]
    ax.scatter([c["fourier"]], [c["bem"]], s=100, color=style.ORANGE, zorder=3,
               label=f"order until 25.10.2026 {' '.join(map(str, res['old_order']))}")
    ax.scatter([n["fourier"]], [n["bem"]], s=100, color=style.AQUA, edgecolor=style.INK, lw=1.2, zorder=4,
               label=f"current order {' '.join(map(str, res['new_order']))}")
    ax.scatter([b["fourier"]], [b["bem"]], s=190, facecolor="none", edgecolor=style.INK, lw=1.6, zorder=3,
               label=f"best in BEM {' '.join(map(str, b['order']))}")
    ax.grid(True); ax.set_xlabel("Fourier model, mean d_n"); ax.set_ylabel("BEM, mean d_n")
    ax.set_title(f"{res['size']}, {res['tiles']} tiles: {res['classes_bem']} of {res['classes_total']} orders in BEM", pad=10)
    ax.legend(loc="lower right", fontsize=10)
    fig.savefig(OUT / f"fig-walls-{name}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", nargs="*", default=list(WALLS))
    ap.add_argument("--cache-only", action="store_true", help="no BEM: summary and figures from the cached band energies")
    a = ap.parse_args()
    workers = min(a.workers, 8)
    p = OUT / "walls.json"
    results = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    t0 = time.time()
    for name in a.only:
        if not rank_wall(name, workers, t0, results, a.cache_only):
            break
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
