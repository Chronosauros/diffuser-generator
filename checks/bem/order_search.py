"""Search for a better 7-tile order of the P1S / X1 deep walls in BEM, starting from the Fourier model.

    python checks/bem/order_search.py [--workers 12] [--cache-only]

The walls grow by prefix (3 tiles = the first three of the order, 5 = the first five, 7 = all), so an order is judged
by the optimiser's joint score J = (S3 + S5 + S7) / 3, each S the score of order_ranking.py (mean d_n over 0, +-30,
+-60 degrees and the bands 800 Hz - 4 kHz, 3 frequencies per band) for that wall size. S3 in BEM is known for every
prefix (order_ranking.py), so the search only runs 5- and 7-tile walls:

1. The Fourier model scores all 5040 orders.
2. Round 1 runs a shortlist in BEM: the top 40 by a hybrid predictor (exact BEM S3, Fourier S5 and S7 mapped onto
   the BEM scale with the linear fit of the 3-tile ranking), the top 25 by the Fourier J, the order published until
   25.10.2026 (OLD below), and 16 uniformly random orders (the only unbiased estimate of how well the Fourier model ranks 5 and 7 tiles).
3. Local search in BEM only: every swap of two tiles in the 3 best orders, then again from a new best.

A wall and its mirror image (tile t upside down is 9 - t) score the same, so 5- and 7-tile runs are cached by mirror
class. Band energies (B, I, R) of every run, float32 in base64, go to out/search-P1S-bem.jsonl as they come (a stopped
run resumes), the summary to out/search-P1S.json, the figure to out/fig-P1S-search.png. With --cache-only nothing is
simulated: the search is replayed from the cached band energies (it takes the same steps, since they are all in the
cache) and the script stops if a run is missing.

The search was run against the order published until 25.10.2026 and found the order that is in the generator now
(CURRENT); both are frozen below, so the numbers do not change when layouts.json does.
"""
from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

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
sys.path.insert(0, str(HERE))
import bem2d  # noqa: E402
import order_ranking as R  # noqa: E402

geo, L = R.geo, R.L
OUT = R.OUT
JSONL = OUT / "search-P1S-bem.jsonl"
# Deep L order (P1S / X1, MK4 / Core One) in model (STL X) order: published until 25.10.2026, and the one in
# reference/out/layouts.json since then. Frozen on purpose (the search compares against the old one).
OLD = (6, 3, 4, 2, 5, 7, 1)
CURRENT = (5, 6, 2, 7, 1, 4, 3)
SEED = 20261025
DEADLINE_S = 3 * 3600


def task(item):
    kind, arg = item
    poly = geo.centre_x(geo.flat(R.SIZE, arg) if kind == "flat" else geo.wall(R.SIZE, list(arg)))
    rows = []
    for f in R.freqs():
        k = 2 * np.pi * f / bem2d.C
        ms = bem2d.mesh(poly, bem2d.element_size(f))
        rows.append(bem2d.farfield(ms, k, bem2d.solve(ms, k, R.INC), R.TH))
    return item, R.band_energy(np.stack(rows)).astype(np.float32)


def key(item):
    kind, arg = item
    return f"flat-{arg}" if kind == "flat" else " ".join(map(str, R.canonical(arg)))


class Cache:
    def __init__(self):
        self.E = {}
        if JSONL.exists():
            for line in JSONL.read_text(encoding="utf-8").splitlines():
                r = json.loads(line)
                self.E[r["key"]] = np.frombuffer(base64.b64decode(r["E"]), np.float32).reshape(r["shape"]).astype(float)
        for line in R.JSONL.read_text(encoding="utf-8").splitlines():      # the 3-tile runs of order_ranking.py
            r = json.loads(line)
            k = "flat-3" if r["order"] == "flat" else " ".join(map(str, R.canonical(tuple(r["order"]))))
            self.E.setdefault(k, np.array(r["E"]))

    def run(self, items, workers, t0, label, cache_only=False):
        todo, seen = [], set()
        for it in items:
            k = key(it)
            if k not in self.E and k not in seen:
                seen.add(k)
                todo.append(it)
        print(f"{label}: {len(todo)} BEM runs ({len(items) - len(todo)} cached)", flush=True)
        if not todo:
            return
        if cache_only:
            raise SystemExit(f"--cache-only: {len(todo)} runs are not in {JSONL.name}, first {todo[0]}")
        pool = get_context("spawn").Pool(workers)
        try:
            with JSONL.open("a", encoding="utf-8") as fh:
                for i, (it, E) in enumerate(pool.imap_unordered(task, todo)):
                    k = key(it)
                    self.E[k] = E.astype(float)
                    fh.write(json.dumps(dict(key=k, shape=list(E.shape), E=base64.b64encode(E.tobytes()).decode())) + "\n")
                    fh.flush()
                    if (i + 1) % 10 == 0 or i + 1 == len(todo):
                        print(f"{time.time() - t0:7.0f}s  {label} {i + 1}/{len(todo)}", flush=True)
                    if time.time() - t0 > DEADLINE_S:
                        raise TimeoutError("search deadline")
        finally:
            pool.terminate()
            pool.join()

    def dn(self, prefix):
        """(B, I) d_n of the wall `prefix` against the flat box of the same width."""
        return R.score(self.E[" ".join(map(str, R.canonical(prefix)))], self.E[f"flat-{len(prefix)}"])


def needs(order):
    return [("wall", order[:5]), ("wall", order)]


def fourier_scores():
    m = geo.fourier(R.SIZE, R.freqs(), R.TH, R.INC)
    E0 = {n: R.band_energy(m.flat(n).transpose(0, 2, 1)) for n in (3, 5, 7)}
    pos = [m.shift(j) for j in range(7)]
    memo = {}

    def S(prefix):
        if prefix not in memo:
            p = sum(m.T[t] * pos[j] for j, t in enumerate(prefix))
            memo[prefix] = R.score(R.band_energy(p.transpose(0, 2, 1)), E0[len(prefix)])
        return memo[prefix]

    return {o: (S(o[:3]), S(o[:5]), S(o)) for o in itertools.permutations(range(1, 8))}


def spearman(a, b):
    ra = np.argsort(np.argsort(-np.asarray(a))); rb = np.argsort(np.argsort(-np.asarray(b)))
    return float(np.corrcoef(ra, rb)[0, 1])


def swaps(order):
    out = []
    for i, j in itertools.combinations(range(7), 2):
        o = list(order); o[i], o[j] = o[j], o[i]
        out.append(tuple(o))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--cache-only", action="store_true", help="no BEM: replay the search from the cached band energies")
    a = ap.parse_args()
    t0 = time.time()
    cache = Cache()
    four = fourier_scores()
    SF = {o: tuple(float(s.mean()) for s in v) for o, v in four.items()}
    JF = {o: sum(v) / 3 for o, v in SF.items()}
    # linear map Fourier -> BEM from the 3-tile ranking (one entry per mirror pair)
    p3 = sorted({R.canonical(o[:3]) for o in four})
    s3b = np.array([float(cache.dn(p).mean()) for p in p3]); s3f = np.array([SF[next(o for o in four if o[:3] == p)][0] for p in p3])
    slope, icpt = np.polyfit(s3f, s3b, 1)
    S3B = {p: float(cache.dn(p).mean()) for p in {o[:3] for o in four}}
    JH = {o: (S3B[o[:3]] + slope * (SF[o][1] + SF[o][2]) + 2 * icpt) / 3 for o in four}
    orders = list(four)
    rng = np.random.default_rng(SEED)
    rand = [orders[i] for i in rng.choice(len(orders), 16, replace=False)]
    top_h = sorted(orders, key=lambda o: -JH[o])[:40]
    top_f = sorted(orders, key=lambda o: -JF[o])[:25]
    source = {}
    for lab, group in (("random", rand), ("hybrid", top_h), ("fourier", top_f), ("old", [OLD])):
        for o in group:
            source.setdefault(o, lab)
    cache.run([("flat", 5), ("flat", 7)] + [n for o in source for n in needs(o)], a.workers, t0, "round 1", a.cache_only)

    def JB(o):
        return float((cache.dn(o[:3]).mean() + cache.dn(o[:5]).mean() + cache.dn(o).mean()) / 3)

    rounds = []
    best = sorted(source, key=lambda o: -JB(o))[:3]
    starts = list(best)
    for rnd in (2, 3):
        cand = [s for st in starts for s in swaps(st) if s not in source]
        for o in cand:
            source.setdefault(o, f"swap round {rnd}")
        cache.run([n for o in cand for n in needs(o)], a.workers, t0, f"round {rnd}", a.cache_only)
        new_best = max(source, key=JB)
        rounds.append(dict(round=rnd, starts=[list(s) for s in starts], evaluated=len(cand), best=list(new_best),
                           best_J=round(JB(new_best), 4)))
        if new_best in starts or JB(new_best) <= max(JB(s) for s in starts) + 1e-9:
            break
        starts = [new_best]
    seconds = time.time() - t0
    if a.cache_only and (OUT / "search-P1S.json").exists():       # keep the time of the run that made the cache
        seconds = json.loads((OUT / "search-P1S.json").read_text(encoding="utf-8")).get("seconds", seconds)
    summarise(cache, four, SF, JF, JH, source, rand, JB, rounds, (slope, icpt), seconds)


def summarise(cache, four, SF, JF, JH, source, rand, JB, rounds, fit, seconds):
    ev = sorted(source, key=lambda o: -JB(o))
    b = geo.qrd.deep_band(R.SIZE)
    bands = R.bands()
    low = [i for i, fc in enumerate(bands) if fc < b["f_period"]]
    high = [i for i, fc in enumerate(bands) if fc >= b["f_period"]]
    inc = {"0": [0], "30": [1, 2], "60": [3, 4]}
    jf_rank = {o: r + 1 for r, o in enumerate(sorted(four, key=lambda o: -JF[o]))}

    def detail(o):
        d = {}
        for n in (3, 5, 7):
            x = cache.dn(o[:n])
            d[f"S{n}"] = round(float(x.mean()), 4)
            d[f"S{n}_f0_2f0"] = round(float(x[low].mean()), 4)
            d[f"S{n}_2f0_7f0"] = round(float(x[high].mean()), 4)
            for k, ii in inc.items():
                d[f"S{n}_inc{k}"] = round(float(x[:, ii].mean()), 4)
            d[f"S{n}_fourier"] = round(SF[o][(3, 5, 7).index(n)], 4)
        d.update(J=round(JB(o), 4), J_fourier=round(JF[o], 4), J_fourier_rank=jf_rank[o], source=source[o])
        return d

    bem_rank = {o: r + 1 for r, o in enumerate(ev)}
    top5 = ev[:5]
    res = dict(
        method=__doc__.strip(), size=R.SIZE, order_until_2026_10_25=list(OLD), current_order=list(CURRENT), seconds=round(seconds),
        bem_runs=len(cache.E), orders_evaluated=len(ev), linear_fit_fourier_to_bem=[round(float(v), 4) for v in fit],
        rounds=rounds,
        spearman_random=dict(S5=spearman([SF[o][1] for o in rand], [cache.dn(o[:5]).mean() for o in rand]),
                             S7=spearman([SF[o][2] for o in rand], [cache.dn(o).mean() for o in rand]),
                             J=spearman([JF[o] for o in rand], [JB(o) for o in rand])),
        spearman_evaluated=dict(J=spearman([JF[o] for o in ev], [JB(o) for o in ev])),
        old_bem_rank=bem_rank[OLD], current_bem_rank=bem_rank[CURRENT], best=list(ev[0]),
        top5_sources=[source[o] for o in top5],
        compare={"old": dict(order=list(OLD), **detail(OLD)), "best": dict(order=list(ev[0]), **detail(ev[0]))},
        evaluated=[dict(order=list(o), rank_bem=bem_rank[o], **detail(o)) for o in ev],
    )
    (OUT / "search-P1S.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    figure(res)
    print(json.dumps({k: res[k] for k in ("orders_evaluated", "bem_runs", "old_bem_rank", "current_bem_rank", "best", "top5_sources",
                                         "spearman_random", "rounds")}))
    for lab in ("old", "best"):
        c = res["compare"][lab]
        print(lab, c["order"], {k: c[k] for k in ("J", "S3", "S5", "S7", "J_fourier", "J_fourier_rank")})


def figure(res):
    import style
    from style import plt
    ev = res["evaluated"]
    colours = {"random": style.FLAT, "hybrid": style.BLUE, "fourier": style.AQUA, "old": style.ORANGE}
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw=dict(width_ratios=[1.15, 1]))
    ax = axs[0]
    for lab, col, name in (("random", style.FLAT, "random orders"), ("hybrid", style.BLUE, "shortlist: BEM S3 + Fourier S5, S7"),
                           ("fourier", style.AQUA, "shortlist: Fourier J"), ("swap", style.MAGENTA, "local search: swaps in BEM")):
        pts = [e for e in ev if (e["source"].startswith("swap") if lab == "swap" else e["source"] == lab)]
        if pts:
            ax.scatter([e["J_fourier"] for e in pts], [e["J"] for e in pts], s=22, color=col, alpha=0.7, lw=0, label=name)
    for lab, mk in (("old", dict(s=110, color=style.ORANGE)), ("best", dict(s=110, facecolor="none", edgecolor=style.INK, lw=1.8))):
        c = res["compare"][lab]
        name = "order until 25.10.2026" if lab == "old" else "best found = current order"
        ax.scatter([c["J_fourier"]], [c["J"]], zorder=4, label=f"{name} {' '.join(map(str, c['order']))}", **mk)
    ax.grid(True)
    ax.set_xlabel("Fourier model, joint score J"); ax.set_ylabel("BEM, joint score J")
    ax.set_title("Evaluated 7-tile orders (J = mean of the 3-, 5- and 7-tile walls)", pad=10, fontsize=13)
    ax.legend(loc="lower right", fontsize=10)
    ax = axs[1]
    pub, best = res["compare"]["old"], res["compare"]["best"]
    labels = ["3 tiles", "5 tiles", "7 tiles", "joint J"]
    keys = ["S3", "S5", "S7", "J"]
    x = np.arange(len(labels))
    ax.bar(x - 0.2, [pub[k] for k in keys], 0.38, color=style.ORANGE, label=f"until 25.10.2026 {' '.join(map(str, pub['order']))}")
    ax.bar(x + 0.2, [best[k] for k in keys], 0.38, color=style.BLUE, label=f"best found (current) {' '.join(map(str, best['order']))}")
    for i, k in enumerate(keys):
        ax.text(x[i] - 0.2, pub[k] + 0.005, f"{pub[k]:.2f}", ha="center", va="bottom", fontsize=10, color=style.INK)
        ax.text(x[i] + 0.2, best[k] + 0.005, f"{best[k]:.2f}", ha="center", va="bottom", fontsize=10, color=style.INK)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("BEM, mean d_n"); ax.set_ylim(0, max(best["S3"], pub["S3"], best["J"]) * 1.25)
    ax.grid(True, axis="y")
    ax.set_title("BEM score per wall size", pad=10, fontsize=13)
    ax.legend(loc="upper right", fontsize=10)
    fig.suptitle("P1S / X1 deep walls: BEM search for a better wall order", fontsize=15)
    fig.savefig(OUT / "fig-P1S-search.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
