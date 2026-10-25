"""Nested wall orders: which 7-tile order serves its smaller walls best, and does the generator use it?

    python checks/bem/order_nested.py [--no-bem] [--workers 6]

A deep wall grows when the owner adds tiles (P1S / X1 class L: 3 -> 5 -> 7 tiles, A1 mini class S: 4 -> 7), and tiles
already stuck on do not move, so the smaller walls are contiguous pieces of the larger one. This script ranks the chains
of such walls from every BEM run in the repository (no new simulation for that part), by the optimiser's score
(order_ranking.py: mean d_n over 0, +-30, +-60 degrees and the bands from f0 to below 7 f0, 3 frequencies per band)
averaged over the wall sizes of the chain:

- prefix chains: 3 tiles = the first three of the order, 5 = the first five, 7 = all (how the generator grows a wall
  now: new tiles are added on one side, the wall order is the one in layouts.json);
- window chains: the smaller wall is any contiguous window of the larger one (the wall grows on either side). Both
  orientations of every order are tried (an order and its mirror image score the same as walls, but their windows differ).

Only chains whose walls were all run in BEM are ranked, so "best" means best among the orders BEM has seen: 145 of the
2544 mirror classes of 7 P1S tiles and 56 of 2544 for the A1 mini (order_walls.py, order_search.py).

Then it checks the orders in reference/out/layouts.json: the walls are prefixes of the order, room_left_to_right is their
reverse (the room list is seen from the room, the model list is STL X, tile order[0] at x = 0), the cached score of every
wall, and unless --no-bem a fresh BEM run of the smallest wall of each class built from room_left_to_right reversed back
into model order, next to the same tiles taken in room order without reversing (the direction mistake). The fresh run
must give the cached number. Summary: out/nested.json, figure out/fig-nested.png.
"""
from __future__ import annotations

import os

for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import base64  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import order_ranking as R  # noqa: E402
import order_walls as W  # noqa: E402

geo = R.geo
# class -> (size, wall sizes of the chain, band-energy files)
CLASSES = {
    "L": ("P1S", (3, 5, 7), ["walls-P1S-bem.jsonl", "search-P1S-bem.jsonl", "ranking-P1S-3-bem.jsonl"]),
    "S": ("A1mini", (4, 7), ["walls-A1mini-bem.jsonl"]),
}
# The orders published until 25.10.2026 and the ones in layouts.json since then (frozen, as in order_walls.py).
OLD, NEW = W.OLD, W.NEW
CLASS_SIZE = {"L": "P1S", "S": "A1mini"}
TOP = 12


def load(files):
    """Band energies by key ("a b c" canonical class, "flat-N"), float64."""
    E = {}
    for name in files:
        for line in (OUT / name).read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            if "key" in r:
                E.setdefault(r["key"], np.frombuffer(base64.b64decode(r["E"]), np.float32).reshape(r["shape"]).astype(float))
            else:                                                       # ranking-P1S-3-bem.jsonl
                o = r["order"]
                k = "flat-3" if o == "flat" else " ".join(map(str, R.canonical(tuple(o))))
                E.setdefault(k, np.array(r["E"], float))
    return E


def score(E, order):
    k, f = " ".join(map(str, R.canonical(tuple(order)))), f"flat-{len(order)}"
    return float(R.score(E[k], E[f]).mean()) if k in E and f in E else None


def best_chain(E, o, sizes, prefix):
    """Best chain of contiguous walls ending in the full order o (None if one of its walls was not run in BEM)."""
    big = len(o)
    best = None

    def windows(w, n):
        return [w[:n]] if prefix else [w[i:i + n] for i in range(len(w) - n + 1)]

    if len(sizes) == 3:
        for w5 in windows(o, sizes[1]):
            for w3 in windows(w5, sizes[0]):
                s = [score(E, w3), score(E, w5), score(E, o)]
                if None not in s and (best is None or np.mean(s) > best[0]):
                    best = (float(np.mean(s)), s, [w3, w5, o])
    else:
        for w4 in windows(o, sizes[0]):
            s = [score(E, w4), score(E, o)]
            if None not in s and (best is None or np.mean(s) > best[0]):
                best = (float(np.mean(s)), s, [w4, o])
    return best


def chains(E, sizes, prefix):
    """One entry per mirror class of 7-tile orders (its best orientation), best first."""
    rows = []
    for k in [k for k in E if not k.startswith("flat") and len(k.split()) == sizes[-1]]:
        c = tuple(map(int, k.split()))
        cands = [b for o in (c, R.mirror(c)) if (b := best_chain(E, o, sizes, prefix))]
        if cands:
            rows.append(max(cands, key=lambda b: b[0]))
    rows.sort(key=lambda r: -r[0])
    return rows


def offsets(ch):
    """Where each wall of the chain starts inside the largest one (0 = at the left end, in model order)."""
    big = tuple(ch[-1])
    return [next(i for i in range(len(big) - len(w) + 1) if big[i:i + len(w)] == tuple(w)) for w in ch]


def entry(row, rank=None):
    mean, s, ch = row
    d = dict(mean=round(mean, 4), scores=[round(x, 4) for x in s], chain=[list(map(int, w)) for w in ch], offsets=offsets(ch))
    if rank is not None:
        d["rank"] = rank
    return d


def chain_of(E, order, sizes, prefix=True):
    """The chain of the walls of one given order (prefixes of it, or the best windows)."""
    return best_chain(E, tuple(order), sizes, prefix)


def rank_of(rows, order):
    k = R.canonical(tuple(order))
    for i, r in enumerate(rows):
        if R.canonical(tuple(r[2][-1])) == k:
            return i + 1
    return None


def fresh_checks(layouts, workers):
    """Fresh BEM of the smallest wall of each class: model order from room_left_to_right reversed, and the naive order."""
    items = []
    for cls, (size, sizes, _) in CLASSES.items():
        n = sizes[0]
        room = tuple(layouts[cls]["room_left_to_right"][str(n)])
        items += [(size, room[::-1]), (size, n), (size, room)]
    with get_context("spawn").Pool(min(workers, len(items))) as pool:
        return dict(pool.map(W.task, items))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-bem", action="store_true", help="skip the fresh BEM runs (cached numbers only)")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    layouts = json.loads((geo.REF / "out" / "layouts.json").read_text(encoding="utf-8"))["deep"]
    res = dict(method=__doc__.strip(), classes={})
    fresh = None if a.no_bem else fresh_checks(layouts, a.workers)
    for cls, (size, sizes, files) in CLASSES.items():
        E = load(files)
        pre, win = chains(E, sizes, True), chains(E, sizes, False)
        old, new = OLD[size], NEW[size]
        c_old, c_new = chain_of(E, old, sizes), chain_of(E, new, sizes)
        lay = layouts[cls]
        walls, room = lay["walls"], lay["room_left_to_right"]
        check = dict(
            order=lay["order"],
            walls_are_prefixes=all(walls[str(n)] == lay["order"][:n] for n in sizes),
            room_is_reversed_walls=all(room[str(n)] == walls[str(n)][::-1] for n in sizes),
            equals_new_order=list(lay["order"]) == list(new),
            cached_scores={str(n): score(E, walls[str(n)]) for n in sizes},
        )
        n0 = sizes[0]
        model = tuple(room[str(n0)])[::-1]
        check["smallest_wall"] = dict(
            tiles=n0, room_left_to_right=room[str(n0)], model_order=list(model), cached=score(E, model),
            cached_room_order_as_model=score(E, tuple(room[str(n0)])))
        if fresh:
            d = lambda arg, n: float(R.score(fresh[(size, arg)].astype(float), fresh[(size, n)].astype(float)).mean())  # noqa: E731
            check["smallest_wall"].update(fresh=round(d(model, n0), 4), fresh_room_order_as_model=round(d(tuple(room[str(n0)]), n0), 4))
        classes7 = [k for k in E if not k.startswith("flat") and len(k.split()) == sizes[-1]]
        res["classes"][cls] = dict(
            size=size, wall_sizes=list(sizes), orders_with_all_walls_in_bem=len(classes7), chains_prefix=len(pre), chains_window=len(win),
            old_prefix=entry(c_old, rank_of(pre, old)), new_prefix=entry(c_new, rank_of(pre, new)),
            top_prefix=[entry(r, i + 1) for i, r in enumerate(pre[:TOP])], top_window=[entry(r, i + 1) for i, r in enumerate(win[:TOP])],
            layouts_check=check)
        print(cls, size, sizes, f"{len(classes7)} orders of {sizes[-1]} tiles, {len(pre)} prefix / {len(win)} window chains")
        print("  old    ", res["classes"][cls]["old_prefix"])
        print("  new    ", res["classes"][cls]["new_prefix"])
        for lab in ("top_prefix", "top_window"):
            for e in res["classes"][cls][lab][:3]:
                print(f"  {lab:10s}", e)
        print("  check  ", check)
    (OUT / "nested.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    figure(res)


def figure(res):
    import style
    from style import plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    for ax, cls, title in ((axs[0], "L", "P1S / X1 (3 -> 5 -> 7 tiles)"), (axs[1], "S", "A1 mini (4 -> 7 tiles)")):
        r = res["classes"][cls]
        other = next(e for e in r["top_window"] if e["chain"] != r["new_prefix"]["chain"])
        bars = [("order until 25.10.2026", r["old_prefix"], style.ORANGE),
                ("current order (prefixes)", r["new_prefix"], style.AQUA),
                ("best other window chain", other, style.BLUE)]
        if cls == "L":
            bars.append(("near tie (prefixes)", r["top_prefix"][1], style.FLAT))
        sizes = r["wall_sizes"]
        x = np.arange(len(sizes) + 1)
        w = 0.8 / len(bars)
        for i, (lab, e, col) in enumerate(bars):
            v = e["scores"] + [e["mean"]]
            xs = x - 0.4 + w * (i + 0.5)
            ax.bar(xs, v, w * 0.92, color=col, label=f"{lab}: {' '.join(map(str, e['chain'][-1]))}")
            for xi, vi in zip(xs, v):
                ax.text(xi, vi + 0.004, f"{vi:.3f}", ha="center", va="bottom", fontsize=8.5, color=style.INK, rotation=90)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n} tiles" for n in sizes] + ["mean"])
        ax.set_ylim(0, 0.72)
        ax.set_ylabel("BEM, mean d_n")
        ax.grid(True, axis="y")
        ax.set_title(title, pad=10, fontsize=13)
        ax.legend(loc="upper left", fontsize=9.5, ncol=1, bbox_to_anchor=(0, 1.0))
    fig.suptitle("Nested wall orders: the walls of a chain are contiguous pieces of the largest one", fontsize=14)
    fig.savefig(OUT / "fig-nested.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
