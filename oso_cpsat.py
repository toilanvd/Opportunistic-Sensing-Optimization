#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exact solver for the Opportunistic Sensing Optimization (OSO) problem
using OR-Tools CP-SAT. Supports both scenarios of the paper.

Reference: V.D. Nguyen, P.L. Nguyen, K. Nguyen, P.T. Do, "Constant approximation
for opportunistic sensing in mobile air quality monitoring system",
Computer Networks 209 (2022) 108646.

Scenarios
---------
general     : each bus route has two distinct paths, at most k turn-on positions
              on each of them (2k in total).
              Input: p q r / n / per route: total_vertices, then for each of the
              two paths: count followed by the coordinates / c / c squares.
              Reference implementation: greedyApproximation_general.cpp

simplified  : the two paths of a route coincide, so a route is a single polyline
              carrying a budget of 2k turn-on positions.
              Input: p q r / n / per route: count followed by the coordinates /
              c / c squares.
              Reference implementation: greedyApproximation_simplified.cpp

In both C++ programs k is the per-path limit and is NOT pre-doubled when read.
The doubling happens inside the algorithm: the simplified DP runs its budget loop
over `for j = 0; j <= 2 * k`, and the general greedy uses turnOnLimit[0] +
turnOnLimit[1] = k + k. This module follows the same convention: pass the same k
as the C++ code and the budget is derived internally.

Model (clean linear reformulation of Section 3.2):

    max  sum_t  y[t]
    s.t. y[t] <= sum over (i,p,j) with A[i,p,j,t]=1 of b[i,p,j]     for all t
         sum_j b[i,p,j] <= budget * a[i]                            for all i,p
         sum_i a[i] <= m
         a, b, y binary

with budget = k in the general scenario and 2k in the simplified one.

Usage
-----
    python3 oso_cpsat.py test.txt --scenario general --max-m 5 --max-k 3
    python3 oso_cpsat.py test.txt --scenario simplified --max-m 60 --max-k 10 --csv out.csv
    python3 oso_cpsat.py test.txt --scenario general --m 3 --k 1
"""

import argparse
import math
import sys
import time

from ortools.sat.python import cp_model

EPS = 1e-9
SCENARIOS = ("general", "simplified")

# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------


class Instance:
    """p x q grid, sensing radius r, n bus routes, c unit critical squares.

    routes[i] is a list of paths; it holds two paths in the general scenario and
    a single one in the simplified scenario. Each path is a list of (x, y)."""

    def __init__(self, p, q, r, routes, squares, scenario):
        self.p = p
        self.q = q
        self.r = r
        self.routes = routes
        self.squares = squares
        self.scenario = scenario
        self.n = len(routes)
        self.c = len(squares)

    def budget(self, k):
        """Turn-on positions allowed on each path of a route."""
        return 2 * k if self.scenario == "simplified" else k


def read_instance(path, scenario):
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    with open(path) as fh:
        tok = fh.read().split()
    pos = 0

    def nxt():
        nonlocal pos
        v = tok[pos]
        pos += 1
        return v

    p, q, r = int(nxt()), int(nxt()), float(nxt())
    n = int(nxt())

    routes = []
    try:
        for _ in range(n):
            if scenario == "general":
                nxt()                   # total vertices of both paths (redundant)
                paths = []
                for _d in range(2):
                    cnt = int(nxt())
                    paths.append([(float(nxt()), float(nxt())) for _ in range(cnt)])
            else:
                cnt = int(nxt())
                paths = [[(float(nxt()), float(nxt())) for _ in range(cnt)]]
            routes.append(paths)

        c = int(nxt())
        squares = [(int(nxt()), int(nxt())) for _ in range(c)]
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"{path}: cannot parse as scenario '{scenario}' ({exc}).\n"
            f"  general    expects, per route: total_vertices, then twice "
            f"(count + coordinates)\n"
            f"  simplified expects, per route: count + coordinates\n"
            f"Try the other --scenario.") from None

    if pos != len(tok):
        raise ValueError(
            f"{path}: parsed as scenario '{scenario}' but consumed only {pos} of "
            f"{len(tok)} tokens. Try the other --scenario.")
    return Instance(p, q, r, routes, squares, scenario)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------


def dist_point_square(px, py, sq):
    """Euclidean distance from a point to the closed unit square
    [sx,sx+1]x[sy,sy+1]; zero inside. Same semantics as disPtsToSqr in C++."""
    sx, sy = sq
    dx = max(sx - px, 0.0, px - (sx + 1.0))
    dy = max(sy - py, 0.0, py - (sy + 1.0))
    return math.hypot(dx, dy)


def observes(px, py, sq, r):
    return dist_point_square(px, py, sq) < r + 1e-9


def _cross_horizontal(p1, p2, y0, xa, xb, out):
    (x1, y1), (x2, y2) = p1, p2
    dx, dy = x2 - x1, y2 - y1
    if abs(dy) <= EPS:
        if abs(y1 - y0) <= EPS and abs(dx) > EPS:
            lo, hi = max(min(x1, x2), xa), min(max(x1, x2), xb)
            if lo <= hi + EPS:
                out.append((lo - x1) / dx)
                out.append((hi - x1) / dx)
        return
    t = (y0 - y1) / dy
    if -EPS <= t <= 1.0 + EPS and xa - EPS <= x1 + t * dx <= xb + EPS:
        out.append(min(max(t, 0.0), 1.0))


def _cross_vertical(p1, p2, x0, ya, yb, out):
    (x1, y1), (x2, y2) = p1, p2
    dx, dy = x2 - x1, y2 - y1
    if abs(dx) <= EPS:
        if abs(x1 - x0) <= EPS and abs(dy) > EPS:
            lo, hi = max(min(y1, y2), ya), min(max(y1, y2), yb)
            if lo <= hi + EPS:
                out.append((lo - y1) / dy)
                out.append((hi - y1) / dy)
        return
    t = (x0 - x1) / dx
    if -EPS <= t <= 1.0 + EPS and ya - EPS <= y1 + t * dy <= yb + EPS:
        out.append(min(max(t, 0.0), 1.0))


def _in_quadrant(px, py, ox, oy, which):
    """which: 0 = NW, 1 = NE, 2 = SE, 3 = SW (as ptsInQuarter in the C++ code)."""
    if which == 0:
        return px <= ox + EPS and py >= oy - EPS
    if which == 1:
        return px >= ox - EPS and py >= oy - EPS
    if which == 2:
        return px >= ox - EPS and py <= oy + EPS
    return px <= ox + EPS and py <= oy + EPS


def _cross_arc(p1, p2, ox, oy, r, which, out):
    (x1, y1), (x2, y2) = p1, p2
    dx, dy = x2 - x1, y2 - y1
    fx, fy = x1 - ox, y1 - oy
    a = dx * dx + dy * dy
    if a <= EPS:
        return
    b = 2.0 * (fx * dx + fy * dy)
    cc = fx * fx + fy * fy - r * r
    disc = b * b - 4.0 * a * cc
    if disc < -EPS:
        return
    sq = math.sqrt(max(disc, 0.0))
    for t in ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)):
        if -EPS <= t <= 1.0 + EPS:
            t = min(max(t, 0.0), 1.0)
            px, py = x1 + t * dx, y1 + t * dy
            if _in_quadrant(px, py, ox, oy, which):
                out.append(t)


def segment_boundary_crossings(p1, p2, sq, r):
    """Parameters t in [0,1] where segment p1p2 crosses the observable boundary
    O(C) of square sq: the Minkowski sum of the square with a disk of radius r
    (4 straight edges + 4 quarter arcs). Mirrors intsSegAndExtSqr in C++."""
    sx, sy = sq
    out = []
    _cross_horizontal(p1, p2, sy - r, sx, sx + 1.0, out)              # bottom
    _cross_horizontal(p1, p2, sy + 1.0 + r, sx, sx + 1.0, out)        # top
    _cross_vertical(p1, p2, sx - r, sy, sy + 1.0, out)                # left
    _cross_vertical(p1, p2, sx + 1.0 + r, sy, sy + 1.0, out)          # right
    _cross_arc(p1, p2, sx, sy + 1.0, r, 0, out)                       # NW corner
    _cross_arc(p1, p2, sx + 1.0, sy + 1.0, r, 1, out)                 # NE corner
    _cross_arc(p1, p2, sx + 1.0, sy, r, 2, out)                       # SE corner
    _cross_arc(p1, p2, sx, sy, r, 3, out)                             # SW corner
    return out


def path_geometry(path):
    seg = [math.dist(path[j], path[j + 1]) for j in range(len(path) - 1)]
    pre = [0.0]
    for s in seg:
        pre.append(pre[-1] + s)
    return seg, pre


def point_at_arclength(path, seg, pre, s):
    if len(path) == 1:
        return path[0]
    j = 0
    while j < len(seg) - 1 and s > pre[j + 1] + EPS:
        j += 1
    if seg[j] <= EPS:
        return path[j]
    t = min(max((s - pre[j]) / seg[j], 0.0), 1.0)
    (x1, y1), (x2, y2) = path[j], path[j + 1]
    return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))


def observable_arclengths(path, sq, r):
    """All arc lengths at which the observable boundary O(sq) is crossed, plus
    the two ends of the path when they observe sq. Empty when sq is unreachable."""
    if len(path) < 2:
        return [0.0] if observes(path[0][0], path[0][1], sq, r) else []
    seg, pre = path_geometry(path)
    L = []
    if observes(path[0][0], path[0][1], sq, r):
        L.append(0.0)
    for j in range(len(seg)):
        for t in segment_boundary_crossings(path[j], path[j + 1], sq, r):
            L.append(pre[j] + t * seg[j])
    if observes(path[-1][0], path[-1][1], sq, r):
        L.append(pre[-1])
    return sorted(L)


# --------------------------------------------------------------------------
# Critical points and coverage
# --------------------------------------------------------------------------


def critical_points_of_path(path, squares, r):
    """Candidate turn-on positions on one path, as arc lengths from its start.

    Theorem 3.1 shows that only the left endpoints of the observable intervals
    need to be considered. Every such endpoint is either the start of the route
    or a crossing of some O(C), so the set below is sufficient for optimality;
    it is the same set the C++ code enumerates."""
    pts = set()
    for sq in squares:
        pts.update(observable_arclengths(path, sq, r))
    ordered, merged = sorted(pts), []
    for s in ordered:
        if not merged or abs(s - merged[-1]) > 1e-7:
            merged.append(s)
    return merged


def check_interval_assumption(inst, samples=200):
    """Section 5.1 assumes that, for every (route, square) pair, the turn-on
    positions observing the square form a single closed segment. The simplified
    C++ code relies on it: it keeps only the hull [L[0], L[-1]] and decides
    coverage by interval membership rather than by distance. Report violations."""
    viol, total = 0, 0
    for i, paths in enumerate(inst.routes):
        for d, path in enumerate(paths):
            if len(path) < 2:
                continue
            seg, pre = path_geometry(path)
            for sq in inst.squares:
                L = observable_arclengths(path, sq, inst.r)
                if not L:
                    continue
                total += 1
                lo, hi = L[0], L[-1]
                for s in range(1, samples):
                    t = lo + (hi - lo) * s / samples
                    x, y = point_at_arclength(path, seg, pre, t)
                    if not observes(x, y, sq, inst.r):
                        viol += 1
                        break
    return viol, total


def build_coverage(inst):
    """cand[i][d] = list of (mask, (x, y), arclength) candidate positions."""
    cand = []
    for paths in inst.routes:
        per_path = []
        for path in paths:
            seg, pre = path_geometry(path) if len(path) > 1 else ([0.0], [0.0])
            entries = []
            for s in critical_points_of_path(path, inst.squares, inst.r):
                px, py = point_at_arclength(path, seg, pre, s)
                mask = 0
                for t, sq in enumerate(inst.squares):
                    if observes(px, py, sq, inst.r):
                        mask |= 1 << t
                if mask:
                    entries.append((mask, (px, py), s))
            per_path.append(entries)
        cand.append(per_path)
    return cand


def reduce_candidates(entries):
    """Drop duplicate and dominated positions on one path, keeping arc-length order.

    If the coverage set of position A is contained in that of B on the same path,
    any solution using A can use B instead, or simply drop A since the budget is
    an upper bound. The reduction therefore preserves the optimum. It also keeps
    the interval structure the simplified DP relies on, because the positions
    covering a square remain a contiguous subsequence."""
    best = {}
    for mask, pt, s in entries:
        if mask not in best:
            best[mask] = (pt, s)
    kept = []
    for m in sorted(best, key=lambda m: -bin(m).count("1")):
        if not any((m & km) == m for km in kept):
            kept.append(m)
    out = [(m, best[m][0], best[m][1]) for m in kept]
    out.sort(key=lambda e: e[2])
    return out


# --------------------------------------------------------------------------
# grOSO (Algorithm 1) with the two submaxSet variants
# --------------------------------------------------------------------------


def _submax_greedy(entries_per_path, limits, covered):
    """Algorithm 3: greedy submaxSet for the general scenario, (1-1/e) of the
    maximum observable set. Ties resolved as in the C++ code (strict >)."""
    picked = [0] * len(entries_per_path)
    acc = 0
    for _ in range(sum(limits)):
        best_gain, best_mask, best_dir = 0, 0, 0
        for d, entries in enumerate(entries_per_path):
            if picked[d] >= limits[d]:
                continue
            for mask, _pt, _s in entries:
                new = mask & ~covered & ~acc
                g = bin(new).count("1")
                if g > best_gain:
                    best_gain, best_mask, best_dir = g, new, d
        picked[best_dir] += 1
        acc |= best_mask
    return acc


def _submax_dp(entries_per_path, limits, covered):
    """Algorithm 2: dynamic programming submaxSet for the simplified scenario.

    Positions are sorted by arc length and a sentinel -infinity is prepended, as
    criticalPts.push_back(-INF) does in the C++ code. Under the single-interval
    assumption a square is covered exactly once, at the first chosen position
    that sees it, so f(i, j) = max(f(i, j-1), max_{u>i} f(u, j-1) + g(i, u))
    computes the true maximum observable set of the route."""
    entries = entries_per_path[0]
    budget = limits[0]
    masks = [0] + [m & ~covered for m, _pt, _s in entries]     # index 0 = -INF
    npts = len(masks)
    pc = bin

    # f[j][i] = best gain using at most j turn-ons among the positions after i
    f = [[0] * npts for _ in range(budget + 1)]
    for j in range(1, budget + 1):
        fj, fp = f[j], f[j - 1]
        for i in range(npts - 1, -1, -1):
            best, mi = fp[i], masks[i]
            for u in range(i + 1, npts):
                v = fp[u] + pc(masks[u] & ~mi).count("1")
                if v > best:
                    best = v
            fj[i] = best

    acc, i, j = 0, 0, budget
    while j > 0:
        if f[j][i] == f[j - 1][i]:
            j -= 1
            continue
        for u in range(i + 1, npts):
            if f[j][i] == f[j - 1][u] + pc(masks[u] & ~masks[i]).count("1"):
                acc |= masks[u]
                i, j = u, j - 1
                break
        else:
            break
    return acc


def greedy_oso(inst, cand, m, k):
    """grOSO (Algorithm 1). Returns (objective, [(route, mask), ...])."""
    submax = _submax_dp if inst.scenario == "simplified" else _submax_greedy
    budget = inst.budget(k)
    limits = [budget] if inst.scenario == "simplified" else [k, k]

    covered, chosen = 0, []
    remaining = set(range(inst.n))
    for _ in range(m):
        best_route, best_mask, best_cnt = None, 0, -1
        for i in sorted(remaining):
            mask = submax(cand[i], limits, covered)
            cnt = bin(mask).count("1")
            if cnt >= best_cnt:                     # ">=" as in the C++ code
                best_route, best_mask, best_cnt = i, mask, cnt
        if best_route is None:
            break
        chosen.append((best_route, best_mask))
        covered |= best_mask
        remaining.discard(best_route)
    return bin(covered).count("1"), chosen


# --------------------------------------------------------------------------
# CP-SAT model
# --------------------------------------------------------------------------


def solve_exact(inst, cand, m, k, time_limit=600.0, workers=8, hint=None,
                log=False):
    model = cp_model.CpModel()
    budget = inst.budget(k)

    a = [model.NewBoolVar(f"a[{i}]") for i in range(inst.n)]
    b, covers = {}, [[] for _ in range(inst.c)]

    for i in range(inst.n):
        for d, entries in enumerate(cand[i]):
            for j, (mask, _pt, _s) in enumerate(entries):
                v = model.NewBoolVar(f"b[{i},{d},{j}]")
                b[(i, d, j)] = v
                mm = mask
                while mm:
                    low = mm & -mm
                    covers[low.bit_length() - 1].append(v)
                    mm ^= low

    y = [model.NewBoolVar(f"y[{t}]") for t in range(inst.c)]

    for t in range(inst.c):
        if covers[t]:
            model.Add(y[t] <= sum(covers[t]))
        else:
            model.Add(y[t] == 0)                    # unreachable square

    for i in range(inst.n):
        for d, entries in enumerate(cand[i]):
            if entries:
                model.Add(sum(b[(i, d, j)] for j in range(len(entries)))
                          <= budget * a[i])

    model.Add(sum(a) <= m)
    model.Maximize(sum(y))

    if hint:
        for i, _mask in hint:
            model.AddHint(a[i], 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit)
    solver.parameters.num_workers = workers
    solver.parameters.log_search_progress = log

    t0 = time.time()
    status = solver.Solve(model)
    elapsed = time.time() - t0
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None

    chosen = []
    for i in range(inst.n):
        if solver.Value(a[i]):
            pos = [(d, pt) for d, entries in enumerate(cand[i])
                   for j, (_m, pt, _s) in enumerate(entries)
                   if solver.Value(b[(i, d, j)])]
            chosen.append((i, pos))

    return {
        "obj": int(round(solver.ObjectiveValue())),
        "bound": int(math.floor(solver.BestObjectiveBound() + 1e-6)),
        "status": solver.StatusName(status),
        "time": elapsed,
        "chosen": chosen,
        "covered": [t for t in range(inst.c) if solver.Value(y[t])],
    }


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


def prepare(inst, dominance=True, verbose=True, check_assumption=True):
    t0 = time.time()
    cand = build_coverage(inst)
    raw = sum(len(e) for r in cand for e in r)
    if dominance:
        cand = [[reduce_candidates(e) for e in r] for r in cand]
    red = sum(len(e) for r in cand for e in r)
    reach = sum(1 for t in range(inst.c)
                if any(mask >> t & 1 for r in cand for e in r for mask, _p, _s in e))
    if verbose:
        print(f"[geometry] candidate positions: {raw} raw -> {red} kept"
              f"   ({time.time() - t0:.2f}s)")
        print(f"[geometry] reachable critical squares: {reach}/{inst.c}"
              f"   (unreachable: {inst.c - reach})")
    if check_assumption:
        viol, total = check_interval_assumption(inst)
        if verbose:
            tag = "holds" if viol == 0 else f"VIOLATED on {viol} pairs"
            print(f"[geometry] Section 5.1 single-interval assumption: {tag}"
                  f" ({total} observable route/square pairs)")
        if viol and inst.scenario == "simplified":
            print("[warn] the simplified C++ code decides coverage by interval "
                  "membership, so its results may overstate coverage here.")
    return cand


def main():
    ap = argparse.ArgumentParser(description="Exact CP-SAT solver for OSO")
    ap.add_argument("instance")
    ap.add_argument("--scenario", choices=SCENARIOS, required=True,
                    help="general: two paths, budget k each; "
                         "simplified: one polyline, budget 2k")
    ap.add_argument("--m", type=int, help="number of sensors (single value)")
    ap.add_argument("--k", type=int, help="turn-on limit per path (single value)")
    ap.add_argument("--max-m", type=int,
                    help="run all m from 1 to max-m (inclusive)")
    ap.add_argument("--max-k", type=int,
                    help="run all k from 1 to max-k (inclusive)")
    ap.add_argument("--time-limit", type=float, default=600.0)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--no-dominance", action="store_true")
    ap.add_argument("--no-hint", action="store_true")
    ap.add_argument("--show-solution", action="store_true")
    ap.add_argument("--csv")
    ap.add_argument("--log", action="store_true")
    args = ap.parse_args()
    t_total = time.time()

    inst = read_instance(args.instance, args.scenario)
    print(f"[input] p={inst.p} q={inst.q} r={inst.r} n={inst.n} c={inst.c}"
          f"   scenario={inst.scenario}")

    cand = prepare(inst, dominance=not args.no_dominance)

    if args.max_m:
        ms = list(range(1, args.max_m + 1))
    elif args.m:
        ms = [args.m]
    else:
        ms = None

    if args.max_k:
        ks = list(range(1, args.max_k + 1))
    elif args.k:
        ks = [args.k]
    else:
        ks = None

    if not ms or not ks:
        ap.error("provide (--m or --max-m) and (--k or --max-k)")

    print()
    print(f"{'m':>3} {'k':>3} {'budget':>6} | {'greedy':>6} {'exact':>6} {'gap%':>6}"
          f" | {'eff%':>6} | {'status':>9} {'time(s)':>8}")
    print("-" * 78)

    rows = []
    for k in ks:
        for m in ms:
            g_obj, g_sol = greedy_oso(inst, cand, m, k)
            res = solve_exact(inst, cand, m, k, time_limit=args.time_limit,
                              workers=args.workers,
                              hint=None if args.no_hint else g_sol, log=args.log)
            if res is None:
                print(f"{m:>3} {k:>3} {inst.budget(k):>6} | {g_obj:>6} {'-':>6}"
                      f" {'-':>6} | {'-':>6} | {'NO SOL':>9}")
                continue
            gap = 100.0 * (res["obj"] - g_obj) / res["obj"] if res["obj"] else 0.0
            eff = 100.0 * g_obj / res["obj"] if res["obj"] else 100.0
            print(f"{m:>3} {k:>3} {inst.budget(k):>6} | {g_obj:>6} {res['obj']:>6}"
                  f" {gap:>6.2f} | {eff:>6.2f} | {res['status']:>9}"
                  f" {res['time']:>8.2f}")
            rows.append((m, k, g_obj, res["obj"], gap, eff, res["status"],
                         res["time"]))
            if args.show_solution:
                for i, pos in res["chosen"]:
                    coords = ", ".join(f"p{d}({x:.2f},{y:.2f})" for d, (x, y) in pos)
                    print(f"        route {i + 1}: {coords}")

    if rows:
        effs = [r[5] for r in rows]
        print("-" * 78)
        print(f"greedy optimal in {sum(1 for r in rows if r[2] == r[3])}/{len(rows)}"
              f" cases | min true efficiency {min(effs):.2f}%"
              f" | mean {sum(effs) / len(effs):.2f}%")

    print(f"[total wall time] {time.time() - t_total:.2f}s")
    if args.csv and rows:
        with open(args.csv, "w") as fh:
            fh.write("scenario,m,k,budget,greedy,exact,gap_percent,"
                     "efficiency_percent,status,time_s\n")
            for r in rows:
                fh.write(f"{inst.scenario},{r[0]},{r[1]},{inst.budget(r[1])},"
                         f"{r[2]},{r[3]},{r[4]:.4f},{r[5]:.4f},{r[6]},{r[7]:.4f}\n")
        print(f"[csv] written to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())