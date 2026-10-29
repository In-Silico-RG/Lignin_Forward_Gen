#!/usr/bin/env python
"""
connectivity_experiment.py

The connectivity gap: the generator emits ONE molecule of n_target units
(DP = n_target), while native MWL is ~20 units across many fragments. The paper
reports "% of total modelled bonds", which raises two questions: (a) how much the
reported inventory depends on that single-molecule choice, and (b) whether the model
over-spends condensed bonds (5-5, 4-O-5) to force single-molecule connectivity,
which would make the 4-O-5 residual a topology artifact.

Two parts, each answering one half honestly.

PART A -- DP invariance (defends the reported quantity). Sweep n_target (the DP) with
the nucleus DENSITY held fixed: n_nuclei scaled proportionally to n_target and
floored at 1 (a molecule needs at least one nucleus; at DP 20 the scaled count is
0.79, so the floor is what keeps the small-DP end in a physical regime rather than
the degenerate no-nucleus start). If the six-linkage FRACTIONS are ~flat across a
12x range of DP, then the quantity the paper reports does not depend on the
single-molecule size -- DP is a simulation size, as the docstring already states.
beta-beta is the one exception by construction (it is set by nucleus count / bonds),
and that is the ceiling story, not a new problem.

PART B -- fragmentation vs condensed content (tests that causal claim directly). At
n_target = 120, dial down oo_scale (the multiplier on the 5-5 / 4-O-5 channel, the
ONLY channel that welds fragments -- beta-coupling never joins two chains). Report
n_frags, the largest fragment, the condensed inventory and chi2. This exposes the
real structural fact: in this model condensation IS the welding, so the only way to
get many fragments (realistic topology) is to suppress 5-5 / 4-O-5 below their
reference. You cannot have native fragmentation AND the reference condensed content
at once -- which is exactly why the paper should claim a linkage inventory under a
forced-linear topology, not native topology. Whether the 4-O-5 over-prediction eases
as fragments appear is read directly from the 4-O-5 column.

Usage:
    python code/connectivity_experiment.py [--reps 150]
Writes results/connectivity_experiment.csv and .log and prints a report.
"""
import argparse
import json
import os
import sys


import numpy as np

import mini_forward_generator as mfg
from calibrate_full import REFERENCE, N, ANCHORS, ALL_KEYS, sigma

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import rel

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

LABEL = {0.0: "softwood", 1.0: "poplar", 3.0: "high-S"}
DP_GRID = [20, 40, 60, 120, 240]
OO_GRID = [0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]


def load_fitted(path=None):
    """Constants from a fitted_params*.json (default: the v1.0.1 three-anchor fit)."""
    with open(path or os.path.join(RESULTS, "fitted_params.json")) as fh:
        raw = json.load(fh)
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def run_cell(sg, base, reps, seed0, n_target=None, n_nuclei=None,
             oo_scale=None, poisson=False):
    p = dict(base)
    if n_nuclei is not None:
        p["n_nuclei"] = float(n_nuclei)
    if oo_scale is not None:
        p["oo_scale"] = float(oo_scale)
    if poisson:
        p["nuclei_poisson"] = 1.0
    nt = n_target if n_target is not None else N
    agg = {k: 0 for k in ALL_KEYS}
    n_frags, max_frag = [], []
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=nt, params=p, seed=seed0 + r + int(sg * 100))
        s = g.run()
        for k in ALL_KEYS:
            agg[k] += g.bonds[k]
        n_frags.append(s["n_frags"])
        max_frag.append(s["max_frag"])
    t = sum(agg.values())
    fr = {k: (agg[k] / t if t else 0.0) for k in ALL_KEYS}
    chi2 = sum(((fr[k] - REFERENCE[sg][k]) / sigma(REFERENCE[sg][k])) ** 2
               for k in ALL_KEYS)
    return {"fr": fr, "n_frags": float(np.mean(n_frags)),
            "max_frag": float(np.mean(max_frag)), "chi2": chi2}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=150)
    ap.add_argument("--params", default=None,
                    help="fitted_params*.json to run at (default: results/fitted_params.json, v1.0.1)")
    ap.add_argument("--tag", default="", help="suffix for the output files")
    args = ap.parse_args()

    base = load_fitted(args.params)
    lam = base["n_nuclei"]
    rows = []
    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("Connectivity experiment (DP-invariance and the welding channel)")
    emit("adopted params; base n_nuclei = %.4f at n_target = %d; reps = %d"
         % (lam, N, args.reps))
    emit("")

    # ---- PART A: DP invariance ------------------------------------------------
    emit("=== PART A: linkage inventory vs DP (nucleus density fixed) ===")
    emit("n_nuclei scaled proportionally to n_target and floored at 1; the generator rounds")
    emit("it to an integer, so the nucleus DENSITY is only approximately constant across DP")
    emit("(the integer actually used is n_nuc).")
    for sg in ANCHORS:
        emit("-- %s (S/G %.1f) --  ref: %s" %
             (LABEL[sg], sg, " ".join("%s %.1f" % (k, 100 * REFERENCE[sg][k])
                                      for k in ALL_KEYS)))
        emit("%6s %8s %7s %8s | %6s %6s %6s %6s %6s %6s | %6s"
             % ("DP", "n_nuc", "nfrags", "maxfrag", "bO4", "b5", "bb", "b1", "5-5",
                "4O5", "chi2"))
        for dp in DP_GRID:
            nn = max(1.0, lam * dp / N)
            d = run_cell(sg, base, args.reps, seed0=31000, n_target=dp,
                         n_nuclei=nn, poisson=False)
            fr = d["fr"]
            emit("%6d %8.2f %7.2f %8.1f | %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f | %6.2f"
                 % (dp, nn, d["n_frags"], d["max_frag"], 100 * fr["bo4"],
                    100 * fr["b5"], 100 * fr["bb"], 100 * fr["b1"], 100 * fr["55"],
                    100 * fr["4o5"], d["chi2"]))
            rows.append({"part": "A", "anchor": LABEL[sg], "sg": sg, "dp": dp,
                         "oo_scale": 1.0, "n_nuclei": nn, "n_frags": d["n_frags"],
                         "max_frag": d["max_frag"], "chi2": d["chi2"],
                         **{k: 100 * fr[k] for k in ALL_KEYS}})
        emit("")

    # ---- PART B: fragmentation vs condensed content ---------------------------
    emit("=== PART B: fragmentation vs condensed inventory (n_target = %d) ===" % N)
    emit("oo_scale multiplies the 5-5 / 4-O-5 channel, the ONLY welder. Lowering it")
    emit("makes fragments -- at the cost of the condensed inventory it is calibrated to.")
    for sg in ANCHORS:
        ref55, ref4o5 = 100 * REFERENCE[sg]["55"], 100 * REFERENCE[sg]["4o5"]
        emit("-- %s (S/G %.1f) --  ref 5-5 = %.1f, ref 4-O-5 = %.1f --"
             % (LABEL[sg], sg, ref55, ref4o5))
        emit("%9s %7s %8s | %6s %6s | %6s"
             % ("oo_scale", "nfrags", "maxfrag", "5-5", "4O5", "chi2"))
        for oo in OO_GRID:
            d = run_cell(sg, base, args.reps, seed0=42000, oo_scale=oo)
            fr = d["fr"]
            emit("%9.3f %7.2f %8.1f | %5.1f %5.1f | %6.2f"
                 % (oo, d["n_frags"], d["max_frag"], 100 * fr["55"], 100 * fr["4o5"],
                    d["chi2"]))
            rows.append({"part": "B", "anchor": LABEL[sg], "sg": sg, "dp": N,
                         "oo_scale": oo, "n_nuclei": lam, "n_frags": d["n_frags"],
                         "max_frag": d["max_frag"], "chi2": d["chi2"],
                         **{k: 100 * fr[k] for k in ALL_KEYS}})
        emit("")

    emit("=== reading ===")
    emit("A: if the six fractions are flat across DP 20-240 (bar beta-beta), the")
    emit("   reported linkage inventory does not depend on the single-molecule size.")
    emit("B: n_frags rises only as oo_scale -> 0, where 5-5 and 4-O-5 collapse below")
    emit("   reference. Native fragmentation and the reference condensed content are")
    emit("   mutually exclusive in this model, because condensation IS the welding.")
    emit("   The 4-O-5 column shows whether its over-prediction is a welding artifact.")

    os.makedirs(RESULTS, exist_ok=True)
    csv_path = os.path.join(RESULTS, "connectivity_experiment%s.csv" % args.tag)
    cols = ["part", "anchor", "sg", "dp", "oo_scale", "n_nuclei", "n_frags",
            "max_frag", "chi2"] + ALL_KEYS
    with open(csv_path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for row in rows:
            fh.write(",".join(("%s" % row[c]) if c in ("part", "anchor")
                              else ("%.4f" % row[c]) for c in cols) + "\n")
    log_path = os.path.join(RESULTS, "connectivity_experiment%s.log" % args.tag)
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    emit("")
    emit("Wrote %s" % rel(csv_path))
    emit("Wrote %s" % rel(log_path))


if __name__ == "__main__":
    main()
