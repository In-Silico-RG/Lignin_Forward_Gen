#!/usr/bin/env python
"""
calibrate_full.py

Full-inventory calibration of the fresh forward generator. Unlike calibrate.py
(which fit only the beta-side-chain distribution and left the condensed channel
unconstrained -> over-condensed softwood, branch 0.31), this fits the ENTIRE
modeled linkage inventory -- beta-O-4, beta-5, beta-beta, beta-1, 5-5, 4-O-5 --
as fractions of TOTAL linkages (the standard NMR reporting convention), and fits
the condensed rate constants k_55, k_4o5 together with k_b5, k_b1, k_bb and
pool_target (k_bo4 = 8 fixed as scale).

One parameter set must fit all three S/G anchors, so any S/G trend comes from the
S-blocks-C5 chemistry, not per-point fitting. Seeds fixed -> deterministic loss.

Run: conda run -n lignin python code/calibrate_full.py
Outputs: results/calibration_full.log, results/mini_forward_full_calibrated.csv
"""
import os
import importlib.util
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "mfg", os.path.join(HERE, "mini_forward_generator.py"))
mfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mfg)

ALL_KEYS = ["bo4", "b5", "bb", "b1", "55", "4o5"]


def _norm(d):
    t = sum(d.values())
    return {k: v / t for k, v in d.items()}


# Reference linkage inventories as % of total (modeled) linkages, from literature
# (method-dependent; see docs/02). Softwood is notably condensed (high 5-5);
# hardwood/high-S is beta-O-4-rich and largely uncondensed.
REFERENCE = {
    0.0: _norm({"bo4": 46.0, "b5": 11.0, "bb": 2.5, "b1": 2.0, "55": 12.0, "4o5": 4.0}),
    1.0: _norm({"bo4": 62.0, "b5": 6.0,  "bb": 4.0, "b1": 1.5, "55": 3.0,  "4o5": 1.5}),
    3.0: _norm({"bo4": 78.0, "b5": 3.0,  "bb": 4.0, "b1": 1.0, "55": 1.0,  "4o5": 0.5}),
}
ANCHORS = [0.0, 1.0, 3.0]
FIT_REPS = 8
N = 120
K_BO4 = 8.0

# fit order: k_b5, k_b1, k_bb, k_55, k_4o5, n_nuclei (pool_target fixed at
# generator default 1.0 -> endwise/linear; beta-beta now set by n_nuclei)
BOUNDS = [("k_b5", 0.05, 40.0), ("k_b1", 0.01, 20.0), ("k_bb", 0.05, 60.0),
          ("k_55", 0.005, 40.0), ("k_4o5", 0.005, 40.0), ("n_nuclei", 1.0, 30.0)]


def full_fracs(bonds):
    t = sum(bonds[k] for k in ALL_KEYS)
    return {k: (bonds[k] / t if t else 0.0) for k in ALL_KEYS}


def model_full(sg, params, reps=FIT_REPS):
    agg = {k: 0 for k in ALL_KEYS}
    branch = []
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=N, params=params, seed=7000 + r + int(sg * 100))
        s = g.run()
        for k in ALL_KEYS:
            agg[k] += g.bonds[k]
        branch.append(s["branch_coeff"])
    return full_fracs(agg), float(np.mean(branch))


def unpack(x):
    vals = np.exp(x)
    params = {"k_bo4": K_BO4}
    for (name, lo, hi), v in zip(BOUNDS, vals):
        params[name] = float(np.clip(v, lo, hi))
    return params


_eval = {"n": 0, "best": np.inf}


def loss(x):
    params = unpack(x)
    L = 0.0
    for sg in ANCHORS:
        m, _ = model_full(sg, params)
        for k in ALL_KEYS:
            L += (m[k] - REFERENCE[sg][k]) ** 2
    _eval["n"] += 1
    if L < _eval["best"]:
        _eval["best"] = L
    if _eval["n"] % 20 == 0:
        print("  eval %3d  loss=%.5f  best=%.5f" % (_eval["n"], L, _eval["best"]))
    return L


def main():
    out_dir = os.path.join(HERE, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)

    # start with lower condensed constants (now decoupled from beta-beta) + more nuclei
    x0 = np.log([2.25, 0.40, 4.0, 0.8, 0.15, 12.0])
    print("Full-inventory calibration (fraction of total linkages)")
    print("anchors:", ANCHORS, " fit reps:", FIT_REPS, " k_bo4 fixed:", K_BO4)
    print("reference inventory (%% of total linkages):")
    for sg in ANCHORS:
        print("   S/G %.1f: %s" % (sg, {k: round(100 * v, 1) for k, v in REFERENCE[sg].items()}))
    print()

    # wider initial simplex (~x1.5 per dim in log space) to force exploration
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy(); pt[i] += 0.4
        simplex.append(pt)
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": 300, "xatol": 1e-2, "fatol": 1e-6,
                            "initial_simplex": np.array(simplex)})
    best = unpack(res.x)
    print("\nfitted params:", {k: round(v, 4) for k, v in best.items()})
    print("final loss: %.5f  (evals: %d)\n" % (res.fun, _eval["n"]))

    print("full-inventory fit vs reference (model | ref, %% of total linkages):")
    print("%5s  %-11s %-11s %-11s %-11s %-11s %-11s  %s" %
          ("S/G", "bO4", "b5", "bb", "b1", "5-5", "4-O-5", "branch"))
    for sg in ANCHORS:
        m, br = model_full(sg, best, reps=30)
        cells = ["%4.1f|%4.1f" % (100 * m[k], 100 * REFERENCE[sg][k]) for k in ALL_KEYS]
        print("%5.1f  %-11s %-11s %-11s %-11s %-11s %-11s  %.3f" % (sg, *cells, br))

    # full S/G sweep with fitted params
    print("\nfull S/G sweep with fitted params (%% of total linkages):")
    print("%5s %7s %6s %6s %6s %6s %6s %7s" %
          ("S/G", "bO4", "b5", "bb", "b1", "5-5", "4o5", "branch"))
    sweep = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
    csv = os.path.join(out_dir, "mini_forward_full_calibrated.csv")
    with open(csv, "w") as fh:
        fh.write("sg_ratio,bo4,b5,bb,b1,55,4o5,branch\n")
        for sg in sweep:
            m, br = model_full(sg, best, reps=30)
            print("%5.2f %6.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %7.3f" %
                  (sg, 100 * m["bo4"], 100 * m["b5"], 100 * m["bb"], 100 * m["b1"],
                   100 * m["55"], 100 * m["4o5"], br))
            fh.write("%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %
                     (sg, m["bo4"], m["b5"], m["bb"], m["b1"], m["55"], m["4o5"], br))
    print("\nfitted params:", {k: round(v, 4) for k, v in best.items()})
    print("Wrote %s" % csv)


if __name__ == "__main__":
    main()
