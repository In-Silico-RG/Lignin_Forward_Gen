#!/usr/bin/env python
"""
calibrate.py

Calibrate the fresh forward generator (mini_forward_generator.py) against real
lignin, using the BETA-SIDE-CHAIN normalization: each beta linkage expressed as a
fraction of (beta-O-4 + beta-5 + beta-beta + beta-1). This is the meaningful
comparison for side-chain coupling selectivity; the condensed C-C linkages (5-5,
4-O-5) use a different normalization and are excluded from this objective.

We fit the side-chain rate constants k_b5, k_b1, k_bb and supply_rate (holding
k_bo4 = 8 as the scale) by minimizing squared error to literature beta-side-chain
distributions at three anchors (softwood, poplar, high-S hardwood) with a
gradient-free Nelder-Mead search. Seeds are fixed so the objective is
deterministic. The S/G dependence is NOT fit per point -- one parameter set must
fit all three, so any trend that emerges comes from the S-blocks-C5 chemistry.

Run: conda run -n lignin python code/calibrate.py
Outputs: results/calibration.log (this stdout) + results/mini_forward_calibrated.csv
"""
import os
import sys

import importlib.util
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import rel

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "mfg", os.path.join(HERE, "mini_forward_generator.py"))
mfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mfg)

# ---- reference beta-side-chain fractions (of bo4+b5+bb+b1) -------------------
# Derived from literature linkage inventories (method-dependent; see
# docs/02_biochemistry_evidence.md). Softwood spruce/pine ~ bo4 46 / b5 11 /
# bb 2.5 / b1 1.5; poplar ~ 62 / 7 / 4 / 1.5; high-S birch/euc ~ 78 / 3 / 4 / 1.
def _norm(d):
    t = sum(d.values())
    return {k: v / t for k, v in d.items()}

REFERENCE = {
    0.0: _norm({"bo4": 46.0, "b5": 11.0, "bb": 2.5, "b1": 1.5}),   # softwood
    1.0: _norm({"bo4": 62.0, "b5": 7.0,  "bb": 4.0, "b1": 1.5}),   # poplar
    3.0: _norm({"bo4": 78.0, "b5": 3.0,  "bb": 4.0, "b1": 1.0}),   # high-S HW
}
ANCHORS = [0.0, 1.0, 3.0]
BETA_KEYS = ["bo4", "b5", "bb", "b1"]
FIT_REPS = 10
N = 120
K_BO4 = 8.0

# search bounds (applied by clipping in the objective)
BOUNDS = {"k_b5": (0.05, 40.0), "k_b1": (0.01, 20.0),
          "k_bb": (0.05, 60.0), "pool_target": (1.5, 10.0)}


def beta_sidechain(bonds):
    t = sum(bonds[k] for k in BETA_KEYS)
    return {k: (bonds[k] / t if t else 0.0) for k in BETA_KEYS}


def model_beta(sg, params, reps=FIT_REPS):
    agg = {k: 0 for k in BETA_KEYS}
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=N, params=params,
                               seed=7000 + r + int(sg * 100))
        g.run()
        for k in BETA_KEYS:
            agg[k] += g.bonds[k]
    return beta_sidechain(agg)


def unpack(x):
    k_b5, k_b1, k_bb, pool = np.exp(x)
    k_b5 = float(np.clip(k_b5, *BOUNDS["k_b5"]))
    k_b1 = float(np.clip(k_b1, *BOUNDS["k_b1"]))
    k_bb = float(np.clip(k_bb, *BOUNDS["k_bb"]))
    pool = float(np.clip(pool, *BOUNDS["pool_target"]))
    return {"k_bo4": K_BO4, "k_b5": k_b5, "k_b1": k_b1, "k_bb": k_bb,
            "pool_target": pool}


_eval = {"n": 0, "best": np.inf}


def loss(x):
    params = unpack(x)
    L = 0.0
    for sg in ANCHORS:
        m = model_beta(sg, params)
        for k in BETA_KEYS:
            L += (m[k] - REFERENCE[sg][k]) ** 2
    _eval["n"] += 1
    if L < _eval["best"]:
        _eval["best"] = L
    if _eval["n"] % 15 == 0:
        print("  eval %3d  loss=%.5f  best=%.5f" % (_eval["n"], L, _eval["best"]))
    return L


def main():
    out_dir = os.path.join(HERE, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)

    x0 = np.log([2.24, 0.39, 2.87, 4.0])   # prior fit + pool_target start
    print("Calibrating side-chain constants (beta-side-chain normalization)")
    print("anchors:", ANCHORS, " fit reps:", FIT_REPS, " k_bo4 fixed:", K_BO4)
    print("reference beta-side-chain (bo4/b5/bb/b1):")
    for sg in ANCHORS:
        print("   S/G %.1f: %s" % (sg, {k: round(v, 3) for k, v in REFERENCE[sg].items()}))
    print()

    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": 90, "xatol": 1e-2, "fatol": 1e-5})
    best = unpack(res.x)
    print("\nfitted params:", {k: round(v, 4) for k, v in best.items()})
    print("final loss: %.5f  (evals: %d)\n" % (res.fun, _eval["n"]))

    # comparison table at anchors (higher-rep for reporting)
    print("beta-side-chain fit vs reference (model | ref), %:")
    print("%5s  %-18s %-18s %-18s %-18s" % ("S/G", "bo4", "b5", "bb", "b1"))
    for sg in ANCHORS:
        m = model_beta(sg, best, reps=30)
        cells = []
        for k in BETA_KEYS:
            cells.append("%4.1f|%4.1f" % (100 * m[k], 100 * REFERENCE[sg][k]))
        print("%5.1f  %-18s %-18s %-18s %-18s" % (sg, *cells))

    # full S/G sweep with fitted params: beta-side-chain bo4 + raw full inventory
    print("\nfull S/G sweep with fitted params (raw %% of all bonds; bO4* = beta-side-chain):")
    print("%5s %7s %7s %6s %6s %6s %6s %6s %7s" %
          ("S/G", "bO4*", "bO4", "b5", "bb", "b1", "55", "4o5", "branch"))
    sweep = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
    csv = os.path.join(out_dir, "mini_forward_calibrated.csv")
    with open(csv, "w") as fh:
        fh.write("sg_ratio,bo4_sidechain,bo4,b5,bb,b1,55,4o5,branch\n")
        for sg in sweep:
            raw = {k: [] for k in ["bo4", "b5", "bb", "b1", "55", "4o5"]}
            bsc, branches = [], []
            for r in range(30):
                g = mfg.MiniForwardKMC(sg, n_target=N, params=best, seed=9000 + r + int(sg * 100))
                s = g.run()
                for k in raw:
                    raw[k].append(s["frac"][k])
                bsc.append(beta_sidechain(g.bonds)["bo4"])
                branches.append(s["branch_coeff"])
            mo = {k: 100 * np.mean(raw[k]) for k in raw}
            b = 100 * np.mean(bsc)
            br = np.mean(branches)
            print("%5.2f %6.1f%% %6.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %7.3f" %
                  (sg, b, mo["bo4"], mo["b5"], mo["bb"], mo["b1"], mo["55"], mo["4o5"], br))
            fh.write("%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %
                     (sg, b, mo["bo4"], mo["b5"], mo["bb"], mo["b1"], mo["55"], mo["4o5"], br))
    print("\nWrote %s" % rel(csv))


if __name__ == "__main__":
    main()
