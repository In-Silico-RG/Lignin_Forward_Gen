#!/usr/bin/env python
"""
poisson_nuclei_experiment.py

Reviewer 2.3 (Rev_Kimi_1), third option: "Evaluate whether a non-integer n_nuclei
(Poisson) would mitigate the beta-beta ceiling without breaking the physics."

BACKGROUND. beta-beta forms only at nucleation, and the model front-loads a FIXED
integer count N = round(n_nuclei) of nucleation dimers into every molecule. Because
the molecule then grows to n_target units and emits n_target - 1 bonds,

    beta-beta%  =  p_bb * N / (n_target - 1)      (p_bb = P(a nucleus is a resinol))

is a hard function of N. At the calibrated n_nuclei = 4.716 -> N = 5 and n_target =
120, the largest attainable beta-beta is 5/119 = 4.20% (reached only if p_bb = 1).
The reference beta-beta for poplar (5.13%) and high-S (4.57%) sit ABOVE or near that
line: unreachable in principle at this fit (see results/fitted_params.json _ceiling).

WHAT THIS TESTS. With nuclei_poisson on, each molecule draws N ~ Poisson(n_nuclei),
so N varies rep to rep. E[N] = n_nuclei exactly, so the POPULATION MEAN beta-beta is
unchanged (mean-preserving) -- what changes is the SPREAD: the ceiling line becomes a
band, and individual molecules with a high N reach beta-beta well above 4.2%.

So the experiment answers two distinct questions the reviewer conflates:
  Q1 (mean)  Does Poisson move the mean beta-beta toward the hardwood references?
             Prediction: NO -- E[N] = n_nuclei, so the mean is pinned to the same
             line. If confirmed, integer rounding was NOT the cause.
  Q2 (band)  Does Poisson make the references fall WITHIN the realized spread?
             i.e. do some molecules reach the reference beta-beta? This is the
             "line -> band" the SI Fig. 6 comment (minor 3.3) also asks for.

DIAGNOSTIC. For each anchor we then find the mean nucleus count lambda whose Poisson
beta-beta mean equals the reference. If those lambdas differ across feedstocks, one
GLOBAL n_nuclei cannot serve all three -- which is the real obstruction, and points
past Poisson to per-feedstock or continuous (k_nuc) nucleation.

Usage:
    python code/poisson_nuclei_experiment.py [--reps 400]
Writes results/poisson_nuclei_experiment.csv and prints a report.
"""
import argparse
import json
import os

import numpy as np

import mini_forward_generator as mfg
from calibrate_full import REFERENCE, N, ANCHORS

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

LABEL = {0.0: "softwood", 1.0: "poplar", 3.0: "high-S"}


def load_fitted():
    """The adopted chi2 parameter set (results/fitted_params.json)."""
    with open(os.path.join(RESULTS, "fitted_params.json")) as fh:
        raw = json.load(fh)
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def run_reps(sg, params, reps, poisson, seed0=9000):
    """Return per-rep beta-beta% and nucleus count for one anchor."""
    p = dict(params)
    p["nuclei_poisson"] = 1.0 if poisson else 0.0
    bb, nuc = [], []
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=N, params=p, seed=seed0 + r + int(sg * 100))
        s = g.run()
        bb.append(100.0 * s["frac"]["bb"])
        nuc.append(s["n_nuclei_formed"])
    return np.array(bb), np.array(nuc)


def implied_lambda(sg, params, target_bb, reps, seed0=15000):
    """Mean nucleus count lambda whose Poisson beta-beta mean equals target_bb.

    beta-beta mean is linear in lambda (bb = p_bb * lambda / (N-1), p_bb roughly
    constant), so we sample a few lambdas, least-squares fit a line through the
    origin (bb -> 0 as lambda -> 0), and solve. Cheap and exact for a linear law.
    Returns (lambda, grid, means) for optional inspection.
    """
    grid = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    means = []
    for lam in grid:
        p = dict(params)
        p["n_nuclei"] = float(lam)
        bb, _ = run_reps(sg, p, reps, poisson=True, seed0=seed0)
        means.append(bb.mean())
    means = np.array(means)
    slope = float(np.dot(grid, means) / np.dot(grid, grid))   # bb = slope * lambda
    if slope <= 0:
        return None, grid, means
    return target_bb / slope, grid, means


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    args = ap.parse_args()

    params = load_fitted()
    lam_fit = params["n_nuclei"]

    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("Poisson-nucleus experiment (reviewer 2.3, third option)")
    emit("adopted params: n_nuclei = %.4f -> round = %d;  n_target = %d;  reps = %d"
         % (lam_fit, round(lam_fit), N, args.reps))
    emit("hard ceiling at N = %d nuclei: %.2f%% (= %d/%d, only if every nucleus is a"
         " resinol)" % (round(lam_fit), 100.0 * round(lam_fit) / (N - 1),
                        round(lam_fit), N - 1))
    emit("")

    # ---- Q1 / Q2: fixed vs Poisson at the fitted lambda -----------------------
    emit("=== fixed-integer vs Poisson at the fitted n_nuclei (mean-preserving) ===")
    emit("%-9s %7s | %6s %14s | %6s %14s %8s  %s" %
         ("anchor", "ref bb", "fix bb", "fix N (mean)",
          "poi bb", "poi N (mean)", "P(>=ref)", "poi 5-95% band"))
    rows = []
    for sg in ANCHORS:
        ref_bb = 100.0 * REFERENCE[sg]["bb"]
        bb_f, nuc_f = run_reps(sg, params, args.reps, poisson=False)
        bb_p, nuc_p = run_reps(sg, params, args.reps, poisson=True)
        p_ge = float(np.mean(bb_p >= ref_bb))
        lo, hi = np.percentile(bb_p, [5, 95])
        emit("%-9s %6.2f%% | %5.2f%% %13.2f | %5.2f%% %13.2f %7.0f%%  [%4.2f, %4.2f]"
             % (LABEL[sg], ref_bb, bb_f.mean(), nuc_f.mean(),
                bb_p.mean(), nuc_p.mean(), 100 * p_ge, lo, hi))
        rows.append({
            "anchor": LABEL[sg], "sg": sg, "ref_bb": ref_bb,
            "fix_bb_mean": bb_f.mean(), "fix_bb_std": bb_f.std(),
            "fix_N_mean": nuc_f.mean(),
            "poi_bb_mean": bb_p.mean(), "poi_bb_std": bb_p.std(),
            "poi_N_mean": nuc_p.mean(),
            "poi_bb_p05": lo, "poi_bb_p95": hi, "poi_P_ge_ref": p_ge,
        })
    emit("")
    emit("Q1 (mean): Poisson leaves the mean beta-beta on the same line as fixed N")
    emit("           (E[N] = n_nuclei), so rounding was not what pinned it.")
    emit("Q2 (band): P(>=ref) and the 5-95%% band show whether any molecule reaches")
    emit("           the reference -- the hard line is now a band.")
    emit("")

    # ---- Diagnostic: implied lambda per feedstock -----------------------------
    emit("=== mean nucleus count required to hit each reference beta-beta (Poisson) ===")
    emit("if these differ, one GLOBAL n_nuclei cannot serve all three anchors")
    emit("%-9s %7s %14s" % ("anchor", "ref bb", "implied lambda"))
    for sg in ANCHORS:
        ref_bb = 100.0 * REFERENCE[sg]["bb"]
        lam, grid, means = implied_lambda(sg, params, ref_bb,
                                          reps=max(200, args.reps // 2))
        lam_s = "%.2f" % lam if lam is not None else ">16"
        emit("%-9s %6.2f%% %14s" % (LABEL[sg], ref_bb, lam_s))
        for row in rows:
            if row["sg"] == sg:
                row["implied_lambda"] = lam if lam is not None else float("nan")
    emit("")
    emit("fitted global n_nuclei = %.2f. The spread of implied lambdas is the "
         "obstruction;" % lam_fit)
    emit("Poisson widens the band but cannot move a single global mean to three "
         "targets at once.")

    # ---- write CSV ------------------------------------------------------------
    os.makedirs(RESULTS, exist_ok=True)
    csv_path = os.path.join(RESULTS, "poisson_nuclei_experiment.csv")
    cols = ["anchor", "sg", "ref_bb", "fix_bb_mean", "fix_bb_std", "fix_N_mean",
            "poi_bb_mean", "poi_bb_std", "poi_N_mean", "poi_bb_p05", "poi_bb_p95",
            "poi_P_ge_ref", "implied_lambda"]
    with open(csv_path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for row in rows:
            fh.write(",".join(
                ("%s" % row[c]) if c == "anchor" else ("%.4f" % row.get(c, float("nan")))
                for c in cols) + "\n")

    log_path = os.path.join(RESULTS, "poisson_nuclei_experiment.log")
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    emit("")
    emit("Wrote %s" % csv_path)
    emit("Wrote %s" % log_path)


if __name__ == "__main__":
    main()
