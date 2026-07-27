"""
sigma_4o5_sensitivity.py

How much of the fit rests on a 4-O-5 reference the literature does not agree on?

The chi2 objective weights each linkage by an estimated measurement sigma,
sigma_k = max(0.05*ref_k, 0.01). For 4-O-5 the FLOOR binds at every anchor: the
reference is 1.5/1.5/0.5 per 100 C9, so 5% of it is far below one point and the
fit is told 4-O-5 is known to +/- 1 percentage point.

That is the tightest claim in the objective, and it is made about the least
certain number in the reference. Capanema et al. 2005 (JAFC 53:9639, Table 1)
tabulate five independent structural determinations of hardwood lignins:

    4-O-5 per 100 Ar    beech  birch  E.globulus  E.grandis(22)  E.grandis(here)
      guaiacyl            1.5    1.0      1.5          6.0            3.0
      syringyl             -     5.5     10.0           -             6.0
      TOTAL               1.5    6.5     11.5          6.0            9.0

A spread of 1.5 to 11.5 per 100 Ar -- nearly 8x -- against 5-5 over the same five
determinations spanning only 2 to 6. 4-O-5 is not merely small, it is the one
linkage whose true abundance is openly disputed: thioacidolysis does not detect
phenolic 4-O-5 at all, and NMR quantification needed purpose-built model
compounds before it worked.

So this script asks the obvious question: if we stop pretending 4-O-5 is known to
one point and give it a sigma matching the literature spread, does anything we
report change? It refits all seven constants at each sigma_4o5 and reports the
fitted values, the per-wood RMSE, and the per-linkage residuals.

Two readings are possible and the script distinguishes them:

  (a) The fit is INSENSITIVE -- the constants and RMSE barely move. Then the
      adopted floor is harmless and the 4-O-5 residual is a real model error.
  (b) The fit MOVES -- relaxing 4-O-5 buys a better fit elsewhere. Then part of
      the reported error is the objective straining against a disputed number,
      and the headline RMSE is pessimistic for the wrong reason.

Usage:
    python code/sigma_4o5_sensitivity.py [--reps 8] [--report-reps 30]
Writes results/sigma_4o5_sensitivity.csv and prints a report.
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calibrate_full as cf  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

# sigma_4o5 in FRACTION of total linkages. 0.01 (one percentage point) is the
# adopted floor. The rest bracket the Capanema Table 1 spread once normalised:
# 1.5-11.5 per 100 Ar over a ~66-88 bond/100 C9 inventory is roughly 2-15% of
# linkages, so even 0.04 is a conservative reading of the disagreement.
SIGMA_GRID = [0.01, 0.02, 0.03, 0.04, 0.06]

# The adopted chi2 optimum (results/fitted_params.json), used as the start for
# every refit so each one begins from the same place and differences are due to
# the objective, not the search path.
#         k_b5    k_b1    k_bb     k_55    k_4o5   n_nuclei  k_bo4_nuc
BASE_X0 = [2.4904, 0.2715, 16.7825, 0.7108, 0.1737, 4.716,    0.9213]


def sigma_for(key, ref_frac, sigma_4o5):
    """Adopted sigma, except 4-O-5 gets the sweep value."""
    if key == "4o5":
        return sigma_4o5
    return max(cf.REL_SIGMA * ref_frac, cf.FLOOR_SIGMA)


def chi2_for(params, sigma_4o5, reps):
    total = 0.0
    for sg in cf.ANCHORS:
        m, _ = cf.model_full(sg, params, reps=reps)
        for k in cf.ALL_KEYS:
            r = cf.REFERENCE[sg][k]
            total += ((m[k] - r) / sigma_for(k, r, sigma_4o5)) ** 2
    return total


def per_wood_rmse(params, reps):
    """RMSE in percentage points, per anchor, and the per-linkage residuals."""
    out = {}
    for sg in cf.ANCHORS:
        m, _ = cf.model_full(sg, params, reps=reps)
        res = {k: 100 * (m[k] - cf.REFERENCE[sg][k]) for k in cf.ALL_KEYS}
        rmse = math.sqrt(sum(v ** 2 for v in res.values()) / len(cf.ALL_KEYS))
        out[sg] = (rmse, res, {k: 100 * m[k] for k in cf.ALL_KEYS})
    return out


def fit_at(sigma_4o5, reps):
    """Refit all seven constants with 4-O-5 given sigma_4o5."""
    n_eval = {"n": 0}

    def loss(x):
        n_eval["n"] += 1
        return chi2_for(cf.unpack(x), sigma_4o5, reps)

    x0 = np.log(BASE_X0)
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy()
        pt[i] += 0.4
        simplex.append(pt)
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": 400, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    return cf.unpack(res.x), res.fun, n_eval["n"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=cf.FIT_REPS,
                    help="replicates per evaluation during the fit")
    ap.add_argument("--report-reps", type=int, default=30,
                    help="replicates used for the reported numbers")
    args = ap.parse_args()

    print("4-O-5 sigma sensitivity")
    print("  adopted objective: sigma_k = max(%.2f*ref_k, %.2f)"
          % (cf.REL_SIGMA, cf.FLOOR_SIGMA))
    print("  4-O-5 reference: %s per 100 C9 -> floor binds at every anchor"
          % [cf.RAW_REFERENCE[sg]["4o5"] for sg in cf.ANCHORS])
    print("  fit reps = %d, report reps = %d\n" % (args.reps, args.report_reps))

    rows = []
    fits = {}
    for s in SIGMA_GRID:
        params, fun, nev = fit_at(s, args.reps)
        wood = per_wood_rmse(params, args.report_reps)
        mean_rmse = sum(w[0] for w in wood.values()) / len(wood)
        fits[s] = (params, wood, mean_rmse)
        tag = "  (adopted)" if s == cf.FLOOR_SIGMA else ""
        print("sigma_4o5 = %.2f (%.0f pt)%s   evals %d   mean RMSE %.3f pts"
              % (s, 100 * s, tag, nev, mean_rmse))
        print("   params: %s" % {k: round(v, 4) for k, v in params.items()
                                 if k != "k_bo4"})
        for sg in cf.ANCHORS:
            rmse, res, model = wood[sg]
            print("   S/G %.1f  RMSE %5.2f   4-O-5 %4.1f (ref %4.1f, %+5.2f)"
                  "   5-5 %4.1f (ref %4.1f, %+5.2f)"
                  % (sg, rmse, model["4o5"], 100 * cf.REFERENCE[sg]["4o5"],
                     res["4o5"], model["55"], 100 * cf.REFERENCE[sg]["55"],
                     res["55"]))
            rows.append({
                "sigma_4o5": s, "sg": sg, "rmse_pts": round(rmse, 4),
                **{"model_" + k: round(model[k], 3) for k in cf.ALL_KEYS},
                **{"resid_" + k: round(res[k], 3) for k in cf.ALL_KEYS},
                **{"param_" + k: round(v, 5) for k, v in params.items()
                   if k != "k_bo4"},
            })
        print()

    base = fits[cf.FLOOR_SIGMA]
    print("CHANGE FROM THE ADOPTED FIT (sigma_4o5 = 0.01)")
    print("%-10s %10s %12s %12s %12s"
          % ("sigma_4o5", "mean RMSE", "d(mean RMSE)", "max |dk|/k", "4-O-5 resid"))
    for s in SIGMA_GRID:
        params, wood, mean_rmse = fits[s]
        rel = max(abs(params[k] - base[0][k]) / base[0][k]
                  for k in params if k != "k_bo4")
        mean_4o5 = sum(wood[sg][1]["4o5"] for sg in cf.ANCHORS) / len(cf.ANCHORS)
        print("%-10.2f %10.3f %+12.3f %11.1f%% %+12.2f"
              % (s, mean_rmse, mean_rmse - base[2], 100 * rel, mean_4o5))

    out = os.path.join(RESULTS, "sigma_4o5_sensitivity.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote %s" % os.path.normpath(out))


if __name__ == "__main__":
    main()
