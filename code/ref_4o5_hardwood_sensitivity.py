"""
ref_4o5_hardwood_sensitivity.py

What happens to the fit if the hardwood 4-O-5 reference is too low?

code/sigma_4o5_sensitivity.py asked whether the WEIGHT on 4-O-5 matters and found
it does not: the optimum is unmoved from sigma = 1 to 4 points. This script asks
the harder question, about the VALUE.

The adopted hardwood 4-O-5 references are 1.5 (poplar) and 0.5 (high-S) per 100
C9. Capanema et al. 2005 (JAFC 53:9639) Table 1 collects five determinations of
hardwood lignin, giving 4-O-5 totals of 1.5, 6.0, 6.5, 9.0 and 11.5 per 100 Ar.
Our value is the lowest of the five. Their own figure for E. grandis, 8-9 per 100
Ar, is not a degradative measurement: it is obtained from quantitative 13C NMR by
difference (etherified guaiacyl minus the 148-144.5 ppm integral, 0.25 - 0.17 =
0.08/Ar) and then cross-checked against permanganate oxidation. So the standing
argument -- that high 4-O-5 figures come from degradative methods that over-count
condensed dimers ~8-fold -- does not reach the strongest counter-example.

The existing --ref-4o5 flag on calibrate_full overrides the SOFTWOOD anchor only.
The disagreement documented above is about hardwoods, so this script scales the
two hardwood 4-O-5 references together, preserving the direction of the S/G trend
while raising its level, and refits all seven constants at each scale.

    scale 1x : poplar 1.5, high-S 0.5   (adopted)
    scale 2x : poplar 3.0, high-S 1.0
    scale 4x : poplar 6.0, high-S 2.0   (poplar at Capanema's E. grandis value)

What to look for: if the fitted constants and the other five residuals barely
move, the choice of 4-O-5 reference is a local matter confined to one linkage. If
they move, the reference uncertainty propagates into everything and the headline
error is resting on a disputed number.

Usage:
    python code/ref_4o5_hardwood_sensitivity.py [--reps 8] [--report-reps 30]
Writes results/ref_4o5_hardwood_sensitivity.csv and prints a report.
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

import calibrate_full as cf  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import rel

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

SCALES = [1.0, 2.0, 4.0]
HARDWOODS = [1.0, 3.0]          # S/G of the two hardwood anchors
#         k_b5    k_b1    k_bb     k_55    k_4o5   n_nuclei  k_bo4_nuc
BASE_X0 = [2.4904, 0.2715, 16.7825, 0.7108, 0.1737, 4.716,    0.9213]


def reference_at(scale):
    """Adopted raw reference with both hardwood 4-O-5 values scaled."""
    raw = {sg: dict(d) for sg, d in cf.RAW_REFERENCE.items()}
    for sg in HARDWOODS:
        raw[sg]["4o5"] = raw[sg]["4o5"] * scale
    return cf.build_reference(raw)


def chi2_against(params, ref, reps):
    total = 0.0
    for sg in cf.ANCHORS:
        m, _ = cf.model_full(sg, params, reps=reps)
        for k in cf.ALL_KEYS:
            r = ref[sg][k]
            total += ((m[k] - r) / cf.sigma(r)) ** 2
    return total


def fit_against(ref, reps):
    def loss(x):
        return chi2_against(cf.unpack(x), ref, reps)

    x0 = np.log(BASE_X0)
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy()
        pt[i] += 0.4
        simplex.append(pt)
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": 400, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    return cf.unpack(res.x), res.fun


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=cf.FIT_REPS)
    ap.add_argument("--report-reps", type=int, default=30)
    args = ap.parse_args()

    print("hardwood 4-O-5 reference sensitivity")
    print("  adopted hardwood 4-O-5 (per 100 C9): poplar %.1f, high-S %.1f"
          % (cf.RAW_REFERENCE[1.0]["4o5"], cf.RAW_REFERENCE[3.0]["4o5"]))
    print("  Capanema 2005 Table 1, five hardwood determinations, per 100 Ar:")
    print("    4-O-5 total = 1.5, 6.0, 6.5, 9.0, 11.5   (ours is the lowest)\n")

    rows, fits = [], {}
    for s in SCALES:
        ref = reference_at(s)
        params, _ = fit_against(ref, args.reps)
        per_wood = {}
        for sg in cf.ANCHORS:
            m, _ = cf.model_full(sg, params, reps=args.report_reps)
            res = {k: 100 * (m[k] - ref[sg][k]) for k in cf.ALL_KEYS}
            rmse = math.sqrt(sum(v ** 2 for v in res.values()) / len(cf.ALL_KEYS))
            per_wood[sg] = (rmse, res, m)
            rows.append({"scale": s, "sg": sg, "rmse_pts": round(rmse, 4),
                         "ref_4o5_pct": round(100 * ref[sg]["4o5"], 3),
                         **{"resid_" + k: round(res[k], 3) for k in cf.ALL_KEYS},
                         **{"param_" + k: round(v, 5)
                            for k, v in params.items() if k != "k_bo4"}})
        mean_rmse = sum(w[0] for w in per_wood.values()) / len(per_wood)
        fits[s] = (params, per_wood, mean_rmse)
        tag = "  (adopted)" if s == 1.0 else ""
        print("scale %.0fx%s   poplar ref %.1f, high-S ref %.1f per 100 C9"
              % (s, tag, cf.RAW_REFERENCE[1.0]["4o5"] * s,
                 cf.RAW_REFERENCE[3.0]["4o5"] * s))
        print("   mean RMSE %.3f pts   params %s"
              % (mean_rmse, {k: round(v, 4) for k, v in params.items()
                             if k != "k_bo4"}))
        for sg in cf.ANCHORS:
            rmse, res, m = per_wood[sg]
            print("   S/G %.1f  RMSE %5.2f  4-O-5 %+6.2f  5-5 %+6.2f  bO4 %+6.2f"
                  % (sg, rmse, res["4o5"], res["55"], res["bo4"]))
        print()

    base = fits[1.0]
    print("CHANGE FROM THE ADOPTED REFERENCE")
    print("%-8s %10s %13s %13s %14s"
          % ("scale", "mean RMSE", "d(mean RMSE)", "max |dk|/k", "mean 4-O-5 res"))
    for s in SCALES:
        params, per_wood, mean_rmse = fits[s]
        rel = max(abs(params[k] - base[0][k]) / base[0][k]
                  for k in params if k != "k_bo4")
        m4 = sum(per_wood[sg][1]["4o5"] for sg in cf.ANCHORS) / len(cf.ANCHORS)
        print("%-8.0f %10.3f %+13.3f %12.1f%% %+14.2f"
              % (s, mean_rmse, mean_rmse - base[2], 100 * rel, m4))

    out = os.path.join(RESULTS, "ref_4o5_hardwood_sensitivity.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote %s" % rel(out))


if __name__ == "__main__":
    main()
