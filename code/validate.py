#!/usr/bin/env python
"""
validate.py

Leave-one-out cross-validation of the forward generator on NATIVE lignin anchors.

We deliberately do NOT validate against industrial pulping lignins (Kraft/Soda/
Organosolv): pulping cleaves beta-O-4 and condenses the polymer, so those samples
are not native lignin and are invalid comparison data for a native-biosynthesis
generator. Independent native full-inventory NMR data with known S/G would be the
ideal external test set; lacking that here, we test generalization internally.

Leave-one-out: for each of the three native wood anchors (softwood S/G 0, poplar 1,
high-S 3), we RE-FIT the coupling constants on the other two and PREDICT the
held-out anchor. If the held-out inventory is predicted well from a fit that never
saw it, the model generalizes rather than merely interpolating fitted points.

Run: conda run -n lignin python code/validate.py
Outputs: results/validation_loocv.log, results/validation_loocv.csv
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

# One source of truth for which objective is primary -- imported, not redefined, so
# the two scripts cannot drift apart and silently validate a fit we do not report.
_cf_spec = importlib.util.spec_from_file_location(
    "cf", os.path.join(HERE, "calibrate_full.py"))
_cf = importlib.util.module_from_spec(_cf_spec)
_cf_spec.loader.exec_module(_cf)
PRIMARY_OBJECTIVE = _cf.PRIMARY_OBJECTIVE

ALL_KEYS = ["bo4", "b5", "bb", "b1", "55", "4o5"]


def _norm(d):
    t = sum(d.values())
    return {k: v / t for k, v in d.items()}


# native literature inventories, % of total linkages (see docs/02)
REFERENCE = {
    # revised 2026-07-16; see calibrate_full.REFERENCE and docs/08 for sourcing
    0.0: _norm({"bo4": 46.0, "b5": 11.0, "bb": 2.5, "b1": 2.0, "55": 4.0, "4o5": 1.0}),
    1.0: _norm({"bo4": 62.0, "b5": 6.0,  "bb": 4.0, "b1": 1.5, "55": 3.0,  "4o5": 1.5}),
    3.0: _norm({"bo4": 78.0, "b5": 3.0,  "bb": 4.0, "b1": 1.0, "55": 1.0,  "4o5": 0.5}),
}
ANCHORS = [0.0, 1.0, 3.0]
FIT_REPS = 8
N = 120
K_BO4 = 8.0

# Must mirror calibrate_full's uncertainty model exactly. See the long comment
# there: sigma_k = max(5% relative, 1 pt floor), the 5% from Bourmaud's HSQC0
# beta-O-4 uncertainties (docs/02 section B), the floor because no minor
# cross-peak integrates better than ~1 pt.
REL_SIGMA = 0.05
FLOOR_SIGMA = 0.01

OBJECTIVE = "sse"     # set by --objective; must match how calibrate_full was run


def sigma(ref_frac):
    return max(REL_SIGMA * ref_frac, FLOOR_SIGMA)
# Must mirror calibrate_full.BOUNDS exactly, or the held-out test would be
# validating a different model from the one we calibrate and report.
BOUNDS = [("k_b5", 0.05, 40.0), ("k_b1", 0.01, 20.0), ("k_bb", 0.05, 60.0),
          ("k_55", 0.005, 40.0), ("k_4o5", 0.005, 40.0), ("n_nuclei", 1.0, 30.0),
          ("k_bo4_nuc", 0.02, 40.0)]


def full_fracs(bonds):
    t = sum(bonds[k] for k in ALL_KEYS)
    return {k: (bonds[k] / t if t else 0.0) for k in ALL_KEYS}


def model_full(sg, params, reps=FIT_REPS):
    agg = {k: 0 for k in ALL_KEYS}
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=N, params=params, seed=7000 + r + int(sg * 100))
        g.run()
        for k in ALL_KEYS:
            agg[k] += g.bonds[k]
    return full_fracs(agg)


def unpack(x):
    vals = np.exp(x)
    params = {"k_bo4": K_BO4}
    for (name, lo, hi), v in zip(BOUNDS, vals):
        params[name] = float(np.clip(v, lo, hi))
    return params


def fit(train_anchors):
    #             k_b5  k_b1  k_bb  k_55  k_4o5 n_nuclei k_bo4_nuc
    x0 = np.log([2.87, 0.28, 5.86, 2.09, 0.30, 5.29, 1.09])
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy(); pt[i] += 0.4; simplex.append(pt)

    def loss(x):
        # Must mirror calibrate_full's objective, or the held-out test would be
        # validating a different fit from the one we report.
        p = unpack(x)
        L = 0.0
        for sg in train_anchors:
            m = model_full(sg, p)
            for k in ALL_KEYS:
                r = REFERENCE[sg][k]
                d = m[k] - r
                L += (d / sigma(r)) ** 2 if OBJECTIVE == "chi2" else d ** 2
        return L

    res = minimize(loss, x0, method="Nelder-Mead",
                   # fatol follows the objective's scale (sse ~0.01, chi2 ~10-100)
                   options={"maxiter": 140, "xatol": 1e-2,
                            "fatol": 1e-3 if OBJECTIVE == "chi2" else 1e-6,
                            "initial_simplex": np.array(simplex)})
    return unpack(res.x)


def main():
    out_dir = os.path.join(HERE, os.pardir, "results")
    # Suffix by objective, as calibrate_full.py does. Without this both objectives
    # wrote to ONE file and whichever ran last silently won -- the CSV then
    # misrepresented which objective produced it, with nothing to reveal the
    # substitution. Found 2026-07-16 by running --objective chi2 over an sse result.
    suffix = "" if OBJECTIVE == PRIMARY_OBJECTIVE else "_" + OBJECTIVE
    csv = os.path.join(out_dir, "validation_loocv%s.csv" % suffix)
    print("Leave-one-out cross-validation on native wood anchors")
    print("(no pulping/industrial lignin used; native inventories only)\n")

    rows = []
    for held in ANCHORS:
        train = [a for a in ANCHORS if a != held]
        params = fit(train)
        pred = model_full(held, params, reps=40)
        ref = REFERENCE[held]
        errs = [(100 * pred[k] - 100 * ref[k]) ** 2 for k in ALL_KEYS]
        rmse = np.sqrt(np.mean(errs))
        rows.append((held, train, pred, ref, rmse))
        print("HELD-OUT S/G=%.1f (fit on %s):" % (held, train))
        print("  %-7s %-11s %-11s %-11s %-11s %-11s %-11s" %
              ("", *ALL_KEYS))
        print("  pred    " + " ".join("%-11.1f" % (100 * pred[k]) for k in ALL_KEYS))
        print("  ref     " + " ".join("%-11.1f" % (100 * ref[k]) for k in ALL_KEYS))
        print("  held-out RMSE: %.1f pts\n" % rmse)

    overall = np.sqrt(np.mean([r[4] ** 2 for r in rows]))
    print("mean held-out RMSE across folds: %.1f pts" % overall)
    with open(csv, "w") as fh:
        fh.write("held_sg,type,predicted,reference\n")
        for held, train, pred, ref, rmse in rows:
            for k in ALL_KEYS:
                fh.write("%.1f,%s,%.2f,%.2f\n" % (held, k, 100 * pred[k], 100 * ref[k]))
    print("Wrote %s" % rel(csv))


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--objective", choices=["sse", "chi2"],
                    default=PRIMARY_OBJECTIVE,
                    help="must match how calibrate_full.py was run, or the "
                         "held-out error would describe a fit we do not report. "
                         "Default follows calibrate_full.PRIMARY_OBJECTIVE, chi2 "
                         "since 2026-07-16.")
    OBJECTIVE = ap.parse_args().objective
    main()
