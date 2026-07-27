"""
hardwood_sg_trends.py

Does the S/G dependence hold up on hardwoods the model never saw?

The calibration uses three complete distributions, one of them a single
high-syringyl hardwood. Balakshin et al. (2016, Holzforschung 70:95-108, Table 4)
measured 12 industrially used hardwood species by quantitative 13C NMR, spanning
S/G 1.15-3.02. None entered the calibration, so the trends they report are an
outside check on the ones the model produces.

What the set CAN test: the direction and strength of the S/G dependence of
beta-O-4, beta-beta, beta-5, and of condensation as a whole.

What it CANNOT test: the six-linkage distribution. They report a single degree of
condensation (DC) and do not resolve 5-5 from 4-O-5, for the same reason everyone
else does not. Their numbers are per 100 aromatic rings; ours are percentages of
the six modelled linkages. Those normalisations differ, so this script compares
the SIGN and the STRENGTH of each trend, never the magnitude of a slope against
a slope.

Usage:
    python code/hardwood_sg_trends.py
Writes results/hardwood_sg_trends.csv and prints a report.
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data")
RESULTS = os.path.join(HERE, os.pardir, "results")

# The calibrated model's own S/G sweep, % of the six modelled linkages.
# Source: results/mini_forward_full_calibrated.csv (code/calibrate_full.py).
MODEL_SWEEP = os.path.join(RESULTS, "mini_forward_full_calibrated.csv")

# Reference anchors, per 100 C9, before normalising (calibrate_full.RAW_REFERENCE).
REF_RAW = {0.0: 46.0, 1.0: 62.0, 3.0: 78.0}


def load_balakshin():
    rows = []
    with open(os.path.join(DATA, "balakshin2016_hardwoods.csv")) as fh:
        for r in csv.DictReader(fh):
            rows.append({"species": r["species"], "sg": float(r["sg_ratio"]),
                         "bo4": float(r["bo4_total"]), "bb": float(r["bb"]),
                         "b5": float(r["b5"]), "dc": float(r["dc"])})
    return rows


def load_model():
    sg, cols = [], {k: [] for k in ("bo4", "b5", "bb", "b1", "55", "4o5")}
    with open(MODEL_SWEEP) as fh:
        for r in csv.DictReader(fh):
            sg.append(float(r["sg_ratio"]))
            for k in cols:
                cols[k].append(100 * float(r[k]))
    return np.array(sg), {k: np.array(v) for k, v in cols.items()}


def trend(x, y):
    slope = np.polyfit(x, y, 1)[0]
    r = float(np.corrcoef(x, y)[0, 1])
    return slope, r, r * r


def main():
    rows = load_balakshin()
    sg = np.array([r["sg"] for r in rows])
    msg, model = load_model()
    # restrict the model sweep to the range the hardwood set actually covers
    lo, hi = sg.min(), sg.max()
    sel = (msg >= lo) & (msg <= hi)
    cond_model = model["55"] + model["4o5"]

    print("Hardwood S/G trends: 12 species (Balakshin 2016) vs the calibrated model")
    print("  measured range S/G %.2f-%.2f, n=%d, quantitative 13C NMR, per 100 Ar"
          % (lo, hi, len(rows)))
    print("  model compared over the SAME S/G range, %% of the six modelled linkages")
    print("  NOTE normalisations differ; compare sign and strength, not magnitude.\n")

    print("%-14s %22s %26s" % ("", "measured (12 species)", "model (same S/G range)"))
    print("%-14s %9s %7s %5s   %9s %7s" %
          ("quantity", "slope", "r", "r2", "slope", "direction"))
    print("-" * 66)

    out = []
    pairs = [("beta-O-4", "bo4", model["bo4"]),
             ("beta-beta", "bb", model["bb"]),
             ("beta-5", "b5", model["b5"]),
             ("condensation", "dc", cond_model)]
    for label, key, mvals in pairs:
        y = np.array([r[key] for r in rows])
        s, r, r2 = trend(sg, y)
        ms, _, _ = trend(msg[sel], mvals[sel])
        agree = "same sign" if np.sign(s) == np.sign(ms) else "OPPOSITE"
        print("%-14s %+9.2f %+7.2f %5.2f   %+9.2f %7s"
              % (label, s, r, r2, ms, agree))
        out.append({"quantity": label, "meas_slope": round(s, 3),
                    "meas_r": round(r, 3), "meas_r2": round(r2, 3),
                    "meas_min": y.min(), "meas_max": y.max(),
                    "model_slope": round(ms, 3), "agreement": agree})

    print("\nThe two quantities that move with S/G in the measured set are beta-O-4")
    print("(up) and condensation (down) -- the pair the S-blocks-C5 mechanism produces.")
    print("beta-beta does not move, which is also what the model does.")

    # The high-syringyl anchor against the measured species nearest it in S/G.
    print("\nHIGH-SYRINGYL ANCHOR")
    top = max(rows, key=lambda r: r["sg"])
    print("  calibration anchor: S/G 3.0, beta-O-4 %.0f per 100 C9" % REF_RAW[3.0])
    print("  nearest measured  : %s, S/G %.2f, beta-O-4 %.0f per 100 Ar"
          % (top["species"], top["sg"], top["bo4"]))
    print("  anchor exceeds it by %+.0f" % (REF_RAW[3.0] - top["bo4"]))
    hi_sg = [r for r in rows if r["sg"] >= 2.5]
    print("  all %d species with S/G >= 2.5: beta-O-4 %s (mean %.1f)"
          % (len(hi_sg), [int(r["bo4"]) for r in hi_sg],
             np.mean([r["bo4"] for r in hi_sg])))
    print("  The model predicts UNDER this anchor, so part of that residual may be")
    print("  the anchor rather than the model.")

    path = os.path.join(RESULTS, "hardwood_sg_trends.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print("\nWrote %s" % os.path.normpath(path))


if __name__ == "__main__":
    main()
