#!/usr/bin/env python
"""
predict.py

Predict the linkage inventory of a new feedstock from its S/G ratio and H
fraction, using the paper's adopted, fitted constants. This is the entry
point an applied user (screening candidate biomass for valorization) actually
wants -- everything else in code/ is the machinery that produced and tested
those constants.

Uses the exact same simulation code as the paper's own numbers
(calibrate_ref.simulate), so a result here is reproducible against the
committed logs, not a reimplementation.

    conda run -n lignin python code/predict.py --sg 1.4 --h 0.03
    conda run -n lignin python code/predict.py --sg 1.4 --h 0.03 --reps 30
    conda run -n lignin python code/predict.py --sg 1.4 --h 0.03 \
        --params results/fitted_params_ref_v11f_sw.json
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reference_set as rs                    # noqa: E402
from calibrate_ref import simulate             # noqa: E402

DEFAULT_PARAMS = os.path.join(HERE, os.pardir, "results", "fitted_params_ref_v11f_sw.json")


def fitted_range():
    """(sg_min, sg_max, h_min, h_max) across the 60 samples actually fitted."""
    rows = rs.load_reference_set(fit_only=True)
    sgs = [r["sg_ratio"] for r in rows]
    hs = [r["h_frac"] or 0.0 for r in rows]
    return min(sgs), max(sgs), min(hs), max(hs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sg", type=float, required=True, help="S/G ratio of the feedstock")
    ap.add_argument("--h", type=float, default=0.0, help="H (p-hydroxyphenyl) fraction, 0-1 (default 0)")
    ap.add_argument("--params", default=DEFAULT_PARAMS,
                     help="fitted_params*.json to run at (default: the paper's adopted v11f_sw fit)")
    ap.add_argument("--reps", type=int, default=30, help="replicate runs (default 30, matches the paper)")
    args = ap.parse_args()

    with open(args.params) as fh:
        params = {k: v for k, v in json.load(fh).items() if not k.startswith("_")}

    row = {"sg_ratio": args.sg, "h_frac": args.h}
    obs = simulate(row, params, args.reps)
    ratios = rs.model_ratios(obs)

    sg_lo, sg_hi, h_lo, h_hi = fitted_range()
    out_of_range = not (sg_lo <= args.sg <= sg_hi) or not (h_lo <= args.h <= h_hi)

    print("Predicted linkage inventory  (S/G = %.2f, H = %.1f%%, %d reps, %s)"
          % (args.sg, 100 * args.h, args.reps, os.path.basename(args.params)))
    print()
    print("  per 100 units          bo4   b5    bb    b1    55    4o5")
    print("  %-20s %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f"
          % ("", obs["bo4"], obs["b5"], obs["bb"], obs["b1"], obs["c55"], obs["c4o5"]))
    print()
    print("  ratio to beta-O-4       b5    bb    b1    55    4o5    dc")
    print("  %-20s %5.3f %5.3f %5.3f %5.3f %5.3f %5.1f"
          % ("", ratios["b5"], ratios["bb"], ratios["b1"], ratios["c55"], ratios["c4o5"], ratios["dc"]))
    print()

    if out_of_range:
        print("WARNING: S/G=%.2f, H=%.2f is OUTSIDE the range of the 60 fitted samples "
              "(S/G %.1f-%.1f, H %.2f-%.2f)." % (args.sg, args.h, sg_lo, sg_hi, h_lo, h_hi))
        print("The model will still run, but this prediction is an extrapolation,")
        print("not interpolation within validated territory -- treat it accordingly.")
    else:
        print("S/G and H are both within the range of the 60 fitted samples "
              "(S/G %.1f-%.1f, H %.2f-%.2f)." % (sg_lo, sg_hi, h_lo, h_hi))
        print("That means your feedstock's composition is interpolation, not")
        print("extrapolation -- it does not mean a sample this different from")
        print("everything in the fit is guaranteed accurate; see the paper's")
        print("Limitations for what the model does and doesn't claim.")


if __name__ == "__main__":
    main()
