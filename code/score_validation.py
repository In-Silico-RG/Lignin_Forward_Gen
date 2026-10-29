#!/usr/bin/env python
"""
score_validation.py

Predict the validation-only samples of data/reference_inventories.csv
(use_fit == 0, use_val == 1) at a fixed parameter set and report the ratio-to-
beta-O-4 residuals, sample by sample. Nothing is fitted. Same simulation and
chi-square as calibrate_ref.py (simulate, reference_set.chi2_sample).

    python code/score_validation.py --params results/fitted_params_ref_v11d_sw.json
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reference_set as rs                    # noqa: E402
from calibrate_ref import simulate, OUT       # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--params", default=os.path.join(OUT, "fitted_params_ref_v11d_sw.json"))
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--only", nargs="*", default=[], help="sample_id prefixes")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    params = {k: v for k, v in json.load(open(a.params)).items() if not k.startswith("_")}
    rows = [r for r in rs.load_reference_set(fit_only=False)
            if r.get("use_fit", "0").strip() == "0" and r.get("use_val", "0").strip() == "1"]
    if a.only:
        rows = [r for r in rows if r["sample_id"].startswith(tuple(a.only))]

    out = {"params": os.path.basename(a.params), "reps": a.reps, "samples": {}}
    print("%-36s %5s %5s  %7s  %s" % ("sample", "S/G", "H", "chi2/n", "measured -> model ratio (z)"))
    for r in rows:
        obs = simulate(r, params, a.reps)
        c, terms = rs.chi2_sample(r, obs)
        meas, mod = rs.measured_ratios(r), rs.model_ratios(obs)
        cells = " ".join("%s %.3f->%.3f (%+.1f)" % (k, meas[k][0], mod[k], z) for k, z in terms.items())
        print("%-36s %5.2f %5.3f  %7.2f  %s" % (r["sample_id"], r["sg_ratio"], r["h_frac"] or 0.0,
                                                c / max(len(terms), 1), cells))
        out["samples"][r["sample_id"]] = {"chi2": c, "z": terms, "model_per100": obs,
                                          "measured_ratio": {k: v[0] for k, v in meas.items()},
                                          "model_ratio": {k: mod[k] for k in terms}}
    path = os.path.join(OUT, "validation_scores%s.json" % a.tag)
    json.dump(out, open(path, "w"), indent=1)
    print("->", os.path.relpath(path))


if __name__ == "__main__":
    main()
