#!/usr/bin/env python
"""
identifiability_ref.py  (rerun R3, docs/17)

Local identifiability of the seven constants on the expanded-set objective. For each
constant, sweep it across its full search range with the others at the adopted values
and measure how far the vector of model observables moves, in units of the model's own
run-to-run noise. The observables are exactly the terms of the objective: the ratios to
beta-O-4 that each of the 60 samples reports (171 components).

  noise sigma_j : across-seed standard deviation of component j at the adopted constants
                  (N_BASES seed bases, FIT_REPS replicates each, i.e. the fitting setup)
  movement      : max_j |r_j(theta) - r_j(theta_hat)| / sigma_j over the sweep grid, and
                  the RMS over components at the extreme grid points

Same reading as SI S5 of the v1.0.1 paper: a constant whose sweep moves the observables
by many sigma is determined by the data; one that stays within a few sigma is not.

    python code/identifiability_ref.py --params results/fitted_params_ref_v60h_sw.json
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reference_set as rs                                  # noqa: E402
import calibrate_ref as cr                                  # noqa: E402
from calibrate_full import BOUNDS, K_BO4, SEED_BASE          # noqa: E402

OUT = os.path.join(HERE, os.pardir, "results")
N_GRID = 9
N_BASES = 16


def load_params(path):
    with open(path) as fh:
        d = json.load(fh)
    p = {k: v for k, v in d.items() if not k.startswith("_")}
    p.setdefault("k_bo4", K_BO4)
    return p


def observable_vector(rows, params, reps, seed_base):
    """Concatenated model ratios, in the order of the samples' measured keys."""
    vec, labels = [], []
    for row in rows:
        obs = cr.simulate(row, params, reps, seed_base=seed_base)
        mr = rs.model_ratios(obs)
        for k in rs.measured_ratios(row):
            vec.append(mr[k])
            labels.append((row["sample_id"], k))
    return np.array(vec), labels


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--params", required=True)
    ap.add_argument("--reps", type=int, default=cr.__dict__.get("FIT_REPS", 8))
    ap.add_argument("--grid", type=int, default=N_GRID)
    ap.add_argument("--bases", type=int, default=N_BASES)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    rows = rs.load_reference_set()
    theta = load_params(args.params)
    lines = []

    def log(s=""):
        lines.append(s); print(s); sys.stdout.flush()

    log("Identifiability sweep on the expanded-set objective")
    log("params: %s" % os.path.relpath(args.params))
    log("adopted: %s" % {k: round(v, 4) for k, v in theta.items()})
    log("samples: %d  reps per evaluation: %d  seed bases for noise: %d  grid points: %d"
        % (len(rows), args.reps, args.bases, args.grid))

    # --- noise: across seed bases at the adopted constants ---------------------------
    vecs = []
    for i in range(args.bases):
        v, labels = observable_vector(rows, theta, args.reps, SEED_BASE + 1000 * i)
        vecs.append(v)
    vecs = np.array(vecs)
    r_hat = vecs[0]                                  # base 7000 = the fitting seeds
    sigma = vecs.std(axis=0, ddof=1)
    sigma = np.where(sigma > 0, sigma, np.nanmin(sigma[sigma > 0]))
    log("\nnoise: %d components; median sigma/|r| = %.3f" % (len(sigma), np.median(sigma / np.maximum(np.abs(r_hat), 1e-9))))

    # --- sweep ------------------------------------------------------------------------
    log("\n%-14s %-16s %10s %10s %10s   verdict" % ("constant", "range", "max|d|/s", "rms lo", "rms hi"))
    rows_out = []
    for name, lo, hi in BOUNDS:
        grid = np.exp(np.linspace(np.log(lo), np.log(hi), args.grid))
        maxd, rms_lo, rms_hi = 0.0, None, None
        per_point = []
        for j, val in enumerate(grid):
            p = dict(theta); p[name] = float(val)
            v, _ = observable_vector(rows, p, args.reps, SEED_BASE)
            d = (v - r_hat) / sigma
            m = float(np.max(np.abs(d))); r = float(np.sqrt(np.mean(d ** 2)))
            per_point.append((float(val), m, r))
            maxd = max(maxd, m)
            if j == 0:
                rms_lo = r
            if j == args.grid - 1:
                rms_hi = r
        verdict = "determined" if maxd > 5 else "weak" if maxd > 2.5 else "UNDETERMINED"
        log("%-14s [%7.3g, %6.3g] %10.1f %10.1f %10.1f   %s" % (name, lo, hi, maxd, rms_lo, rms_hi, verdict))
        for val, m, r in per_point:
            rows_out.append({"constant": name, "value": val, "max_abs_d": m, "rms_d": r})

    log("\nReading: max|d|/sigma is the largest movement of any observable, in units of its")
    log("own run-to-run noise, anywhere in the constant's search range; rms lo/hi are the")
    log("root-mean-square movements at the two ends of the range.")

    tag = args.tag
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    with open(os.path.join(OUT, "identifiability_ref%s.log" % tag), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, "identifiability_ref%s.csv" % tag), "w") as fh:
        fh.write("constant,value,max_abs_d,rms_d\n")
        for r in rows_out:
            fh.write("%s,%.6g,%.4f,%.4f\n" % (r["constant"], r["value"], r["max_abs_d"], r["rms_d"]))
    print("Wrote results/identifiability_ref%s.log" % tag)


if __name__ == "__main__":
    main()
