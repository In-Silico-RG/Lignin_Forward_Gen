#!/usr/bin/env python
"""
calibrate_ref.py

Calibrate the seven coupling constants against the expanded reference set
(data/reference_inventories.csv) with the ratio-to-beta-O-4 objective of
reference_set.py. Each sample is simulated at its own S/G and H fraction; the
chi-square sums over the quantities each sample reports. Design: docs/14.

Same optimiser, bounds, scale (k_bo4 = 8), DP (120), fixed seeds and replicate
counts as calibrate_full.py, so the two fits differ only in the data and the
observable.

    conda run -n lignin python code/calibrate_ref.py [--source-weight] [--reps 8]
        [--hold-out SAMPLE_ID ...]     (leave-one-sample-out folds)

Outputs: results/calibration_ref[<tag>].log and results/fitted_params_ref[<tag>].json
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mini_forward_generator as mfg          # noqa: E402
import reference_set as rs                    # noqa: E402
from calibrate_full import BOUNDS, K_BO4, N, SEED_BASE   # noqa: E402

OUT = os.path.join(HERE, os.pardir, "results")

# start from the adopted v1.0.1 constants (results/fitted_params.json)
X0 = {"k_b5": 2.4904, "k_b1": 0.2715, "k_bb": 16.7825, "k_55": 0.7108,
      "k_4o5": 0.1737, "n_nuclei": 4.716, "k_bo4_nuc": 0.9213}

# Optional eighth parameter: donor-type specificity of chain-extension beta-5
# (mini_forward_generator "b5_s_donor"). Tested like alpha-O-4 and oligo-beta-beta
# were: the eight-parameter model nests the seven-parameter one (b5_s_donor = 1), so
# start AT the seven-parameter optimum (--start) and ask whether the fit takes it.
B5_DONOR_BOUND = ("b5_s_donor", 0.01, 10.0)
# The two eighth parameters already tested at v1.0.1 (SI S7, S11), redone here on the
# expanded set at the adopted constants: the alpha-O-4 chain-joining ether and the
# re-radicalised oligomer beta-beta. Same bounds as calibrate_full.py.
AO4_BOUND = ("k_ao4", 0.005, 40.0)
OLIGOBB_BOUND = ("oligo_bb", 1e-4, 5.0)

# Generator parameters held fixed for the whole fit (ablations): --fix key=value.
FIXED = {}


def unpack(x):
    params = {"k_bo4": K_BO4}
    for (name, lo, hi), v in zip(BOUNDS, np.exp(x)):
        params[name] = float(np.clip(v, lo, hi))
    return params


def simulate(row, params, reps, seed_base=SEED_BASE):
    inp = rs.sim_inputs(row)
    p = dict(params)
    p.update(FIXED)
    p["pct_h"] = inp["pct_h"]
    gens = []
    # seed offset by S/G as in calibrate_full (int(100*sg)) plus the H fraction so
    # two samples with equal S/G but different H do not share draws
    off = int(round(100 * inp["sg"])) + int(round(1000 * inp["pct_h"]))
    for r in range(reps):
        g = mfg.MiniForwardKMC(inp["sg"], n_target=N, params=p, seed=seed_base + r + off)
        g.run()
        gens.append(g)
    return rs.model_observables(gens)


def total_chi2(params, rows, weights, reps):
    tot, per = 0.0, {}
    for row, w in zip(rows, weights):
        obs = simulate(row, params, reps)
        c, terms = rs.chi2_sample(row, obs)
        per[row["sample_id"]] = (c, terms, obs)
        tot += w * c
    return tot, per


def source_weights(rows, by_source):
    # the needs of the many (Happs' 22 samples) are made equal to the needs
    # of the few (the anchors) -- no one source outvotes the rest
    if not by_source:
        return [1.0] * len(rows)
    n = {}
    for r in rows:
        n[r["source"]] = n.get(r["source"], 0) + 1
    return [1.0 / n[r["source"]] for r in rows]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--eval-reps", type=int, default=30)
    ap.add_argument("--source-weight", action="store_true",
                    help="give each source unit total weight instead of each sample")
    ap.add_argument("--hold-out", nargs="*", default=[],
                    help="sample_ids excluded from the fit and predicted afterwards")
    ap.add_argument("--tag", default="")
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--start", default=None,
                    help="fitted_params_*.json whose constants replace the default start")
    ap.add_argument("--with-b5-donor", action="store_true",
                    help="fit b5_s_donor as an eighth parameter (donor-type beta-5 specificity)")
    ap.add_argument("--with-ao4", action="store_true",
                    help="fit k_ao4 as an eighth parameter (alpha-O-4 chain-joining ether)")
    ap.add_argument("--with-oligobb", action="store_true",
                    help="fit oligo_bb as an eighth parameter (re-radicalised oligomer beta-beta)")
    ap.add_argument("--fix", action="append", default=[], metavar="KEY=VALUE",
                    help="hold a generator parameter fixed for the whole fit, e.g. "
                         "--fix end_restricted=0 or --fix nucleation_on=0 --fix pool_target=2")
    args = ap.parse_args()
    if args.start:
        with open(args.start) as fh:
            X0.update({k: v for k, v in json.load(fh).items()
                       if k in X0 or k in ("b5_s_donor", "k_ao4", "oligo_bb")})
    if args.with_b5_donor:
        BOUNDS.append(B5_DONOR_BOUND)
        X0.setdefault("b5_s_donor", 1.0)
    if args.with_ao4:
        BOUNDS.append(AO4_BOUND)
        X0.setdefault("k_ao4", 0.01)          # nearly off: the nested model starts at the 7-param optimum
    if args.with_oligobb:
        BOUNDS.append(OLIGOBB_BOUND)
        X0.setdefault("oligo_bb", 1e-4)
    for kv in args.fix:
        k, v = kv.split("=")
        FIXED[k] = float(v)

    rows_all = rs.load_reference_set()
    rows = [r for r in rows_all if r["sample_id"] not in args.hold_out]
    held = [r for r in rows_all if r["sample_id"] in args.hold_out]
    weights = source_weights(rows, args.source_weight)
    tag = args.tag or ("_sw" if args.source_weight else "") + \
        ("_ho-" + "-".join(args.hold_out) if args.hold_out else "")
    os.makedirs(OUT, exist_ok=True)
    log_path = os.path.join(OUT, "calibration_ref%s.log" % tag)
    lines = []

    def log(s=""):
        lines.append(s); print(s); sys.stdout.flush()

    log("Calibration on the expanded reference set (ratios to beta-O-4)")
    log("samples in fit: %d  held out: %s  reps: %d  N: %d  seed base: %d  source-weighted: %s"
        % (len(rows), args.hold_out or "-", args.reps, N, SEED_BASE, args.source_weight))
    log("generator: %s" % os.path.relpath(mfg.__file__, os.path.join(HERE, os.pardir)))
    log("free parameters: %s  start: %s" % ([b[0] for b in BOUNDS], {k: round(v, 4) for k, v in X0.items()}))
    if FIXED:
        log("held fixed (ablation): %s" % FIXED)

    state = {"n": 0, "best": np.inf}

    def loss(x):
        L, _ = total_chi2(unpack(x), rows, weights, args.reps)
        state["n"] += 1
        state["best"] = min(state["best"], L)
        if state["n"] % 20 == 0:
            log("  eval %3d  chi2=%.3f  best=%.3f" % (state["n"], L, state["best"]))
        return L

    x0 = np.log([X0[name] for name, _, _ in BOUNDS])
    simplex = [x0] + [x0 + np.eye(len(x0))[i] * 0.4 for i in range(len(x0))]
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": args.maxiter, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    best = unpack(res.x)
    log("\nfitted params: %s" % {k: round(v, 4) for k, v in best.items()})
    log("minimised chi2 = %.3f (%d evals)" % (res.fun, state["n"]))

    chi2, per = total_chi2(best, rows, weights, args.eval_reps)
    log("\nchi2 at %d reps: %.3f over %d samples, %d terms"
        % (args.eval_reps, chi2, len(rows), sum(len(t) for _, t, _ in per.values())))
    log("%-28s %5s %6s  %s" % ("sample", "S/G", "chi2", "z per ratio (model - measured)/sigma"))
    for row in rows:
        c, terms, obs = per[row["sample_id"]]
        log("%-28s %5.2f %6.2f  %s" % (row["sample_id"], row["sg_ratio"], c,
                                      " ".join("%s:%+.2f" % (k, z) for k, z in terms.items())))
    if held:
        log("\nheld-out samples (predicted with the fitted constants):")
        for row in held:
            obs = simulate(row, best, args.eval_reps)
            c, terms = rs.chi2_sample(row, obs)
            log("%-28s %5.2f %6.2f  %s" % (row["sample_id"], row["sg_ratio"], c,
                                          " ".join("%s:%+.2f" % (k, z) for k, z in terms.items())))
    log("\nmodel per 100 units at the fitted constants (30 reps):")
    log("%-28s %5s %6s %5s %5s %5s %5s %5s %5s" % ("sample", "S/G", "bo4", "b5", "bb", "b1", "5-5", "4O5", "DC"))
    for row in rows + held:
        o = per[row["sample_id"]][2] if row["sample_id"] in per else simulate(row, best, args.eval_reps)
        log("%-28s %5.2f %6.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f" % (
            row["sample_id"], row["sg_ratio"], o["bo4"], o["b5"], o["bb"], o["b1"], o["c55"], o["c4o5"], o["dc"]))

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, "fitted_params_ref%s.json" % tag), "w") as fh:
        json.dump({"_note": "ratio-to-beta-O-4 fit on data/reference_inventories.csv; see calibration_ref%s.log" % tag,
                   "_held_out": args.hold_out, "_source_weighted": args.source_weight,
                   "_chi2_eval": chi2, **best}, fh, indent=1)
    print("Wrote", os.path.relpath(log_path))


if __name__ == "__main__":
    main()
