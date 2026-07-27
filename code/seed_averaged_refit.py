#!/usr/bin/env python
"""
seed_averaged_refit.py  --  fit the expected inventory, not one noise realisation.

What param_recovery_experiment establishes: the ADOPTED parameters generalise across
seeds (chi2 27+/-3), but the FITTING PROCEDURE itself is noisy -- it minimises a chi2
evaluated at a SINGLE seed base (7000) with only FIT_REPS=8 KMC runs per anchor. A
stochastic objective built from 8 samples has real Monte-Carlo variance, so the
located optimum is partly tuned to seed 7000's particular noise realisation. That is
a property of the PROCEDURE, not the model, and it is worth fixing.

The fix is the statistically correct objective: fit the model's EXPECTED inventory,
estimated by averaging over several DECORRELATED seed streams (M seed bases x reps
each), instead of one. This drives the Monte-Carlo variance of the objective down by
~sqrt(M) and removes any tuning to a single realisation.

This script:
  1. Refits all 7 constants under the seed-AVERAGED chi2 objective.
  2. Compares adopted vs seed-averaged params HEAD-TO-HEAD on HELD-OUT seed bases
     that were never used during either fit -- the only fair test of generalisation.
  3. Decides, on held-out chi2, whether to ADOPT the seed-averaged params.

Nothing is written to fitted_params.json by this script; it prints a recommendation.
Run:  conda run -n lignin python code/seed_averaged_refit.py
Outputs: results/seed_averaged_refit.log, results/seed_averaged_refit.csv
"""
import os
import json
import importlib.util
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("cf", os.path.join(HERE, "calibrate_full.py"))
cf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cf)

PARAM_ORDER = ["k_b5", "k_b1", "k_bb", "k_55", "k_4o5", "n_nuclei", "k_bo4_nuc"]

# Seed bases used to BUILD the objective (the fit sees these). Decorrelated streams.
FIT_SEED_BASES = [7000, 17000, 27000, 37000]
FIT_REPS = 8                       # per seed base -> 32 effective reps per anchor
# Seed bases held OUT of the fit entirely -- the fair generalisation test.
HELDOUT_SEED_BASES = [101000, 202000, 303000, 404000, 505000]
EVAL_REPS = 40                     # per held-out base, for a low-noise score


def adopted_params():
    raw = json.load(open(os.path.join(HERE, os.pardir, "results", "fitted_params.json")))
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def avg_inventory(sg, params, seed_bases, reps):
    """Model inventory averaged over several decorrelated seed streams."""
    accum = {k: 0.0 for k in cf.ALL_KEYS}
    for sb in seed_bases:
        cf.SEED_BASE = sb
        m, _ = cf.model_full(sg, params, reps=reps)
        for k in cf.ALL_KEYS:
            accum[k] += m[k]
    n = len(seed_bases)
    return {k: accum[k] / n for k in cf.ALL_KEYS}


def chi2_over(params, seed_bases, reps):
    """chi2 of the seed-averaged expected inventory against the reference."""
    chi2 = 0.0
    for sg in cf.ANCHORS:
        m = avg_inventory(sg, params, seed_bases, reps)
        for k in cf.ALL_KEYS:
            r = cf.REFERENCE[sg][k]
            chi2 += ((m[k] - r) / cf.sigma(r)) ** 2
    return chi2


def run_avg_fit(emit):
    """Nelder-Mead on the seed-AVERAGED chi2, started at the adopted optimum."""
    ap = adopted_params()
    x0 = np.log([ap[p] for p in PARAM_ORDER])

    state = {"n": 0, "best": np.inf}

    def loss(x):
        params = cf.unpack(x)
        L = chi2_over(params, FIT_SEED_BASES, FIT_REPS)
        state["n"] += 1
        if L < state["best"]:
            state["best"] = L
        if state["n"] % 25 == 0:
            emit("  eval %3d  chi2_avg=%.3f  best=%.3f" % (state["n"], L, state["best"]))
        return L

    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy(); pt[i] += 0.3
        simplex.append(pt)
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": 350, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    return cf.unpack(res.x), state["n"]


def heldout_scores(params):
    """Per-held-out-base chi2 (EVAL_REPS each), plus mean/std -- generalisation."""
    scores = []
    for sb in HELDOUT_SEED_BASES:
        cf.SEED_BASE = sb
        c2, _ = cf.scores(params, reps=EVAL_REPS)
        scores.append(c2)
    arr = np.array(scores)
    return arr, arr.mean(), arr.std()


def main():
    out_dir = os.path.join(HERE, os.pardir, "results")
    lines = []

    def emit(s=""):
        print(s); lines.append(s)

    emit("Seed-averaged refit -- fit stability under a noisy objective")
    emit("fit objective: chi2 of inventory AVERAGED over seed bases %s x %d reps (= %d eff.)"
         % (FIT_SEED_BASES, FIT_REPS, len(FIT_SEED_BASES) * FIT_REPS))
    emit("held-out bases (never seen by any fit): %s x %d reps" % (HELDOUT_SEED_BASES, EVAL_REPS))
    emit("")

    ap = adopted_params()
    emit("adopted params (single-seed reps=8 fit): "
         + json.dumps({k: round(ap[k], 4) for k in PARAM_ORDER}))

    emit("")
    emit("Refitting under the seed-averaged objective (started AT the adopted optimum)...")
    new, nev = run_avg_fit(emit)
    emit("seed-averaged params: " + json.dumps({k: round(new[k], 4) for k in PARAM_ORDER}))
    emit("  (%d evals)" % nev)
    emit("")

    # parameter movement, split by identifiability class
    emit("%-11s %10s %10s %8s   %s" % ("param", "adopted", "seed-avg", "ratio", "class"))
    for p in PARAM_ORDER:
        a, b = ap[p], new[p]
        ratio = b / a if a else float("nan")
        cls = "UNDET" if p in {"k_bb", "k_bo4_nuc"} else "det"
        emit("%-11s %10.4f %10.4f %8.2f   %s" % (p, a, b, ratio, cls))
    emit("")
    emit("Expectation: the four determined constants barely move; any large movement")
    emit("should sit in the two UNDET directions (k_bb, k_bo4_nuc), which the data do")
    emit("not pin regardless of seed. That pattern = the crack was Monte-Carlo noise in")
    emit("the objective, not a real relocation of the identified optimum.")
    emit("")

    # ---- head-to-head on HELD-OUT seeds ------------------------------------
    emit("HELD-OUT generalisation (chi2 on bases neither fit ever saw):")
    a_arr, a_mu, a_sd = heldout_scores(ap)
    n_arr, n_mu, n_sd = heldout_scores(new)
    emit("  adopted  : " + " ".join("%6.2f" % c for c in a_arr)
         + "   mean %.2f +/- %.2f" % (a_mu, a_sd))
    emit("  seed-avg : " + " ".join("%6.2f" % c for c in n_arr)
         + "   mean %.2f +/- %.2f" % (n_mu, n_sd))
    emit("")
    better = n_mu < a_mu
    delta = a_mu - n_mu
    emit("held-out mean chi2: %s by %.2f (%.1f%%)"
         % ("seed-avg BETTER" if better else "adopted better/equal",
            abs(delta), 100 * abs(delta) / a_mu))
    emit("")

    # ---- recommendation ----------------------------------------------------
    # Adopt the seed-averaged fit only if it is at least as good held-out AND the
    # determined constants have not been distorted (a fit that "wins" by moving an
    # UNDET direction into a noise pocket is not an improvement).
    det_stable = all(abs(new[p] / ap[p] - 1) < 0.20
                     for p in ["k_b5", "k_55", "k_4o5", "n_nuclei"])
    if better and det_stable:
        rec = ("ADOPT the seed-averaged params: they generalise at least as well AND "
               "the determined constants are unmoved (<20%). The procedure, not just "
               "the result, is now seed-robust.")
    elif not better and det_stable:
        rec = ("KEEP the adopted params: the seed-averaged fit does NOT beat them "
               "held-out, and the determined constants are unmoved -- i.e. the single-"
               "seed optimum was already at the seed-robust location. The crack was "
               "cosmetic: reps=8 added variance to the OBJECTIVE but did not mislocate "
               "the identified optimum. Report the seed-averaged held-out chi2 as the "
               "stability check.")
    else:
        rec = ("INSPECT: a determined constant moved >20%. Report both fits; do not "
               "silently adopt.")
    emit("RECOMMENDATION: " + rec)

    # ---- write artifacts ---------------------------------------------------
    with open(os.path.join(out_dir, "seed_averaged_refit.csv"), "w") as fh:
        fh.write("param,adopted,seed_avg,ratio,class\n")
        for p in PARAM_ORDER:
            a, b = ap[p], new[p]
            cls = "UNDET" if p in {"k_bb", "k_bo4_nuc"} else "det"
            fh.write("%s,%.6f,%.6f,%.4f,%s\n" % (p, a, b, (b / a if a else float("nan")), cls))
        fh.write("\nheldout_base,adopted_chi2,seedavg_chi2\n")
        for sb, ca, cn in zip(HELDOUT_SEED_BASES, a_arr, n_arr):
            fh.write("%d,%.4f,%.4f\n" % (sb, ca, cn))
        fh.write("mean,%.4f,%.4f\n" % (a_mu, n_mu))
        fh.write("std,%.4f,%.4f\n" % (a_sd, n_sd))

    with open(os.path.join(out_dir, "seed_averaged_refit.log"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    emit("")
    emit("wrote results/seed_averaged_refit.csv, results/seed_averaged_refit.log")

    # machine-readable summary for a later adoption decision
    json.dump({
        "adopted": {k: ap[k] for k in PARAM_ORDER},
        "seed_avg": {k: new[k] for k in PARAM_ORDER},
        "heldout_adopted_mean": a_mu, "heldout_adopted_std": a_sd,
        "heldout_seedavg_mean": n_mu, "heldout_seedavg_std": n_sd,
        "determined_stable": det_stable, "seedavg_better": bool(better),
    }, open(os.path.join(out_dir, "seed_averaged_refit_summary.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
