#!/usr/bin/env python
"""
ao4_multifragment_experiment.py

Rev_Kimi_1 concern 2.4: the alpha-O-4 8th parameter was rejected (k_ao4 -> ~0, no chi2
gain) at the calibrated single-molecule regime. The reviewer's hypothesis is that this
is an artifact of the single-molecule restriction -- "the model is over-determined in
its connectivity objective" -- and that alpha-O-4 might earn its place if the model
were allowed multiple fragments. This tests that directly.

alpha-O-4 in this model is a WELDER: a quinone methide (left by a beta-O-4 donor) is
captured by a phenolic O4 on a DIFFERENT fragment, joining the two by an ether without
spending a condensed C-C bond. So it competes with the 5-5 / 4-O-5 channel for the
job of joining fragments -- but it is NOT one of the six measured linkages, and joining
two fragments does not change any of the six inventory FRACTIONS. That is the crux.

DESIGN. For each connectivity regime -- oo_scale from 1 (calibrated single molecule)
down to 0 (multi-fragment, since the connectivity experiment showed welding is bistable:
n_frags = 1 for any oo_scale > 0, and n_frags = n_nuclei only at oo_scale = 0) -- we
REFIT the seven constants PLUS k_ao4 (8 params), holding oo_scale fixed. Letting the
other constants (notably k_55) compensate for reduced welding is what makes this a fair
test rather than the held-fixed sweep. We then read, per regime:
  - best chi2 (does any regime beat the single-molecule baseline?)
  - the adopted k_ao4 (is alpha-O-4 taken up when fragments exist?)
  - n_frags (did we actually reach a multi-fragment regime?)
  - the 4-O-5 residual and the 5-5 / 4-O-5 levels (does alpha-O-4 relieve the 4-O-5
    over-prediction, or does killing the condensed channel just lose 5-5 and 4-O-5?)

The expected reading, if alpha-O-4's rejection is robust (not a single-molecule
artifact): every regime that keeps the condensed inventory (oo_scale > 0) stays one
molecule with k_ao4 rejected; the only genuinely multi-fragment regime (oo_scale = 0)
loses 5-5 and 4-O-5 entirely -- which alpha-O-4 cannot replace, because it is not a
condensed linkage -- so chi2 collapses there regardless of k_ao4. alpha-O-4 is then
unidentified by the linkage inventory in EVERY regime, single-molecule or not.

Usage:
    python code/ao4_multifragment_experiment.py [--fit-reps 8] [--eval-reps 30]
Writes results/ao4_multifragment_experiment.csv and .log.
"""
import argparse
import os

import numpy as np
from scipy.optimize import minimize

import mini_forward_generator as mfg
import calibrate_full as cf

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

OO_REGIMES = [1.0, 0.1, 0.01, 0.0]
# calibrated 7-param chi2 optimum, plus k_ao4 seeded near its lower bound so the
# nested 8-param model can only improve (same nesting logic as calibrate_full).
BASE_X0 = [2.4904, 0.2715, 16.7825, 0.7108, 0.1737, 4.716, 0.9213, 0.01]


def n_frags_at(sg, params, oo_scale, reps=30):
    p = dict(params)
    p["oo_scale"] = oo_scale
    fr = []
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=cf.N, params=p, seed=7000 + r + int(sg * 100))
        s = g.run()
        fr.append(s["n_frags"])
    return float(np.mean(fr))


def refit_at(oo_scale, fit_reps, eval_reps):
    """Refit [7 constants + k_ao4] with oo_scale held. Returns a result dict."""
    # configure calibrate_full for an 8-param (with-ao4) chi2 fit at this oo_scale
    cf.OBJECTIVE = "chi2"
    cf.OO_SCALE_FIXED = oo_scale
    cf.FIT_REPS = fit_reps
    if not any(b[0] == "k_ao4" for b in cf.BOUNDS):
        cf.WITH_AO4 = True
        cf.BOUNDS.append(cf.AO4_BOUND)
    cf._eval = {"n": 0, "best": np.inf}

    x0 = np.log(BASE_X0)
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy(); pt[i] += 0.4
        simplex.append(pt)
    res = minimize(cf.loss, x0, method="Nelder-Mead",
                   options={"maxiter": 400, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    best = cf.unpack(res.x)
    chi2, sse = cf.scores(best, reps=eval_reps)

    # per-linkage model vs reference (at eval reps) + 4-O-5 residual in sigma
    rows_z = {}
    inv = {}
    for sg in cf.ANCHORS:
        m, _ = cf.model_full(sg, best, reps=eval_reps)
        inv[sg] = m
        rows_z[sg] = {k: (m[k] - cf.REFERENCE[sg][k]) / cf.sigma(cf.REFERENCE[sg][k])
                      for k in cf.ALL_KEYS}
    nfr = {sg: n_frags_at(sg, best, oo_scale, reps=eval_reps) for sg in cf.ANCHORS}
    return {"oo_scale": oo_scale, "best": best, "chi2": chi2, "sse": sse,
            "inv": inv, "z": rows_z, "nfr": nfr}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-reps", type=int, default=8)
    ap.add_argument("--eval-reps", type=int, default=30)
    args = ap.parse_args()

    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("alpha-O-4 in a multi-fragment regime (reviewer 2.4)")
    emit("refit [7 constants + k_ao4] at each fixed oo_scale; fit reps %d, eval reps %d"
         % (args.fit_reps, args.eval_reps))
    emit("baseline = the adopted single-molecule fit: chi2 25.44, k_ao4 rejected"
         " (0.0149)")
    emit("")

    results = []
    for oo in OO_REGIMES:
        emit("--- refitting at oo_scale = %.2f ---" % oo)
        r = refit_at(oo, args.fit_reps, args.eval_reps)
        results.append(r)
        emit("  done: chi2 = %.2f  k_ao4 = %.4f" % (r["chi2"], r["best"]["k_ao4"]))
    emit("")

    LAB = {0.0: "softwood", 1.0: "poplar", 3.0: "high-S"}
    emit("=== summary: does alpha-O-4 earn its place as fragments appear? ===")
    emit("%9s %8s %9s | %-22s | %s"
         % ("oo_scale", "chi2", "k_ao4", "n_frags (sw/pop/hiS)",
            "5-5 & 4-O-5 model|ref (sw/pop/hiS)"))
    for r in results:
        nfr = "/".join("%.1f" % r["nfr"][sg] for sg in cf.ANCHORS)
        c55 = "/".join("%.1f" % (100 * r["inv"][sg]["55"]) for sg in cf.ANCHORS)
        c4o5 = "/".join("%.1f" % (100 * r["inv"][sg]["4o5"]) for sg in cf.ANCHORS)
        ref55 = "/".join("%.1f" % (100 * cf.REFERENCE[sg]["55"]) for sg in cf.ANCHORS)
        ref4o5 = "/".join("%.1f" % (100 * cf.REFERENCE[sg]["4o5"]) for sg in cf.ANCHORS)
        emit("%9.2f %8.2f %9.4f | %-22s | 5-5 %s|%s  4O5 %s|%s"
             % (r["oo_scale"], r["chi2"], r["best"]["k_ao4"], nfr,
                c55, ref55, c4o5, ref4o5))
    emit("")
    emit("4-O-5 residual (z, model-ref in sigma) per regime -- the over-prediction "
         "alpha-O-4 was meant to relieve:")
    emit("%9s %10s %10s %10s" % ("oo_scale", "sw 4O5 z", "pop 4O5 z", "hiS 4O5 z"))
    for r in results:
        emit("%9.2f %10.2f %10.2f %10.2f"
             % (r["oo_scale"], r["z"][0.0]["4o5"], r["z"][1.0]["4o5"],
                r["z"][3.0]["4o5"]))
    emit("")
    emit("=== reading ===")
    emit("If chi2 never beats the single-molecule baseline and k_ao4 is never adopted")
    emit("with a chi2 gain, alpha-O-4's rejection is robust to topology: its only")
    emit("action -- joining fragments -- is invisible to the six-linkage inventory, so")
    emit("no inventory fit adopts it, multi-fragment or not. The oo_scale = 0 regime")
    emit("reaches real fragments but loses 5-5 and 4-O-5 (which alpha-O-4 cannot make),")
    emit("so it cannot be the better model the reviewer hypothesized.")

    os.makedirs(RESULTS, exist_ok=True)
    csv_path = os.path.join(RESULTS, "ao4_multifragment_experiment.csv")
    with open(csv_path, "w") as fh:
        cols = (["oo_scale", "chi2", "sse", "k_ao4"]
                + ["nfrags_%s" % LAB[sg] for sg in cf.ANCHORS]
                + ["m55_%s" % LAB[sg] for sg in cf.ANCHORS]
                + ["m4o5_%s" % LAB[sg] for sg in cf.ANCHORS]
                + ["z4o5_%s" % LAB[sg] for sg in cf.ANCHORS])
        fh.write(",".join(cols) + "\n")
        for r in results:
            vals = ([r["oo_scale"], r["chi2"], r["sse"], r["best"]["k_ao4"]]
                    + [r["nfr"][sg] for sg in cf.ANCHORS]
                    + [100 * r["inv"][sg]["55"] for sg in cf.ANCHORS]
                    + [100 * r["inv"][sg]["4o5"] for sg in cf.ANCHORS]
                    + [r["z"][sg]["4o5"] for sg in cf.ANCHORS])
            fh.write(",".join("%.4f" % v for v in vals) + "\n")
    log_path = os.path.join(RESULTS, "ao4_multifragment_experiment.log")
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    emit("")
    emit("Wrote %s" % csv_path)
    emit("Wrote %s" % log_path)


if __name__ == "__main__":
    main()
