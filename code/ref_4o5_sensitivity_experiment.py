#!/usr/bin/env python
"""
ref_4o5_sensitivity_experiment.py  --  reviewer concern 3.2

The adopted softwood 4-O-5 reference is 1.0 / 100 C9 -- the only directly NMR-measured
figure (Yue 2016/2017), against thioacidolysis values that inflate condensed dimers
~8x. The fit still leaves a systematic 4-O-5 over-prediction (~ +2 sigma). Reviewer 3.2
asks the fair question: is that residual real, or is the model just fitting a reference
that is itself too low? If the reference were 2.0 or 3.0, would the fitted constants
move a lot (=> the reference drives the fit) or barely (=> 4-O-5 is nearly unidentified
and the residual is not a reference artefact)?

This RE-RUNS THE WHOLE CHI2 FIT at softwood 4-O-5 in {1.0 (adopted), 2.0, 3.0}, using
calibrate_full's own reference builder, loss and optimiser (same x0, simplex, tol), so
the only thing that changes between rows is the one reference number. It tabulates, per
row: the 7 fitted constants, chi2/sse, and the 4-O-5 model value and residual z at
every anchor. Nothing is adopted; the adopted fit (fitted_params.json) is unchanged.

Run:     conda run -n lignin python code/ref_4o5_sensitivity_experiment.py
Outputs: results/ref_4o5_sensitivity.csv, results/ref_4o5_sensitivity.log
"""
import os
import importlib.util
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "cf", os.path.join(HERE, "calibrate_full.py"))
cf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cf)

REF_VALUES = [1.0, 2.0, 3.0]     # softwood 4-O-5 per-100-C9 (1.0 = adopted)
EVAL_REPS = 30
PARAM_ORDER = ["k_bo4", "k_bo4_nuc", "k_b5", "k_b1", "k_bb", "k_55", "k_4o5", "n_nuclei"]


def fit_at(ref_4o5):
    """Rebuild the reference with this softwood 4-O-5 and re-run the chi2 fit.

    Mirrors calibrate_full.main() exactly: same start point, simplex, tolerances.
    Returns (params, chi2, sse).
    """
    raw = {sg: dict(d) for sg, d in cf.RAW_REFERENCE.items()}
    raw[0.0]["4o5"] = ref_4o5
    cf.REFERENCE = cf.build_reference(raw)      # loss/scores/model_full read this global
    cf.OBJECTIVE = "chi2"
    cf._eval = {"n": 0, "best": np.inf}

    base_x0 = [2.90, 0.26, 5.35, 2.66, 0.21, 6.0, 0.9]
    x0 = np.log(base_x0)
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy(); pt[i] += 0.4
        simplex.append(pt)
    res = minimize(cf.loss, x0, method="Nelder-Mead",
                   options={"maxiter": 400, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    best = cf.unpack(res.x)
    chi2, sse = cf.scores(best, reps=EVAL_REPS)
    return best, chi2, sse


def main():
    out_dir = os.path.join(HERE, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("4-O-5 reference sensitivity: re-fit chi2 at softwood 4-O-5 in {1.0, 2.0, 3.0}")
    emit("(1.0 = adopted, the only NMR-measured value; 2.0/3.0 probe the reviewer's worry")
    emit(" that the +2sigma 4-O-5 residual is just an underestimated reference.)")
    emit("Only the softwood 4-O-5 number changes between rows; fit machinery is identical.")
    emit("")

    results = []
    for rv in REF_VALUES:
        best, chi2, sse = fit_at(rv)
        # 4-O-5 model value and residual z at each anchor, under THIS row's reference
        model_4o5, z_4o5 = {}, {}
        for sg in cf.ANCHORS:
            m, _ = cf.model_full(sg, best, reps=EVAL_REPS)
            r = cf.REFERENCE[sg]["4o5"]
            model_4o5[sg] = m["4o5"]
            z_4o5[sg] = (m["4o5"] - r) / cf.sigma(r)
        results.append((rv, best, chi2, sse, model_4o5, z_4o5))

    # ---- fitted constants table --------------------------------------------
    emit("FITTED CONSTANTS vs softwood 4-O-5 reference:")
    emit("%-10s " % "ref 4-O-5" + "".join("%10s" % p.replace("k_", "") for p in PARAM_ORDER)
         + "%9s%9s" % ("chi2", "sse"))
    for rv, best, chi2, sse, _, _ in results:
        row = "%-10.1f " % rv + "".join("%10.4f" % best[p] for p in PARAM_ORDER)
        row += "%9.2f%9.5f" % (chi2, sse)
        emit(row)
    emit("")

    # relative shift of each constant from the adopted (row 0)
    emit("Constant shift relative to the adopted 1.0 fit (x-fold):")
    base = results[0][1]
    emit("%-10s " % "ref 4-O-5" + "".join("%10s" % p.replace("k_", "") for p in PARAM_ORDER))
    for rv, best, _, _, _, _ in results:
        row = "%-10.1f " % rv + "".join("%10.2f" % (best[p] / base[p]) for p in PARAM_ORDER)
        emit(row)
    emit("")

    # ---- 4-O-5 fit quality table -------------------------------------------
    emit("4-O-5 model value (%) and residual z at each anchor, per reference:")
    emit("%-10s  %-16s %-16s %-16s" % ("ref 4-O-5",
         "softwood(S/G0)", "poplar(S/G1)", "high-S(S/G3)"))
    for rv, _, _, _, model_4o5, z_4o5 in results:
        cells = ["%4.1f%% (z=%+.2f)" % (100 * model_4o5[sg], z_4o5[sg]) for sg in cf.ANCHORS]
        emit("%-10.1f  %-16s %-16s %-16s" % (rv, *cells))
    emit("")

    # ---- verdict ------------------------------------------------------------
    k4o5 = [r[1]["k_4o5"] for r in results]
    z0 = [r[5][0.0] for r in results]
    emit("READING:")
    emit("  * k_4o5 moves %.2fx -> %.2fx over the reference range (1->3): %s"
         % (k4o5[1] / k4o5[0], k4o5[2] / k4o5[0],
            "the constant DOES follow the reference" if k4o5[2] / k4o5[0] > 1.5
            else "the constant barely moves -- 4-O-5 is near-unidentified"))
    emit("  * softwood 4-O-5 residual z: %+.2f (ref 1.0) -> %+.2f (ref 2.0) -> %+.2f (ref 3.0)"
         % (z0[0], z0[1], z0[2]))
    emit("    If |z| shrinks as the reference rises, the +sigma residual at 1.0 is")
    emit("    consistent with an underestimated reference (reviewer's hypothesis);")
    emit("    the point is that this costs nothing elsewhere only if the OTHER linkages'")
    emit("    fits are unchanged -- check chi2 and the other constants above.")
    emit("")

    # ---- CSV ----------------------------------------------------------------
    csv_path = os.path.join(out_dir, "ref_4o5_sensitivity.csv")
    with open(csv_path, "w") as fh:
        fh.write("ref_4o5," + ",".join(PARAM_ORDER) + ",chi2,sse,"
                 + "model4o5_sg0,z_sg0,model4o5_sg1,z_sg1,model4o5_sg3,z_sg3\n")
        for rv, best, chi2, sse, model_4o5, z_4o5 in results:
            fh.write("%.1f," % rv + ",".join("%.4f" % best[p] for p in PARAM_ORDER))
            fh.write(",%.2f,%.5f" % (chi2, sse))
            for sg in cf.ANCHORS:
                fh.write(",%.4f,%.3f" % (model_4o5[sg], z_4o5[sg]))
            fh.write("\n")

    log_path = os.path.join(out_dir, "ref_4o5_sensitivity.log")
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    emit("wrote %s" % os.path.relpath(csv_path, os.path.join(HERE, os.pardir)))
    emit("wrote %s" % os.path.relpath(log_path, os.path.join(HERE, os.pardir)))


if __name__ == "__main__":
    main()
