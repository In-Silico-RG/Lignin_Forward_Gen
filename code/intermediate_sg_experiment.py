#!/usr/bin/env python
"""
intermediate_sg_experiment.py  --  reviewer concern 2.2

The model is calibrated on THREE S/G anchors (0, 1, 3). Reviewer 2.2 asks whether
the claim that it "predicts the complete native linkage inventory from monomer
composition (S/G) alone ... across feedstocks" is supported, given n=3. The ask is
specific and modest:
  (a) predictions at >=2-3 INTERMEDIATE S/G not used in calibration (0.5, 1.5, 2.0),
      shown to be smooth and physically reasonable;
  (b) an explicit statement of the S/G validity range -- interpolation only, or
      extrapolation too? what about S/G > 3, tension/compression wood, mutants?

This script answers (a) and supplies the evidence for (b). It runs the ADOPTED fit
(results/fitted_params.json, chi2 = 25.44) over a dense S/G grid with error bars
(std over seeds), and reports three things:

  1. TABLE + CSV of the inventory at every S/G, anchors and intermediates flagged.
  2. SMOOTHNESS: for each of the six linkages, is the curve monotone across the
     interpolation window [0,3]? does every intermediate lie between its bracketing
     neighbours (no oscillation)? Reported per linkage.
  3. PHYSICAL LIMIT: the S-blocks-C5 rule (mini_forward_generator: S is methoxylated
     at C5, so S cannot donate/accept a C5 bond) forces beta-5, 5-5 and 4-O-5 -- the
     three C5-dependent linkages -- toward ZERO as S/G -> inf, while the C5-INDEPENDENT
     linkages (beta-O-4, beta-beta, beta-1) survive. This is a MECHANISM-level
     guarantee, not a fitted trend, and it is the physical-reasonableness test the
     reviewer asks for. We push the sweep to S/G = 4, 6, 10 and a near-limit S/G = 49
     (98% syringyl) and check the C5 channels collapse and track the G fraction
     1/(1+S/G).

No parameters are refit. This is a prediction/validity experiment on the frozen fit.

Run:     conda run -n lignin python code/intermediate_sg_experiment.py
Outputs: results/intermediate_sg_experiment.csv, results/intermediate_sg_experiment.log
"""
import os
import json
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "mfg", os.path.join(HERE, "mini_forward_generator.py"))
mfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mfg)

ALL_KEYS = ["bo4", "b5", "bb", "b1", "55", "4o5"]
LABELS = {"bo4": "beta-O-4", "b5": "beta-5", "bb": "beta-beta",
          "b1": "beta-1", "55": "5-5", "4o5": "4-O-5"}
# C5-dependent linkages: require a free C5, which syringyl does NOT have. These are
# the ones the S-blocks-C5 mechanism must drive to zero as S/G -> inf.
C5_DEPENDENT = ["b5", "55", "4o5"]
C5_FREE = ["bo4", "bb", "b1"]

N = 120           # matches calibrate_full.N (one 120-unit molecule)
REPS = 60         # more than the fit's 30 -- we want error bars, not just a point

ANCHORS = {0.0, 1.0, 3.0}
# interpolation grid (the reviewer's 0.5/1.5/2.0 are here) + extrapolation tail
INTERP = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]
EXTRAP = [4.0, 6.0, 10.0, 49.0]      # 49 -> 98% syringyl, a proxy for S/G -> inf
GRID = INTERP + EXTRAP


def load_params():
    with open(os.path.join(HERE, os.pardir, "results", "fitted_params.json")) as fh:
        raw = json.load(fh)
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def full_fracs(bonds):
    t = sum(bonds[k] for k in ALL_KEYS)
    return {k: (bonds[k] / t if t else 0.0) for k in ALL_KEYS}


def model_full_stats(sg, params, reps=REPS):
    """Mean and std of each linkage fraction over `reps` seeds, plus mean branch."""
    per_rep = {k: [] for k in ALL_KEYS}
    branch = []
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=N, params=params, seed=7000 + r + int(sg * 100))
        s = g.run()
        fr = full_fracs(g.bonds)
        for k in ALL_KEYS:
            per_rep[k].append(fr[k])
        branch.append(s["branch_coeff"])
    mean = {k: float(np.mean(per_rep[k])) for k in ALL_KEYS}
    std = {k: float(np.std(per_rep[k])) for k in ALL_KEYS}
    return mean, std, float(np.mean(branch))


def main():
    params = load_params()
    out_dir = os.path.join(HERE, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "intermediate_sg_experiment.log")
    csv_path = os.path.join(out_dir, "intermediate_sg_experiment.csv")

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("Intermediate-S/G predictions on the adopted fit (chi2 = 25.44)")
    emit("params: " + json.dumps({k: round(v, 4) for k, v in params.items()}))
    emit("N = %d units, reps = %d (error bars = std over seeds)" % (N, REPS))
    emit("anchors (calibration): S/G = 0, 1, 3.  Everything else is a PREDICTION.")
    emit("")

    # ---- compute the whole grid --------------------------------------------
    means, stds, branches = {}, {}, {}
    for sg in GRID:
        m, sd, br = model_full_stats(sg, params)
        means[sg], stds[sg], branches[sg] = m, sd, br

    # ---- table --------------------------------------------------------------
    emit("Inventory over S/G grid (% of the six modelled linkages, mean +/- std):")
    hdr = "%6s %3s  " % ("S/G", "") + "".join("%-14s" % LABELS[k] for k in ALL_KEYS) + "branch  1/(1+S/G)"
    emit(hdr)
    for sg in GRID:
        m, sd = means[sg], stds[sg]
        tag = "A" if sg in ANCHORS else ("X" if sg in EXTRAP else "i")
        cells = "".join("%5.1f+-%3.1f  " % (100 * m[k], 100 * sd[k]) for k in ALL_KEYS)
        emit("%6.2f %3s  %s%6.3f  %7.3f" % (sg, tag, cells, branches[sg], 1.0 / (1.0 + sg)))
    emit("  legend: A = calibration anchor, i = interpolated prediction, X = extrapolation")
    emit("")

    # ---- (2) smoothness / monotonicity over the interpolation window --------
    emit("SMOOTHNESS over the interpolation window S/G in [0,3]:")
    emit("  A linkage has a real S/G TREND if its end-to-end change (S/G 0 -> 3)")
    emit("  exceeds the run-to-run scatter (single-seed std); the trend is SMOOTH if no")
    emit("  intermediate step reverses that direction by more than that scatter.")
    interp_sorted = sorted(INTERP)
    for k in ALL_KEYS:
        seq = np.array([means[sg][k] for sg in interp_sorted])
        diffs = np.diff(seq)
        noise = float(np.mean([stds[sg][k] for sg in interp_sorted]))
        overall = seq[-1] - seq[0]
        if abs(overall) <= noise:
            verdict = "FLAT (end-to-end change %.1f pts < scatter %.1f) -- no S/G trend" % (
                100 * overall, 100 * noise)
        else:
            direction = np.sign(overall)
            reversals = diffs[np.sign(diffs) != direction]
            big_rev = reversals[np.abs(reversals) > noise]
            trend = "UP" if direction > 0 else "DOWN"
            if len(big_rev) == 0:
                verdict = "monotone %s and SMOOTH (%.1f -> %.1f pts, no reversal > scatter)" % (
                    trend, 100 * seq[0], 100 * seq[-1])
            else:
                verdict = "%s overall but has a reversal > scatter -- CHECK (%.1f -> %.1f pts)" % (
                    trend, 100 * seq[0], 100 * seq[-1])
        emit("  %-9s : %s" % (LABELS[k], verdict))
    emit("")

    # ---- (3) physical limit: C5 channels collapse as S/G -> inf -------------
    emit("PHYSICAL LIMIT (S-blocks-C5): the three C5-dependent linkages must vanish")
    emit("as S/G -> inf (syringyl has no free C5); the C5-free linkages must survive.")
    emit("  C5-dependent: beta-5, 5-5, 4-O-5   |   C5-free: beta-O-4, beta-beta, beta-1")
    m0, mInf = means[0.0], means[49.0]
    for k in C5_DEPENDENT:
        emit("  %-9s : S/G 0 = %5.1f%%  ->  S/G 49 = %5.2f%%   (%.0f%% of its softwood value)"
             % (LABELS[k], 100 * m0[k], 100 * mInf[k], 100 * mInf[k] / max(m0[k], 1e-9)))
    emit("  sum of C5-dependent: %5.1f%% -> %5.2f%%" %
         (100 * sum(m0[k] for k in C5_DEPENDENT), 100 * sum(mInf[k] for k in C5_DEPENDENT)))
    emit("  sum of C5-free     : %5.1f%% -> %5.1f%%" %
         (100 * sum(m0[k] for k in C5_FREE), 100 * sum(mInf[k] for k in C5_FREE)))
    emit("")

    # correlation of each C5-dependent channel with the G fraction 1/(1+S/G)
    emit("Does each C5 channel track G availability g = 1/(1+S/G)? (Pearson r over grid)")
    g_frac = np.array([1.0 / (1.0 + sg) for sg in GRID])
    for k in C5_DEPENDENT:
        y = np.array([means[sg][k] for sg in GRID])
        r = float(np.corrcoef(g_frac, y)[0, 1])
        emit("  %-9s vs g : r = %+.3f" % (LABELS[k], r))
    emit("")

    # ---- CSV ----------------------------------------------------------------
    with open(csv_path, "w") as fh:
        fh.write("sg_ratio,role,g_frac," +
                 ",".join("%s,%s_std" % (k, k) for k in ALL_KEYS) + ",branch\n")
        for sg in GRID:
            role = "anchor" if sg in ANCHORS else ("extrap" if sg in EXTRAP else "interp")
            m, sd = means[sg], stds[sg]
            fh.write("%.2f,%s,%.4f," % (sg, role, 1.0 / (1.0 + sg)))
            fh.write(",".join("%.4f,%.4f" % (m[k], sd[k]) for k in ALL_KEYS))
            fh.write(",%.4f\n" % branches[sg])

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    emit("wrote %s" % os.path.relpath(csv_path, os.path.join(HERE, os.pardir)))
    emit("wrote %s" % os.path.relpath(log_path, os.path.join(HERE, os.pardir)))


if __name__ == "__main__":
    main()
