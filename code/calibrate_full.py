#!/usr/bin/env python
"""
calibrate_full.py

Full-inventory calibration of the fresh forward generator. Unlike calibrate.py
(which fit only the beta-side-chain distribution and left the condensed channel
unconstrained -> over-condensed softwood, branch 0.31), this fits the ENTIRE
modeled linkage inventory -- beta-O-4, beta-5, beta-beta, beta-1, 5-5, 4-O-5 --
as fractions of TOTAL linkages (the standard NMR reporting convention), and fits
the condensed rate constants k_55, k_4o5 together with k_b5, k_b1, k_bb and
pool_target (k_bo4 = 8 fixed as scale).

One parameter set must fit all three S/G anchors, so any S/G trend comes from the
S-blocks-C5 chemistry, not per-point fitting. Seeds fixed -> deterministic loss.

Run: conda run -n lignin python code/calibrate_full.py
Outputs: results/calibration_full.log, results/mini_forward_full_calibrated.csv
"""
import os
import importlib.util
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "mfg", os.path.join(HERE, "mini_forward_generator.py"))
mfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mfg)

ALL_KEYS = ["bo4", "b5", "bb", "b1", "55", "4o5"]

# PRIMARY OBJECTIVE, changed sse -> chi2 on 2026-07-16. The primary writes the
# unsuffixed result files; the other writes _<name>.
#
# Why chi2 is primary, having been the sensitivity check:
#   1. sse is the wrong likelihood here and always was. Absolute squared error on
#      the linkage fractions implicitly asserts that beta-O-4 at ~60-89% and 4-O-5
#      at ~0.6-1.5% are measured to the SAME absolute precision. They are not, so
#      sse lets beta-O-4 dominate the fit by sheer magnitude. chi2 weights each
#      linkage by its estimated measurement uncertainty, which is the assumption we
#      can actually defend.
#   2. At the corrected softwood reference (5-5 = 4.0), chi2 wins on BOTH metrics:
#      chi2 25.44 vs 49.34, and sse 0.00574 vs 0.00663 -- i.e. the chi2 optimum is
#      better BY SSE than the sse optimum, at eval reps. That means the sse fit was
#      not converged to its own global optimum; it had overfit the 8-rep fitting
#      noise. A sensitivity check beating the primary on the primary's own metric is
#      not a sensitivity check.
#   3. It fixes beta-beta. Under sse, beta-beta was the worst-fit linkage at every
#      anchor (z = -1.99/-3.14/-2.58); under chi2 it is within ~1.6 sigma
#      (-0.29/-1.57/-1.01). The "model cannot make low 5-5 and high beta-beta at
#      once" conclusion was an sse artifact.
#   4. Held-out error is lower: LOOCV mean 3.4 vs 3.8 pts.
#
# NOTE, and it is the important one: this does NOT reinstate the withdrawn
# k_bo4_nuc claim. chi2 fits k_bo4_nuc = 0.9213, which sits almost exactly on the
# Gani-implied ~0.9 -- and that is a THIRD coincidence, not corroboration. The same
# data and the same reference give 2.2031 under sse. A parameter whose value moves
# 2.4x when only the objective changes is not determined by the data. See
# ../Lignin_BioChem docs/01 Exp 12: sweeping k_bo4_nuc over its whole prior range
# (2000x) moves the modelled inventory by at most 2.5 sigma of simulator noise.
PRIMARY_OBJECTIVE = "chi2"


def _norm(d):
    t = sum(d.values())
    return {k: v / t for k, v in d.items()}


# Reference linkage inventories, raw per-100-C9, normalised to % of the SIX modelled
# linkages. See docs/08_reference_provenance.md for the sourcing of every number.
#
# REVISED 2026-07-16: softwood 5-5 12.0 -> 4.0 and 4-O-5 4.0 -> 1.0.
#
# The old softwood values had NO provenance -- they entered in the generator's first
# commit annotated "from literature (see docs/02)", and docs/02 contains no 5-5
# number at all. They are degradative-method values, and the degradative methods are
# specifically unreliable for condensed linkages:
#   * Thioacidolysis OVER-REPRESENTS condensed dimers. Yue, Lu, Regner, Sun & Ralph
#     (2017) show 4-O-5 at 6.8-9.1% of identified dimers against "no more than
#     approximately 1% of lignin units" -- an ~8x inflation -- and note the products
#     are "prevalent in the thioacidolysis dimers beyond their actual proportion".
#   * Permanganate oxidation is "mostly qualitative, since the low yields of
#     degradation acids make it difficult to estimate the quantitative distribution
#     of structural units" (Sette, Wechselberger & Crestini 2011).
#   * The 5-5 dimer is 19.9-22.5% of identified thioacidolysis dimers, but "no
#     estimate of the conversion yield from the eight-membered ring dibenzodioxocins
#     is yet available" (Yue 2017) -- so that share was never convertible to a
#     per-100-C9 figure in the first place.
# Chang et al. (2020) defend 20-28% of monolignols in biphenyl (=10-14 bonds/100 C9)
# on UV + permanganate + 13C + thioacidolysis/31P; of those, permanganate is
# qualitative, thioacidolysis is inflated, and Capanema (2004) concedes 13C
# quantification of condensed moieties "was less precise than for side chain
# structures".
#
# The replacements are specific NMR measurements:
#   5-5  = 4.0/100 C9 -- Argyropoulos et al. (2002) measure dibenzodioxocin at
#          3.7/100 C9 (31P after alkaline degradation) and 4.3/100 C9 (DFRC), two
#          independent specific methods; Saito (2021) and Yue (2017) both report
#          dibenzodioxocin as "the main biphenyl structure present" in softwood.
#   4-O-5 = 1.0/100 C9 -- Yue et al. (2016, 2017): authenticated by NMR against
#          synthesised models at "no more than approximately 1% of lignin units".
#
# CAVEAT, and it is not small: "the main biphenyl" is not "all biphenyl". If
# dibenzodioxocin is ~65% of softwood biphenyl, the true 5-5 is ~6/100 C9, which is
# also what Chang's softwood:hardwood ratio (~2.4x) implies given HSQC-anchored
# hardwood values. The defensible range is 4-6. We adopt 4.0 as the only directly
# measured figure and report 4-6 sensitivity. What is excluded is 12: it violates
# the beta-beta conservation bound (B/C = 0.16 against a measured p_bb of 0.43-0.57;
# see ../Lignin_BioChem), implies a 4.0x softwood:poplar contrast against Chang's own
# 2.4x, and rests on the methods above.
#
# POPLAR AND HIGH-S ARE UNCHANGED. Their 5-5 (3.8%, 1.1% of linkages) already agree
# with hardwood HSQC dibenzodioxocin (Rencoret 2009 Paulownia 2%; Rosado 2023 orange
# 3%), so the defect is specific to the guaiacyl case -- which is where the
# degradative methods were applied and where condensed linkages are abundant enough
# to inflate. Their provenance is nonetheless undocumented and should be sourced.
#
# NOTE the inventory is INCOMPLETE, for every anchor. Summed, softwood is 66.5
# bonds/100 C9 => 33.5 fragments/100 units => DP ~3, against a reported MWL DP ~20
# (~95 bonds/100 C9 needed). ~30 linkages/100 C9 are missing from the six-linkage
# scheme entirely: noncyclic alpha-O-4 (~4/100), beta-6/alpha-6 (~3/100), alkyl-O-alkyl
# at alpha/gamma (Balakshin 2020, Zhu 2022), and LCC bonds. Normalising to the six
# modelled linkages is self-consistent, but the model's TOPOLOGY cannot be right if
# the missing linkages join chains -- and alpha-O-4 does. See k_ao4.
# Raw per-100-C9 inventories, before normalising to the six modelled linkages. Kept
# separate from REFERENCE so a single raw number can be overridden (see --ref-4o5,
# reviewer 3.2) and the whole reference rebuilt consistently.
RAW_REFERENCE = {
    0.0: {"bo4": 46.0, "b5": 11.0, "bb": 2.5, "b1": 2.0, "55": 4.0, "4o5": 1.0},
    1.0: {"bo4": 62.0, "b5": 6.0,  "bb": 4.0, "b1": 1.5, "55": 3.0,  "4o5": 1.5},
    3.0: {"bo4": 78.0, "b5": 3.0,  "bb": 4.0, "b1": 1.0, "55": 1.0,  "4o5": 0.5},
}


def build_reference(raw=None):
    raw = raw if raw is not None else RAW_REFERENCE
    return {sg: _norm(d) for sg, d in raw.items()}


REFERENCE = build_reference()
ANCHORS = [0.0, 1.0, 3.0]
FIT_REPS = 8
N = 120
K_BO4 = 8.0

# ---- per-linkage uncertainties: the objective is a chi-square ----------------
# The objective used to be a plain sum of squared errors on the fractions, which
# silently asserts that EVERY linkage is measured to the same ABSOLUTE precision
# -- that a 1-pt error on 4-O-5 (0.6% of linkages at high-S) matters exactly as
# much as a 1-pt error on beta-O-4 (89%). It does not, and the consequence was
# real: beta-O-4's errors are 5-7 pts and beta-beta's 1-2, so beta-O-4 dominated
# the loss ~20x and the optimizer simply traded beta-beta away (docs/06 Exp 6).
#
# We weight instead by an estimated measurement sigma per linkage:
#
#   sigma_k = max(REL_SIGMA * ref_k, FLOOR_SIGMA)
#
#   REL_SIGMA = 0.05 -- Bourmaud/Luterbacher 2024 whole-cell-wall HSQC0 report
#     beta-O-4 as 48.0 +/- 2.4, 60.5 +/- 3.0, 72.5 +/- 3.6, 78.9 +/- 3.9 mol%
#     (docs/02 section B): ~5% relative, consistently. We assume the other
#     linkages are integrated to comparable RELATIVE precision.
#   FLOOR_SIGMA = 0.01 -- a floor of 1 percentage point. Without it a purely
#     relative sigma would hand the tiny channels enormous leverage (4-O-5 at
#     0.6% would need to be fit to +/-0.03 pts), i.e. we would over-fit the
#     smallest and least certain numbers in the reference. No HSQC integration of
#     a minor cross-peak is better than ~1 pt.
#
# Effect: beta-O-4 gets sigma ~4.5 pts, the minor channels ~1 pt. So the fit now
# cares about each linkage in proportion to how well it is actually known, rather
# than in proportion to how big it is.
#
# NOTE: the loss VALUE is a chi-square and does NOT compare to the earlier
# sum-of-squares numbers (0.00916, 0.00995). We report both for every fit.
REL_SIGMA = 0.05
FLOOR_SIGMA = 0.01


def sigma(ref_frac):
    return max(REL_SIGMA * ref_frac, FLOOR_SIGMA)

# fit order: k_b5, k_b1, k_bb, k_55, k_4o5, n_nuclei, k_bo4_nuc (pool_target
# fixed at generator default 1.0 -> endwise/linear).
#
# k_bo4_nuc is the beta-O-4 constant at the MONOMER-MONOMER (nucleation) step,
# fitted separately from the chain-extension K_BO4 = 8 scale. It has to be free:
# nucleation is a different reaction, and Gani (docs/04) puts mon-mon beta-O-4 at
# 6.3 kcal/mol against beta-beta's 5.2 -- so beta-O-4 does not dominate there the
# way it dominates chain growth. Reusing K_BO4 = 8 at nucleation forced
# P(beta-beta | G-G) to 0.32 when the DFT barriers imply ~0.57, and that was the
# beta-beta deficit. Fitting it replaces `bb_s_boost` (a flat x4 fudge on S-S
# beta-beta), which is gone: the S/G trend now comes from Gani's pair
# specificity via mini_forward_generator._pair_mod.
# WITH_AO4 (set by --with-ao4) appends k_ao4 as an EIGHTH free parameter. alpha-O-4
# is the non-radical ether route joining two chains (a quinone methide capturing a
# phenol on another fragment); see mini_forward_generator k_ao4 and a5a7170. It is a
# candidate fix for the 4-O-5 over-prediction, which is this fit's dominant residual.
# ⚠ It is an 8th parameter on 18 numbers in a model where 2 of 7 are ALREADY
# unidentifiable (../Lignin_BioChem docs/01 Exp 12). Adding it must be judged on
# whether it is itself identified, not on whether it improves the fit -- an
# unidentified parameter improves any fit.
WITH_AO4 = False
AO4_BOUND = ("k_ao4", 0.005, 40.0)
AO4_X0 = 1.0

# WITH_KNUC (set by --with-knuc) appends k_nuc as an EIGHTH free parameter. k_nuc is
# CONTINUOUS nucleation during growth (two fresh monomers introduced and coupled by
# the monomer-monomer constants), the mechanism reviewer 2.3 asks for to decouple
# beta-beta from the fixed nucleation burst. The sweep (code/knuc_nuclei_experiment.py)
# showed it CAN lift beta-beta past the 4.20% burst ceiling to every reference, but
# that a single GLOBAL k_nuc at the calibrated params does not lower summed chi2
# (minimised at k_nuc = 0). This refit is the decisive test: with the other six
# constants free to absorb the monomer-monomer beta-5/4-O-5/5-5 that continuous
# nucleation also emits, does k_nuc earn its place? Same standard as k_ao4 -- judged
# on whether it is identified and improves the HELD-OUT error, not the training loss.
WITH_KNUC = False
KNUC_BOUND = ("k_nuc", 0.001, 40.0)

# WITH_OLIGOBB (set by --with-oligobb) appends oligo_bb as an EIGHTH free parameter:
# beta-beta via RE-RADICALISED oligomer chain ends (monomer-oligomer and oligomer-
# oligomer resinols), the mechanism a referee says the fixed-nucleation model wrongly
# forbids, making the beta-beta ceiling (allegedly) a self-inflicted artifact. The
# channel demonstrably lifts beta-beta far past the burst ceiling (oligo_bb = 1 gives
# poplar beta-beta ~59%). The decisive test, same standard as k_ao4 / k_nuc: does the
# refit ACCEPT it (lower summed chi2 AND held-out error, fixing the poplar/high-S
# beta-beta deficit) or collapse it toward the floor? Collapse => the ceiling is
# chemistry, not an artifact.
WITH_OLIGOBB = False
OLIGOBB_BOUND = ("oligo_bb", 0.0001, 5.0)

# OO_SCALE_FIXED (set by --oo-scale, or by ao4_multifragment_experiment.py) forces
# oo_scale -- the multiplier on the 5-5 / 4-O-5 (welding) channel -- to a fixed value
# for every model run in the fit, instead of leaving it at the generator default 1.0.
# It is NOT a fitted parameter; it is held so that a fit can be run in a chosen
# connectivity regime. oo_scale = 1 is the calibrated single-molecule regime;
# oo_scale = 0 is the multi-fragment regime (n_frags = n_nuclei). See the connectivity
# experiment (bistable welding) and reviewer concern 2.4.
OO_SCALE_FIXED = None

# REF_4O5_TAG (set by --ref-4o5) is only a filename tag so a sensitivity run does not
# overwrite the adopted result files; the actual override happens on REFERENCE.
REF_4O5_TAG = None

BOUNDS = [("k_b5", 0.05, 40.0), ("k_b1", 0.01, 20.0), ("k_bb", 0.05, 60.0),
          ("k_55", 0.005, 40.0), ("k_4o5", 0.005, 40.0), ("n_nuclei", 1.0, 30.0),
          ("k_bo4_nuc", 0.02, 40.0)]


def full_fracs(bonds):
    t = sum(bonds[k] for k in ALL_KEYS)
    return {k: (bonds[k] / t if t else 0.0) for k in ALL_KEYS}


# SEED_BASE offsets the per-rep RNG seeds. Default 7000 = the calibrated/adopted seed
# set. Varying it (see code/param_recovery_experiment.py) fits a DIFFERENT noise
# realisation, which is how the "fixed-seed fit chases one noise realisation" critique
# is tested: if the fitted optimum barely moves across seed bases, it is not overfit.
SEED_BASE = 7000


def model_full(sg, params, reps=FIT_REPS):
    if OO_SCALE_FIXED is not None:
        params = dict(params)
        params["oo_scale"] = OO_SCALE_FIXED
    agg = {k: 0 for k in ALL_KEYS}
    branch = []
    for r in range(reps):
        g = mfg.MiniForwardKMC(sg, n_target=N, params=params, seed=SEED_BASE + r + int(sg * 100))
        s = g.run()
        for k in ALL_KEYS:
            agg[k] += g.bonds[k]
        branch.append(s["branch_coeff"])
    return full_fracs(agg), float(np.mean(branch))


def unpack(x):
    vals = np.exp(x)
    params = {"k_bo4": K_BO4}
    for (name, lo, hi), v in zip(BOUNDS, vals):
        params[name] = float(np.clip(v, lo, hi))
    return params


_eval = {"n": 0, "best": np.inf}


def scores(params, reps=FIT_REPS):
    """(chi2, sse) for one parameter set.

    chi2 is what we minimise; sse is the old absolute sum-of-squares, carried
    along only so that fits under the two objectives can be compared.
    """
    chi2, sse = 0.0, 0.0
    for sg in ANCHORS:
        m, _ = model_full(sg, params, reps=reps)
        for k in ALL_KEYS:
            r = REFERENCE[sg][k]
            d = m[k] - r
            chi2 += (d / sigma(r)) ** 2
            sse += d ** 2
    return chi2, sse


OBJECTIVE = "sse"     # set by --objective; see main()


def loss(x):
    chi2, sse = scores(unpack(x))
    L = chi2 if OBJECTIVE == "chi2" else sse
    _eval["n"] += 1
    if L < _eval["best"]:
        _eval["best"] = L
    if _eval["n"] % 20 == 0:
        print("  eval %3d  %s=%.5g  best=%.5g" % (_eval["n"], OBJECTIVE, L, _eval["best"]))
    return L


def main(objective="sse"):
    global OBJECTIVE
    OBJECTIVE = objective
    out_dir = os.path.join(HERE, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)
    suffix = "" if objective == PRIMARY_OBJECTIVE else "_" + objective
    if WITH_AO4:
        suffix += "_ao4"
    if WITH_KNUC:
        suffix += "_knuc"
    if WITH_OLIGOBB:
        suffix += "_oligobb"
    if REF_4O5_TAG is not None:
        suffix += "_ref4o5-%s" % REF_4O5_TAG

    # start from the previous fit, plus k_bo4_nuc seeded at the Gani-implied
    # monomer-monomer ratio: beta-O-4 6.3 vs beta-beta 5.2 kcal/mol -> beta-O-4
    # is ~6x SLOWER than beta-beta at this step, i.e. k_bo4_nuc ~ k_bb/6, nowhere
    # near the chain-extension K_BO4 = 8 it used to inherit.
    #                k_b5  k_b1  k_bb  k_55  k_4o5 n_nuclei k_bo4_nuc
    base_x0 = [2.90, 0.26, 5.35, 2.66, 0.21, 6.0, 0.9]
    if WITH_AO4:
        # Start the 8-param fit AT the reported 7-param optimum with k_ao4 nearly
        # off. The 8-param model NESTS the 7-param one (k_ao4 -> its lower bound), so
        # from here the fit can only improve or stay -- any worse result is an
        # optimizer failure, not evidence about alpha-O-4. Starting instead from the
        # generic x0 gave chi2 36.75 against the 7-param 25.44, which is impossible
        # for a nested model and diagnosed exactly that: Nelder-Mead does not
        # navigate 8 dims from a cold start here.
        base_x0 = [2.4904, 0.2715, 16.7825, 0.7108, 0.1737, 4.716, 0.9213, 0.01]
    if WITH_KNUC:
        # Same nesting logic as WITH_AO4: start AT the 7-param chi2 optimum with
        # k_nuc at its lower bound, so the 8-param model can only improve. A worse
        # result would be an optimizer failure, not evidence about k_nuc.
        base_x0 = [2.4904, 0.2715, 16.7825, 0.7108, 0.1737, 4.716, 0.9213, 0.001]
    if WITH_OLIGOBB:
        # Same nesting logic: start AT the 7-param chi2 optimum with oligo_bb at its
        # floor (channel off), so the 8-param model NESTS the 7-param one and can only
        # improve. Any worse result is optimizer failure, not evidence about oligo_bb.
        base_x0 = [2.4904, 0.2715, 16.7825, 0.7108, 0.1737, 4.716, 0.9213, 0.0001]
    x0 = np.log(base_x0)
    print("Full-inventory calibration (fraction of total linkages)")
    print("anchors:", ANCHORS, " fit reps:", FIT_REPS, " k_bo4 fixed:", K_BO4)
    print("reference inventory (%% of total linkages):")
    for sg in ANCHORS:
        print("   S/G %.1f: %s" % (sg, {k: round(100 * v, 1) for k, v in REFERENCE[sg].items()}))
    print()

    # wider initial simplex (~x1.5 per dim in log space) to force exploration
    simplex = [x0]
    for i in range(len(x0)):
        pt = x0.copy(); pt[i] += 0.4
        simplex.append(pt)
    res = minimize(loss, x0, method="Nelder-Mead",
                   # fatol must follow the objective's SCALE: sse lives at ~0.01,
                   # chi2 at ~10-100 over 18 terms. One value cannot serve both --
                   # 1e-6 on chi2 just runs to maxiter, 1e-3 on sse stops early.
                   options={"maxiter": 400, "xatol": 1e-2,
                            "fatol": 1e-3 if objective == "chi2" else 1e-6,
                            "initial_simplex": np.array(simplex)})
    best = unpack(res.x)
    chi2, sse = scores(best, reps=30)
    print("\nfitted params:", {k: round(v, 4) for k, v in best.items()})
    print("minimised: %s = %.5g  (evals: %d)" % (objective, res.fun, _eval["n"]))
    print("this fit scores:  sse = %.5f   chi2 = %.2f" % (sse, chi2))
    print("  Both are printed for every fit because they are different objectives")
    print("  and their VALUES do not compare. Judge a fit by the per-linkage")
    print("  residuals below and by validate.py's held-out error, not by a loss.\n")

    print("full-inventory fit vs reference (model | ref, %% of total linkages):")
    print("%5s  %-11s %-11s %-11s %-11s %-11s %-11s  %s" %
          ("S/G", "bO4", "b5", "bb", "b1", "5-5", "4-O-5", "branch"))
    for sg in ANCHORS:
        m, br = model_full(sg, best, reps=30)
        cells = ["%4.1f|%4.1f" % (100 * m[k], 100 * REFERENCE[sg][k]) for k in ALL_KEYS]
        print("%5.1f  %-11s %-11s %-11s %-11s %-11s %-11s  %.3f" % (sg, *cells, br))

    # Residuals in units of sigma -- this is what the fit now actually sees, and
    # it is the honest way to read which linkages are and are not fit.
    print("\nresiduals in sigma ((model - ref) / sigma; |z| < 1 is within measurement):")
    print("%5s  %7s %7s %7s %7s %7s %7s" % ("S/G", "bO4", "b5", "bb", "b1", "5-5", "4-O-5"))
    for sg in ANCHORS:
        m, _ = model_full(sg, best, reps=30)
        zs = [(m[k] - REFERENCE[sg][k]) / sigma(REFERENCE[sg][k]) for k in ALL_KEYS]
        print("%5.1f  %s" % (sg, " ".join("%+7.2f" % z for z in zs)))

    # full S/G sweep with fitted params
    print("\nfull S/G sweep with fitted params (%% of total linkages):")
    print("%5s %7s %6s %6s %6s %6s %6s %7s" %
          ("S/G", "bO4", "b5", "bb", "b1", "5-5", "4o5", "branch"))
    sweep = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
    csv = os.path.join(out_dir, "mini_forward_full_calibrated%s.csv" % suffix)
    with open(csv, "w") as fh:
        fh.write("sg_ratio,bo4,b5,bb,b1,55,4o5,branch\n")
        for sg in sweep:
            m, br = model_full(sg, best, reps=30)
            print("%5.2f %6.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %7.3f" %
                  (sg, 100 * m["bo4"], 100 * m["b5"], 100 * m["bb"], 100 * m["b1"],
                   100 * m["55"], 100 * m["4o5"], br))
            fh.write("%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %
                     (sg, m["bo4"], m["b5"], m["bb"], m["b1"], m["55"], m["4o5"], br))
    print("\nfitted params:", {k: round(v, 4) for k, v in best.items()})
    print("Wrote %s" % csv)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--objective", choices=["sse", "chi2"], default=PRIMARY_OBJECTIVE,
                    help="'chi2' (default since 2026-07-16, and the reported fit): "
                         "weight each linkage by its estimated measurement sigma "
                         "(see REL_SIGMA/FLOOR_SIGMA). 'sse': absolute squared error "
                         "on the linkage fractions -- implicitly assumes every "
                         "linkage is measured to the same absolute precision, so "
                         "beta-O-4 dominates by magnitude alone. sse is retained as "
                         "the sensitivity check. See PRIMARY_OBJECTIVE for why they "
                         "swapped: at the corrected softwood reference chi2 beats "
                         "sse on BOTH metrics, including sse's own.")
    ap.add_argument("--with-ao4", action="store_true",
                    help="fit k_ao4 as an 8th parameter: the alpha-O-4 ether route "
                         "that joins two chains without spending a condensed bond. "
                         "Candidate fix for the 4-O-5 residual.")
    ap.add_argument("--with-knuc", action="store_true",
                    help="fit k_nuc as an 8th parameter: continuous nucleation during "
                         "growth, the mechanism that decouples beta-beta from the "
                         "fixed nucleation burst (reviewer 2.3). Candidate fix for the "
                         "beta-beta ceiling.")
    ap.add_argument("--with-oligobb", action="store_true",
                    help="fit oligo_bb as an 8th parameter: beta-beta via re-radicalised "
                         "OLIGOMER chain ends (monomer-oligomer and oligomer-oligomer "
                         "resinols). The referee mechanism that would lift the beta-beta "
                         "ceiling if the fit accepts it. Judged on whether it lowers "
                         "held-out error or collapses to the floor.")
    ap.add_argument("--oo-scale", type=float, default=None,
                    help="hold oo_scale (the 5-5/4-O-5 welding-channel multiplier) "
                         "fixed at this value for every model run, instead of the "
                         "generator default 1.0. Use with --with-ao4 to fit alpha-O-4 "
                         "in a chosen connectivity regime (reviewer 2.4); oo_scale=0 "
                         "is multi-fragment, oo_scale=1 the calibrated single molecule.")
    ap.add_argument("--ref-4o5", type=float, default=None,
                    help="override the SOFTWOOD (S/G 0) 4-O-5 reference (raw per-100-C9, "
                         "default 1.0) and rebuild the normalised reference. The adopted "
                         "1.0 is the only directly NMR-measured figure; thioacidolysis "
                         "inflates it ~8x. Reviewer 3.2 asks what the fitted constants do "
                         "if the reference is higher -- run 2.0 and 3.0. Writes a _ref4o5-<v> "
                         "suffix so it never overwrites the adopted result files.")
    _a = ap.parse_args()
    if _a.ref_4o5 is not None:
        raw = {sg: dict(d) for sg, d in RAW_REFERENCE.items()}
        raw[0.0]["4o5"] = _a.ref_4o5
        globals()["REFERENCE"] = build_reference(raw)
        globals()["REF_4O5_TAG"] = ("%g" % _a.ref_4o5).replace(".", "p")
    if _a.with_ao4:
        globals()["WITH_AO4"] = True
        BOUNDS.append(AO4_BOUND)
    if _a.with_knuc:
        globals()["WITH_KNUC"] = True
        BOUNDS.append(KNUC_BOUND)
    if _a.with_oligobb:
        globals()["WITH_OLIGOBB"] = True
        BOUNDS.append(OLIGOBB_BOUND)
    if _a.oo_scale is not None:
        globals()["OO_SCALE_FIXED"] = _a.oo_scale
    main(objective=_a.objective)
