#!/usr/bin/env python
"""
cross_common.py

Shared configuration + statistics for the cross-algorithm validation
(manuscript section 3.1). Imported by:

    cross_lgs_matched.py             (LGS, sequence-controlled enumeration)
    cross_kmc_matched.py             (Lignin-KMC, incremental)
    cross_ligningraphs_incremental.py(LigninGraphs, incremental)
    cross_wallclock_stability.py     (is t_gen a measurement of chemistry?)
    cross_algorithm_analysis.py      (aggregation + claim audit)

RECONSTRUCTION NOTICE
---------------------
The original scripts for this analysis were lost with a deleted session
scratchpad (see docs/03 and docs/05 "Known losses"); only a prose summary of the
results survived. These scripts are a RE-IMPLEMENTATION from that summary plus
the mechanism documented in docs/02 section D and docs/03. They are not a
byte-recovery of the originals, and interpretive choices are flagged inline as
RECONSTRUCTION CHOICE. Numbers they produce must be reported as what these
scripts do, not asserted to equal the archived values.

Note that the archived prose was WRONG about LGS's mechanism -- it is an
enumerator, not a rejection sampler (see cross_lgs_matched.py, read from the LGS
source). A "whole-structure rejection port" built on the archived description
lived here until 2026-07-15 and was removed: it mirrored a mechanism no generator
has, so its outputs described a fiction. It is in git history (886e252, d06b3b9)
if the reasoning is ever needed.

Kept python 3.6 compatible on purpose: the `lignin_kmc` conda env is py3.6, so
no f-strings and no walrus operator anywhere in this module.
"""
import os

import numpy as np

# ---- the matched parameterization (docs/03 "Setup") -------------------------
# Three feedstock configurations run against every generator with identical
# targets, so that any difference between generators is the generator, not the
# parameterization. These targets are the OLD (dead) manuscript's targets; they
# are kept verbatim because the point of section 3.1 is to audit that manuscript's
# claim on its own terms, and they must NOT be "corrected" here. They are, note,
# chemically off in DIFFERENT directions per config, which is itself worth knowing:
# docs/02 section A puts SW beta-O-4 = 0.54 some 4-10 pts ABOVE the literature
# range (44-50%), while HW = 0.69 at S/G 1.5 sits ~13 pts BELOW what that
# chemistry natively makes (~82%).
CONFIGS = [
    # label, S/G ratio, target beta-O-4 fraction
    ("SW", 0.0, 0.54),   # softwood  -- pure guaiacyl
    ("HW", 1.5, 0.69),   # hardwood  -- syringyl-rich
    ("GR", 0.8, 0.77),   # grass
]

DP = 10           # degree of polymerization (monomers per structure)
N_REPS = 30       # replicates per configuration
BASE_SEED = 20260714

LABELS = [c[0] for c in CONFIGS]

# ---- per-archetype target linkage profiles ----------------------------------
# Shared by the two BACKWARD arms (rejection port, LigninGraphs), which must be
# handed the same target to be comparable. Lignin-KMC is forward and takes none.
#
# The old manuscript specified only a beta-O-4 target per config; the remaining
# channels' target weights must come from somewhere.
#
# RECONSTRUCTION CHOICE: anchor on this project's own NATIVE reference
# inventories (docs/06 Exp 4c, % of total linkages) indexed by S/G, and
# interpolate the non-beta-O-4 shares linearly in S/G.
#
#  * It must be PER-ARCHETYPE. Real LGS gives every config its own `linkages:`
#    block -- confirmed in every YAML of data/LGS_runs/POTLS_LGS_configs.zip. An
#    earlier version of these scripts used ONE softwood-derived residual for all
#    three configs, which asked S-rich HW to make ~39% of its non-bo4 mass as
#    5-5. S blocks C5, so that is impossible and HW never accepted (>20 000
#    attempts). A target must at least be MAKEABLE by the chemistry handed it.
#  * Anchoring on measured native lignins gives makeability without us
#    hand-picking numbers per config.
#  * It does NOT bake in the result: these shares only make the non-bo4 channels
#    reachable. Whether softwood is the rejection maximum still depends on where
#    each config's beta-O-4 TARGET falls in the distribution its own chemistry
#    produces (docs/02 section D).
#
# The SW/HW/GR profiles the lost scripts actually used are NOT archived (only the
# oil-palm POTLS configs survived) -- this is a reconstruction, not a recovery.
BO4, B5, B1, BB, C55, C4O5 = "bo4", "b5", "b1", "bb", "55", "4o5"
LINKAGES = [BO4, B5, B1, BB, C55, C4O5]
NON_BO4 = [B5, BB, B1, C55, C4O5]

NATIVE_REF = {
    # S/G: {linkage: % of total linkages}  -- docs/06 Exp 4c reference column
    0.0: {BO4: 59.0, B5: 14.0, BB: 3.0, B1: 3.0, C55: 16.0, C4O5: 5.0},  # softwood
    1.0: {BO4: 80.0, B5:  8.0, BB: 5.0, B1: 2.0, C55:  4.0, C4O5: 2.0},  # poplar
    3.0: {BO4: 89.0, B5:  3.0, BB: 5.0, B1: 1.0, C55:  1.0, C4O5: 1.0},  # high-S
}


def sg_to_pct_s(sg_ratio):
    """Convert an S/G ratio to the syringyl fraction of monomers."""
    return sg_ratio / (1.0 + sg_ratio)


def residual_share(sg_ratio):
    """Non-beta-O-4 channel shares (summing to 1) at this S/G, interpolated
    linearly between the NATIVE_REF anchors and clamped outside their range."""
    anchors = sorted(NATIVE_REF)

    def shares_at(sg):
        ref = NATIVE_REF[sg]
        tot = sum(ref[k] for k in NON_BO4)
        return dict((k, ref[k] / tot) for k in NON_BO4)

    if sg_ratio <= anchors[0]:
        return shares_at(anchors[0])
    if sg_ratio >= anchors[-1]:
        return shares_at(anchors[-1])
    for lo, hi in zip(anchors[:-1], anchors[1:]):
        if lo <= sg_ratio <= hi:
            w = (sg_ratio - lo) / (hi - lo)
            a, b = shares_at(lo), shares_at(hi)
            return dict((k, a[k] + w * (b[k] - a[k])) for k in NON_BO4)


def target_profile(target_bo4, sg_ratio):
    """The config's full target linkage profile: the manuscript's beta-O-4 target
    plus a chemistry-appropriate residual for this S/G."""
    w = {BO4: target_bo4}
    for k, share in residual_share(sg_ratio).items():
        w[k] = (1.0 - target_bo4) * share
    return w


def results_dir():
    """Absolute path of the repo's results/ directory (created if absent)."""
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, os.pardir, "results")
    out = os.path.abspath(out)
    if not os.path.isdir(out):
        os.makedirs(out)
    return out


def write_csv(path, cols, rows):
    """Write a list-of-dicts to CSV with an explicit column order."""
    with open(path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r[c]
                if isinstance(v, float):
                    vals.append("%.6f" % v)
                else:
                    vals.append(str(v))
            fh.write(",".join(vals) + "\n")


# ---- statistics -------------------------------------------------------------

def one_way_anova(groups):
    """One-way ANOVA over a list of 1-D sample arrays.

    Returns (F, p, eta_squared). eta^2 = SS_between / SS_total is the effect
    size reported in docs/03 alongside the p-value.
    """
    from scipy import stats

    groups = [np.asarray(g, dtype=float) for g in groups]
    F, p = stats.f_oneway(*groups)

    grand = np.concatenate(groups)
    grand_mean = grand.mean()
    ss_between = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
    ss_total = ((grand - grand_mean) ** 2).sum()
    eta_sq = ss_between / ss_total if ss_total > 0 else 0.0
    return float(F), float(p), float(eta_sq)


def welch_ttest(a, b):
    """Welch's t-test (unequal variances). Returns (t, p)."""
    from scipy import stats

    t, p = stats.ttest_ind(np.asarray(a, dtype=float),
                           np.asarray(b, dtype=float),
                           equal_var=False)
    return float(t), float(p)


def describe(label, values):
    """One formatted summary line for a metric within one configuration."""
    v = np.asarray(values, dtype=float)
    return ("%-3s n=%3d  mean=%10.3f  sd=%9.3f  median=%9.3f  min=%9.3f  max=%10.3f"
            % (label, len(v), v.mean(), v.std(ddof=1), np.median(v), v.min(), v.max()))


def kruskal(groups):
    """Kruskal-Wallis rank test. Returns (H, p)."""
    from scipy import stats

    H, p = stats.kruskal(*[np.asarray(g, dtype=float) for g in groups])
    return float(H), float(p)


def mannwhitney(a, b):
    """Mann-Whitney U, two-sided. Returns (U, p)."""
    from scipy import stats

    U, p = stats.mannwhitneyu(np.asarray(a, dtype=float),
                              np.asarray(b, dtype=float),
                              alternative="two-sided")
    return float(U), float(p)


def report_anova(name, groups_by_label, log):
    """Run + print/log the full comparison of one metric across configs.

    Reports BOTH a mean-based ANOVA/Welch and a rank-based Kruskal-Wallis/
    Mann-Whitney. The rank tests are not optional decoration: the LigninGraphs
    cost metrics are heavily right-skewed and censored at i_max (median ~10-19 vs
    max 501), so a mean-based ANOVA there is dominated by a handful of censored
    runs and can report an effect the medians do not support. Read the rank test
    first whenever `n_cens` is non-zero or mean >> median.

    groups_by_label: list of (label, values) in CONFIGS order.
    log: callable taking one string.
    Returns (F, p_anova, eta^2, p_kruskal).
    """
    log("")
    log("--- %s ---" % name)
    for label, vals in groups_by_label:
        log("  " + describe(label, vals))

    F, p, eta = one_way_anova([v for _, v in groups_by_label])
    log("  ANOVA (means): F=%.3f  p=%.3g  eta^2=%.3f" % (F, p, eta))
    for i in range(len(groups_by_label)):
        for j in range(i + 1, len(groups_by_label)):
            li, vi = groups_by_label[i]
            lj, vj = groups_by_label[j]
            t, pt = welch_ttest(vi, vj)
            log("    %s vs %s: Welch t=%7.3f  p=%.3g" % (li, lj, t, pt))

    H, pk = kruskal([v for _, v in groups_by_label])
    log("  Kruskal-Wallis (ranks): H=%.3f  p=%.3g" % (H, pk))
    for i in range(len(groups_by_label)):
        for j in range(i + 1, len(groups_by_label)):
            li, vi = groups_by_label[i]
            lj, vj = groups_by_label[j]
            U, pu = mannwhitney(vi, vj)
            log("    %s vs %s: Mann-Whitney p=%.3g" % (li, lj, pu))

    return F, p, eta, pk
