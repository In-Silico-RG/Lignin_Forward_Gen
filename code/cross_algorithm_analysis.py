#!/usr/bin/env python
"""
cross_algorithm_analysis.py

Cross-algorithm validation, aggregation + audit (manuscript section 3.1).

Reads whatever the three arms wrote into results/ and produces:
  1. a per-arm, per-metric summary with one-way ANOVA + pairwise Welch tests;
  2. an AUDIT of every claim recorded in docs/03 against what the rewritten
     scripts actually produce.

Point (2) is the reason this script exists. The original scripts were lost with a
deleted session scratchpad, so docs/03's numbers are currently unbacked by any
code in this repo -- and manuscript section 3.1 cites them. This script states,
per claim, whether the re-implementation reproduces it. Where it does not, the
manuscript cannot keep citing the archived number.

Runs on the base env (numpy + scipy). The three arms must be run first:
    python code/cross_rejection_port.py
    conda run -n lignin_kmc python code/cross_kmc_matched.py
    conda run -n ligningraphs python code/cross_ligningraphs_incremental.py
    python code/cross_algorithm_analysis.py
"""
import os

import numpy as np

from cross_common import LABELS, results_dir, report_anova


def read_csv(path):
    """Minimal CSV reader -> dict of column -> list (floats where possible)."""
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    cols = lines[0].split(",")
    data = dict((c, []) for c in cols)
    for line in lines[1:]:
        for c, v in zip(cols, line.split(",")):
            try:
                data[c].append(float(v))
            except ValueError:
                data[c].append(v)
    return data


def group_by_config(data, metric):
    """[(label, values)] in LABELS order."""
    out = []
    for label in LABELS:
        vals = [v for cfg, v in zip(data["config"], data[metric]) if cfg == label]
        out.append((label, vals))
    return out


# Each entry: (arm, csv basename, [(metric, human name)], archived claim text)
ARMS = [
    ("LGS (enumeration)", "cross_lgs_matched.csv",
     [("structures", "structures enumerated"), ("wall_s", "wall-clock seconds")],
     "docs/03: t_gen shows a non-monotonic peak at HARDWOOD (HW >> GR > SW). "
     "NOTE: LGS enumerates orderings of a target multiset -- it does not reject "
     "on bond profile (docs/03's description of it is wrong; see "
     "cross_lgs_matched.py). Its cost is expected to track the TARGET "
     "composition, not S/G. n=3 and it is not deterministic, so the stats below "
     "are indicative only; read the ordering."),
    ("Lignin-KMC (incremental)", "cross_kmc_matched.csv",
     [("wall_s", "wall-clock seconds"), ("events", "KMC events fired")],
     "docs/03: wall-clock REVERSED (HW fastest); event-count HW & GR both > SW "
     "(p<1e-4) but HW indistinguishable from GR -> HW-unique peak NOT reproduced."),
    ("LigninGraphs (incremental)", "cross_ligningraphs_incremental.csv",
     [("wall_s", "wall-clock seconds"), ("steps", "Metropolis steps")],
     "docs/03: NULL -- no significant difference (p=0.58 t_gen, p=0.63 steps); "
     "huge variance; the clean n=5 pilot pattern vanished at n=30."),
]


def main():
    out_dir = results_dir()
    log_path = os.path.join(out_dir, "cross_algorithm_analysis.log")
    lines = []

    def log(s):
        lines.append(s)
        print(s)

    log("=" * 74)
    log("CROSS-ALGORITHM VALIDATION -- aggregate analysis (manuscript section 3.1)")
    log("=" * 74)
    log("")
    log("These are RECONSTRUCTED scripts (originals lost with a deleted session")
    log("scratchpad). Every number below is what the rewritten code produces now.")
    log("The audit at the end compares them to the archived docs/03 claims.")

    findings = []

    for arm, basename, metrics, claim in ARMS:
        path = os.path.join(out_dir, basename)
        data = read_csv(path)
        log("")
        log("=" * 74)
        log("ARM: %s" % arm)
        log("  source: results/%s" % basename)
        log("  archived claim -- %s" % claim)
        if data is None:
            log("  !! NOT RUN (missing %s) -- run the arm first" % basename)
            continue
        for metric, human in metrics:
            groups = group_by_config(data, metric)
            F, p, eta, pk = report_anova(human, groups, log)
            meds = dict((lab, float(np.median(v))) for lab, v in groups)
            findings.append((arm, human, meds, p, pk))

    # ---- audit ------------------------------------------------------------
    log("")
    log("=" * 74)
    log("AUDIT -- reconstruction vs. archived claim")
    log("=" * 74)
    log("")
    log("Ordering is by MEDIAN, and the verdict follows the RANK test -- these")
    log("cost distributions are skewed and censored, so means mislead (see")
    log("cross_common.report_anova).")
    log("")
    log("%-28s %-26s %14s %10s %10s  %s"
        % ("arm", "metric", "ordering(med)", "p_anova", "p_rank", "verdict"))
    for arm, human, meds, p, pk in findings:
        order = " > ".join(sorted(meds, key=lambda k: -meds[k]))
        sig = "sig" if pk < 0.05 else "NULL"
        flag = "  <- means/ranks DISAGREE" if (p < 0.05) != (pk < 0.05) else ""
        log("%-28s %-26s %14s %10.3g %10.3g  %s%s"
            % (arm[:28], human[:26], order, p, pk, sig, flag))

    log("")
    log("Read this table against the archived claims printed per arm above.")
    log("A claim is only safe to cite in the manuscript if the reconstruction")
    log("reproduces both its ORDERING and its significance.")

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nWrote %s" % log_path)


if __name__ == "__main__":
    main()
