#!/usr/bin/env python
"""
cross_ligningraphs_incremental.py

Cross-algorithm validation, arm 3 of 3: LigninGraphs (Wang et al. 2022) on the
matched SW/HW/GR parameterization (manuscript section 3.1).

LigninGraphs is a BACKWARD generator of the Metropolis kind: it is handed a
target linkage + monomer distribution and walks a graph, accepting/rejecting
SINGLE MONOMER ADDITIONS by a Metropolis criterion on the distance between the
current structure's metrics and the target. Like Lignin-KMC and unlike LGS, it
builds ONE structure incrementally and never discards a partial polymer -- so its
cost is again not the same unit of work as LGS's whole-structure rejection.

We record:
  wall_s -- wall-clock seconds per structure (the dead manuscript's t_gen proxy)
  steps  -- Metropolis steps taken (the mechanism-native unit of work)

RECONSTRUCTION NOTICE
---------------------
The original script was lost with a deleted session scratchpad; only the docs/03
summary survived ("NULL -- no significant difference (p=0.58 t_gen, p=0.63
steps); huge variance; the clean n=5 pilot pattern vanished at n=30"). This is a
re-implementation from that summary, not a byte-recovery.

Run inside the `ligningraphs` conda env (numpy<2, networkx<3, pysmiles==1.1.2):
    conda run -n ligningraphs python code/cross_ligningraphs_incremental.py
"""
import os
import sys
import time

import numpy as np

from ligning.optimization import Trajectory

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_common import (CONFIGS, DP, N_REPS, BASE_SEED, sg_to_pct_s,
                          results_dir, write_csv, report_anova, target_profile)
from repo_paths import rel

# ligning.rules.linkage_names order -- the Trajectory linkage vector must follow it
LINKAGE_NAMES = ['4-O-5', 'alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1']
MONOMER_NAMES = ['H', 'G', 'S']

TMETRO = 30      # docs/03: the matched incremental runs used Tmetro = 30
I_MAX = 500      # LigninGraphs' own default inner-loop cap

# Map cross_common's linkage keys onto ligning's names, so both backward arms are
# handed the SAME per-archetype target (cross_common.target_profile).
# alpha-O-4 is set to 0: its Lignin-KMC barrier is 20.7 kcal/mol for every pair
# (docs/04), i.e. effectively never formed, so giving it target weight would ask
# LigninGraphs to hit a linkage the chemistry excludes.
COMMON_TO_LIGNING = {
    "bo4": 'beta-O-4', "b5": 'beta-5', "b1": 'beta-1',
    "bb": 'beta-beta', "55": '5-5', "4o5": '4-O-5',
}


def linkage_vector(target_bo4, sg_ratio):
    """Target linkage distribution in ligning's canonical order.

    Uses the per-archetype profile: an earlier version applied one
    softwood-derived residual to every config, which asks S-rich HW for ~39% of
    its non-beta-O-4 mass as 5-5 -- impossible, since S blocks C5. See
    cross_common.target_profile and docs/03.
    """
    prof = target_profile(target_bo4, sg_ratio)
    d = dict((COMMON_TO_LIGNING[k], v) for k, v in prof.items())
    d['alpha-O-4'] = 0.0
    return [d[name] for name in LINKAGE_NAMES]


def monomer_vector(sg_ratio):
    """Target monomer distribution [H, G, S].

    H is held at 0 for all three configs. The matched parameterization specifies
    only S/G, and giving grass H units here would make GR the one arm running
    different chemistry -- docs/02 section F is explicit that every generator's
    grass model is incomplete, which is a separate finding, not something to
    smuggle into this comparison.
    """
    pct_s = sg_to_pct_s(sg_ratio)
    return [0.0, 1.0 - pct_s, pct_s]


def run_one(sg_ratio, target_bo4, seed):
    """One incremental Metropolis polymerization. Returns (wall_s, steps, size)."""
    traj = Trajectory(
        linkage_distribution_input=linkage_vector(target_bo4, sg_ratio),
        monomer_distribution_input=monomer_vector(sg_ratio),
        Tmetro=TMETRO,
        expected_size=DP,
        max_size=DP,
        distribution_scaling=1.0,
        verbose=False,
    )
    t0 = time.time()
    polymer, distance, monomer_count, i_step = traj.run_MCMC(rseed=seed, i_max=I_MAX)
    wall = time.time() - t0
    return wall, int(i_step), int(monomer_count)


def main():
    out_dir = results_dir()
    log_path = os.path.join(out_dir, "cross_ligningraphs_incremental.log")
    lines = []

    def log(s):
        lines.append(s)
        print(s)

    log("Cross-algorithm validation -- LigninGraphs, matched params (incremental)")
    log("RECONSTRUCTION of a lost script; see docstring + docs/03.")
    log("params: DP=%d Tmetro=%d i_max=%d n_reps=%d base_seed=%d"
        % (DP, TMETRO, I_MAX, N_REPS, BASE_SEED))
    log("")
    log("%-3s %6s %9s %12s %13s %12s"
        % ("cfg", "S/G", "target", "wall_s(mean)", "steps(mean)", "size(mean)"))

    rows = []
    by_label_wall = []
    by_label_steps = []

    for ci, (label, sg, target) in enumerate(CONFIGS):
        walls, steps, sizes = [], [], []
        for rep in range(N_REPS):
            seed = BASE_SEED + 1000 * ci + rep
            w, s, size = run_one(sg, target, seed)
            walls.append(w)
            steps.append(s)
            sizes.append(size)
            rows.append({
                "config": label, "sg_ratio": sg, "target_bo4": target,
                "rep": rep, "seed": seed, "wall_s": w, "steps": s,
                "monomer_count": size,
            })
        by_label_wall.append((label, walls))
        by_label_steps.append((label, steps))
        log("%-3s %6.2f %9.2f %12.4f %13.1f %12.1f"
            % (label, sg, target, np.mean(walls), np.mean(steps), np.mean(sizes)))

    write_csv(os.path.join(out_dir, "cross_ligningraphs_incremental.csv"),
              ["config", "sg_ratio", "target_bo4", "rep", "seed",
               "wall_s", "steps", "monomer_count"],
              rows)

    log("")
    log("=" * 72)
    log("STATISTICS")
    report_anova("wall-clock seconds (the t_gen proxy)", by_label_wall, log)
    report_anova("Metropolis steps (native unit of work)", by_label_steps, log)

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nWrote %s" % rel(os.path.join(out_dir, "cross_ligningraphs_incremental.csv")))
    print("Wrote %s" % rel(log_path))


if __name__ == "__main__":
    main()
