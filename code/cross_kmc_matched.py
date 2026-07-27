#!/usr/bin/env python
"""
cross_kmc_matched.py

Cross-algorithm validation, arm 2 of 3: Lignin-KMC on the matched SW/HW/GR
parameterization (manuscript section 3.1).

THE MECHANISM MISMATCH (this is a finding, not a nuisance)
----------------------------------------------------------
LGS rejects WHOLE candidate structures: one "run" = many complete structures
built and discarded, so its cost metric counts structures. Lignin-KMC builds ONE
structure INCREMENTALLY by firing kinetic events and never discards a partial
polymer. The two do not share a unit of work, so per-structure wall-clock is not
an apples-to-apples comparison between them. Worse, Lignin-KMC is a FORWARD
generator: it has no beta-O-4 target to hit, so the dead manuscript's targets
(54/69/77%) cannot even be supplied to it. "Matching" the parameterization can
therefore only mean matching S/G and DP. That asymmetry is exactly why a
generation-cost metric cannot be compared across generator classes -- which is
the point of section 3.1.

We record two cost metrics:
  wall_s  -- wall-clock seconds per structure (the dead manuscript's t_gen proxy)
  events  -- number of KMC events fired (the mechanism-native unit of work)

RECONSTRUCTION NOTICE
---------------------
The original script was lost with a deleted session scratchpad; only the docs/03
summary survived ("wall-clock REVERSED (HW fastest); KMC event-count partially
supportive (HW & GR both > SW, p<1e-4, but HW indistinguishable from GR -- the
HW-unique peak NOT reproduced)"). This is a re-implementation from that summary,
not a byte-recovery. Report what it produces.

RATE-SET PROVENANCE (docs/04): ligninkmc ships two rate sets -- DEF_RXN_RATES and
MANUSCRIPT_RATES, which differ by ~150x. We derive rates from
DEF_E_BARRIER_KCAL_MOL via calc_rates() at 298.15 K, matching
forward_emergence_sweep.py. Any reuse of these timings MUST state this.

Run inside the `lignin_kmc` conda env (python 3.6, ligninkmc 0.2.2):
    conda run -n lignin_kmc python code/cross_kmc_matched.py
"""
import os
import sys
import time

import numpy as np

from ligninkmc.create_lignin import (calc_rates, create_initial_monomers,
                                     create_initial_events, create_initial_state)
from ligninkmc.kmc_functions import run_kmc, analyze_adj_matrix
from ligninkmc.kmc_common import (Event, GROW, ADJ_MATRIX, BONDS, TIME,
                                  BO4, DEF_E_BARRIER_KCAL_MOL)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_common import (CONFIGS, DP, N_REPS, BASE_SEED, sg_to_pct_s,
                          results_dir, write_csv, report_anova)

TEMP = 298.15        # K -- the DFT barriers' reference state
INI_MONOS = 2
ADD_RATE = 1.0       # monomers/s -- slow vs coupling => endwise regime
SIM_TIME = 1.0e9     # s -- large, so growth stops on DP, not on time


def run_one(pct_s, sg_ratio, rates, seed):
    """One incremental KMC polymerization to DP monomers.

    Returns (wall_seconds, n_events, bo4_fraction).
    """
    np.random.seed(seed)
    monomer_draw = np.random.rand(INI_MONOS)

    t0 = time.time()
    initial_monomers = create_initial_monomers(pct_s, monomer_draw)
    initial_events = create_initial_events(initial_monomers, rates)
    initial_state = create_initial_state(initial_events, initial_monomers)
    if DP > INI_MONOS:
        initial_events.append(Event(GROW, [], rate=ADD_RATE))
    result = run_kmc(rates, initial_state, initial_events, n_max=DP,
                     t_max=SIM_TIME, random_seed=seed, sg_ratio=sg_ratio)
    wall = time.time() - t0

    # run_kmc appends one entry to its TIME list per fired event, so the length
    # of that list is the mechanism-native event count.
    n_events = len(result[TIME])

    summary = analyze_adj_matrix(result[ADJ_MATRIX])
    bonds = summary[BONDS]
    total = sum(bonds.values())
    f_bo4 = bonds[BO4] / float(total) if total else 0.0
    return wall, n_events, f_bo4


def main():
    out_dir = results_dir()
    log_path = os.path.join(out_dir, "cross_kmc_matched.log")
    lines = []

    def log(s):
        lines.append(s)
        print(s)

    rates = calc_rates(TEMP, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

    log("Cross-algorithm validation -- Lignin-KMC, matched params (incremental)")
    log("RECONSTRUCTION of a lost script; see docstring + docs/03.")
    log("params: DP=%d ini=%d add_rate=%g temp=%.2fK n_reps=%d"
        % (DP, INI_MONOS, ADD_RATE, TEMP, N_REPS))
    log("rates: calc_rates(DEF_E_BARRIER_KCAL_MOL) -- NOT MANUSCRIPT_RATES (docs/04)")
    log("NOTE: forward generator -- the 54/69/77%% beta-O-4 targets CANNOT be")
    log("      supplied; only S/G and DP are matched. See docstring.")
    log("")
    log("%-3s %6s %12s %14s %12s" % ("cfg", "S/G", "wall_s(mean)", "events(mean)", "bo4(mean)"))

    rows = []
    by_label_wall = []
    by_label_events = []

    for ci, (label, sg, target) in enumerate(CONFIGS):
        pct_s = sg_to_pct_s(sg)
        walls, evs, bo4s = [], [], []
        for rep in range(N_REPS):
            seed = BASE_SEED + 1000 * ci + rep
            w, n_ev, f_bo4 = run_one(pct_s, sg, rates, seed)
            walls.append(w)
            evs.append(n_ev)
            bo4s.append(f_bo4)
            rows.append({
                "config": label, "sg_ratio": sg, "target_bo4": target,
                "rep": rep, "seed": seed, "wall_s": w, "events": n_ev,
                "bo4_frac": f_bo4,
            })
        by_label_wall.append((label, walls))
        by_label_events.append((label, evs))
        log("%-3s %6.2f %12.4f %14.1f %12.3f"
            % (label, sg, np.mean(walls), np.mean(evs), np.mean(bo4s)))

    write_csv(os.path.join(out_dir, "cross_kmc_matched.csv"),
              ["config", "sg_ratio", "target_bo4", "rep", "seed",
               "wall_s", "events", "bo4_frac"],
              rows)

    log("")
    log("=" * 72)
    log("STATISTICS")
    report_anova("wall-clock seconds (the t_gen proxy)", by_label_wall, log)
    report_anova("KMC events fired (native unit of work)", by_label_events, log)

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nWrote %s" % os.path.join(out_dir, "cross_kmc_matched.csv"))
    print("Wrote %s" % log_path)


if __name__ == "__main__":
    main()
