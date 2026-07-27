#!/usr/bin/env python
"""
cross_wallclock_stability.py

Is wall-clock generation time (the dead manuscript's t_gen proxy) a property of
the chemistry, or of the machine?

This is a controlled repeat of the Lignin-KMC arm: the SAME seeded, deterministic
code, run R times over. The kinetic trajectories are identical every time -- same
seeds, same events, same structures. Anything that moves between repeats is
measurement noise, not chemistry.

We contrast two cost metrics per repeat:
  events  -- KMC events fired (the mechanism-native unit of work)
  wall_s  -- wall-clock seconds (t_gen, the metric the dead manuscript used)

and ask, per repeat, whether the SW/HW/GR comparison would have been called
significant. If the verdict flips across repeats of identical code, then t_gen
cannot support a claim about lignin.

Motivation: two runs of cross_kmc_matched.py on this machine gave wall-clock
rank-test p=0.155 (null) and p=0.000375 (significant, softwood slowest) while the
event counts were bit-identical. This script tests whether that was a fluke.

Run inside the `lignin_kmc` conda env (python 3.6, ligninkmc 0.2.2):
    conda run -n lignin_kmc python code/cross_wallclock_stability.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_common import (CONFIGS, N_REPS, BASE_SEED, sg_to_pct_s,
                          results_dir, write_csv, kruskal)
from cross_kmc_matched import run_one, TEMP

from ligninkmc.create_lignin import calc_rates
from ligninkmc.kmc_common import DEF_E_BARRIER_KCAL_MOL

N_REPEATS = 5


def main():
    out_dir = results_dir()
    log_path = os.path.join(out_dir, "cross_wallclock_stability.log")
    lines = []

    def log(s):
        lines.append(s)
        print(s)

    rates = calc_rates(TEMP, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

    log("Wall-clock stability -- identical seeded code, repeated %d times" % N_REPEATS)
    log("Same seeds => identical trajectories. Anything that moves is the machine.")
    log("")

    rows = []
    wall_ps, event_ps = [], []

    log("%-7s | %-34s | %s" % ("", "wall-clock seconds (t_gen)", "KMC events fired"))
    log("%-7s | %8s %8s %8s %7s | %6s %6s %6s %10s"
        % ("repeat", "SW", "HW", "GR", "p", "SW", "HW", "GR", "p"))

    for rep_i in range(N_REPEATS):
        walls_by_cfg, events_by_cfg = [], []
        for ci, (label, sg, target) in enumerate(CONFIGS):
            pct_s = sg_to_pct_s(sg)
            walls, evs = [], []
            for rep in range(N_REPS):
                seed = BASE_SEED + 1000 * ci + rep
                w, n_ev, f_bo4 = run_one(pct_s, sg, rates, seed)
                walls.append(w)
                evs.append(n_ev)
                rows.append({"repeat": rep_i, "config": label, "rep": rep,
                             "seed": seed, "wall_s": w, "events": n_ev})
            walls_by_cfg.append(walls)
            events_by_cfg.append(evs)

        _, pw = kruskal(walls_by_cfg)
        _, pe = kruskal(events_by_cfg)
        wall_ps.append(pw)
        event_ps.append(pe)
        log("%-7d | %8.4f %8.4f %8.4f %7.4f | %6.1f %6.1f %6.1f %10.3g"
            % (rep_i,
               np.mean(walls_by_cfg[0]), np.mean(walls_by_cfg[1]), np.mean(walls_by_cfg[2]), pw,
               np.mean(events_by_cfg[0]), np.mean(events_by_cfg[1]), np.mean(events_by_cfg[2]), pe))

    write_csv(os.path.join(out_dir, "cross_wallclock_stability.csv"),
              ["repeat", "config", "rep", "seed", "wall_s", "events"], rows)

    log("")
    log("=" * 72)
    wall_sig = ["sig" if p < 0.05 else "NULL" for p in wall_ps]
    event_sig = ["sig" if p < 0.05 else "NULL" for p in event_ps]
    log("wall-clock verdict per repeat: %s" % ", ".join(wall_sig))
    log("  p range: %.3g to %.3g  (%.0fx spread)"
        % (min(wall_ps), max(wall_ps), max(wall_ps) / max(min(wall_ps), 1e-300)))
    log("events    verdict per repeat: %s" % ", ".join(event_sig))
    log("  p range: %.3g to %.3g" % (min(event_ps), max(event_ps)))
    log("")
    if len(set(wall_sig)) > 1:
        log("=> t_gen FLIPS its verdict across repeats of identical, seeded code.")
        log("   The structures are the same every time; only the clock moved. A")
        log("   wall-clock 'complexity' signal is therefore not a property of the")
        log("   chemistry, and cannot support a claim about lignin.")
    else:
        log("=> t_gen did not flip in %d repeats (verdict: %s). It is still a")
        log("   machine measurement, but this run does not demonstrate instability."
            % (N_REPEATS, wall_sig[0]))
    if len(set(event_sig)) == 1:
        log("=> The event count, by contrast, is stable: identical counts and the")
        log("   same verdict every repeat -- it measures the algorithm, not the host.")

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nWrote %s" % os.path.join(out_dir, "cross_wallclock_stability.csv"))
    print("Wrote %s" % log_path)


if __name__ == "__main__":
    main()
