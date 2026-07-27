#!/usr/bin/env python
"""
supply_regime_sweep.py

Follow-up to forward_emergence_sweep.py. The emergence sweep (pure endwise,
add_rate=1) reproduced the S/G->beta-O-4 relationship for hardwoods but
under-produced softwood beta-O-4 and grossly over-produced beta-5, with the
condensed/branch linkages (5-5, 4-O-5) essentially absent.

This probes whether that is a REGIME artifact or a PARAMETER fact: vary the
monomer-addition rate (the endwise <-> bulk supply knob) over many orders of
magnitude at fixed S/G. Low add_rate = endwise (monomer adds onto oligomer);
high add_rate = bulk (all monomers present at once, monomer-monomer and
oligomer-oligomer coupling enabled -> can form 5-5 / 4-O-5 branching).

Run: conda run -n lignin_kmc python code/supply_regime_sweep.py
Outputs: results/supply_regime_sweep.csv + printed table.
"""
import os
import sys

import numpy as np

from ligninkmc.create_lignin import (calc_rates, create_initial_monomers,
                                      create_initial_events, create_initial_state)
from ligninkmc.kmc_functions import run_kmc, analyze_adj_matrix
from ligninkmc.kmc_common import (Event, GROW, ADJ_MATRIX, BONDS, CHAIN_BRANCH_COEFF,
                                  BO4, BB, B5, B1, C5O4, AO4, C5C5, DEF_E_BARRIER_KCAL_MOL)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import rel

TEMP = 298.15
MAX_MONOS = 100
INI_MONOS = 2
SIM_TIME = 1.0e12
N_REPS = 15
BASE_SEED = 20260714

# fixed compositions: softwood (pure G) and a hardwood (S/G=1)
SG_POINTS = [0.0, 1.0]
# supply knob spanning strongly-endwise -> strongly-bulk
ADD_RATES = [1.0, 1e2, 1e4, 1e6, 1e8]

BOND_KEYS = [BO4, BB, B5, B1, C5O4, AO4, C5C5]


def run_one(sg_ratio, add_rate, rates, seed):
    pct_s = sg_ratio / (1.0 + sg_ratio)
    np.random.seed(seed)
    mons = create_initial_monomers(pct_s, np.random.rand(INI_MONOS))
    evs = create_initial_events(mons, rates)
    st = create_initial_state(evs, mons)
    if MAX_MONOS > INI_MONOS:
        evs.append(Event(GROW, [], rate=add_rate))
    res = run_kmc(rates, st, evs, n_max=MAX_MONOS, t_max=SIM_TIME,
                  random_seed=seed, sg_ratio=sg_ratio)
    summ = analyze_adj_matrix(res[ADJ_MATRIX])
    bonds = summ[BONDS]
    tot = sum(bonds.values())
    fr = {k: (bonds[k] / tot if tot else 0.0) for k in BOND_KEYS}
    bc = summ[CHAIN_BRANCH_COEFF]
    branch = float(np.mean(list(bc.values()))) if bc else 0.0
    return fr, tot, branch


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "supply_regime_sweep.csv")
    rates = calc_rates(TEMP, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

    print("Supply-regime sweep (endwise <-> bulk) at fixed S/G")
    print("params: max_monos=%d n_reps=%d  bonds shown as %% of total\n" % (MAX_MONOS, N_REPS))
    print("%4s %9s %7s %6s %6s %6s %6s %6s %7s %7s"
          % ("S/G", "add_rate", "bO4", "b5", "bb", "55", "5o4", "b1", "branch", "nbond"))

    rows = []
    for sg in SG_POINTS:
        for ar in ADD_RATES:
            acc = {k: [] for k in BOND_KEYS}
            branches, nbonds = [], []
            for rep in range(N_REPS):
                seed = BASE_SEED + rep + int(sg * 100) + int(np.log10(ar) * 7)
                fr, tot, branch = run_one(sg, ar, rates, seed)
                for k in BOND_KEYS:
                    acc[k].append(fr[k])
                branches.append(branch); nbonds.append(tot)
            m = {k: np.mean(acc[k]) for k in BOND_KEYS}
            row = {"sg_ratio": sg, "add_rate": ar,
                   "bo4": m[BO4], "b5": m[B5], "bb": m[BB], "c55": m[C5C5],
                   "c5o4": m[C5O4], "b1": m[B1], "ao4": m[AO4],
                   "branch": np.mean(branches), "n_bonds": np.mean(nbonds)}
            rows.append(row)
            print("%4.1f %9.0e %6.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %7.3f %7.1f"
                  % (sg, ar, 100 * m[BO4], 100 * m[B5], 100 * m[BB], 100 * m[C5C5],
                     100 * m[C5O4], 100 * m[B1], row["branch"], row["n_bonds"]))

    cols = ["sg_ratio", "add_rate", "bo4", "b5", "bb", "c55", "c5o4", "b1", "ao4",
            "branch", "n_bonds"]
    with open(out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join("%.5f" % r[c] for c in cols) + "\n")
    print("\nWrote %s" % rel(out_csv))


if __name__ == "__main__":
    main()
