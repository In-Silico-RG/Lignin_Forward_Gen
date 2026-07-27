#!/usr/bin/env python
"""
forward_emergence_sweep.py

Central experiment for the "forward generation with the right chemistry" paper.

Question: does Lignin-KMC (a FORWARD, mechanistic generator driven only by
DFT-derived radical-coupling rate constants) reproduce the real S/G -> beta-O-4
relationship and the softwood beta-O-4 "ceiling" as EMERGENT properties --
i.e. WITHOUT being told any target bond fraction?

Method: sweep the syringyl fraction from pure-G (softwood) to high-S (hardwood),
run n replicate kinetic Monte Carlo polymerizations at each point under an
endwise (slow monomer supply) regime, and measure the emergent bond distribution.

Run inside the `lignin_kmc` conda env (python 3.6, ligninkmc 0.2.2):
    conda run -n lignin_kmc python code/forward_emergence_sweep.py

Outputs: results/emergence_sweep.csv  and a printed summary table.
"""
import os
import numpy as np

from ligninkmc.create_lignin import (calc_rates, create_initial_monomers,
                                      create_initial_events, create_initial_state)
from ligninkmc.kmc_functions import run_kmc, analyze_adj_matrix
from ligninkmc.kmc_common import (Event, GROW, ADJ_MATRIX, BONDS, CHAIN_LEN,
                                  CHAIN_BRANCH_COEFF, BO4, BB, B5, B1, C5O4, AO4,
                                  C5C5, DEF_E_BARRIER_KCAL_MOL)

# ---- experiment parameters -------------------------------------------------
TEMP = 298.15          # K  (matches the DFT barrier reference state)
MAX_MONOS = 100        # monomers per polymer
INI_MONOS = 2          # seed monomers present at t=0
ADD_RATE = 1.0         # monomers/s -- SLOW vs coupling rates => endwise regime
SIM_TIME = 1.0e9       # s  -- large so growth is governed by MAX_MONOS, not time
N_REPS = 20
BASE_SEED = 20260714

# S/G ratios spanning softwood -> hardwood. sg=0 is pure guaiacyl (softwood).
SG_RATIOS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]

# Literature anchors (whole-cell-wall HSQC0, Bourmaud/Luterbacher 2024) for context
NMR_ANCHORS = "SW black pine ~48%; poplar ~60%; birch ~72%; high-S poplar ~79%"

BOND_KEYS = [BO4, BB, B5, B1, C5O4, AO4, C5C5]
C5_LINKAGES = [B5, C5O4, C5C5]   # linkages that REQUIRE a free C5 (G/H only)


def run_one(sg_ratio, rates, seed):
    """Run a single KMC polymerization; return (bond_fractions, n_bonds, branch_coeff)."""
    pct_s = sg_ratio / (1.0 + sg_ratio)
    np.random.seed(seed)
    monomer_draw = np.random.rand(INI_MONOS)
    initial_monomers = create_initial_monomers(pct_s, monomer_draw)
    initial_events = create_initial_events(initial_monomers, rates)
    initial_state = create_initial_state(initial_events, initial_monomers)
    if MAX_MONOS > INI_MONOS:
        initial_events.append(Event(GROW, [], rate=ADD_RATE))
    result = run_kmc(rates, initial_state, initial_events, n_max=MAX_MONOS,
                     t_max=SIM_TIME, random_seed=seed, sg_ratio=sg_ratio)
    summary = analyze_adj_matrix(result[ADJ_MATRIX])
    bonds = summary[BONDS]
    total = sum(bonds.values())
    fracs = {k: (bonds[k] / total if total else 0.0) for k in BOND_KEYS}
    # mean branching coefficient across oligomers (0 = perfectly linear)
    bc = summary[CHAIN_BRANCH_COEFF]
    branch = float(np.mean(list(bc.values()))) if bc else 0.0
    return fracs, total, branch


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "emergence_sweep.csv")

    rates = calc_rates(TEMP, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

    rows = []
    print("Lignin-KMC forward emergence sweep (endwise regime)")
    print("params: max_monos=%d ini=%d add_rate=%g temp=%.2fK n_reps=%d"
          % (MAX_MONOS, INI_MONOS, ADD_RATE, TEMP, N_REPS))
    print("NMR anchors: %s\n" % NMR_ANCHORS)
    header = ("S/G", "pct_S", "bO4_mean", "bO4_sd", "C5link_mean", "bb_mean",
              "b5_mean", "55_mean", "branch", "n_bonds")
    print("%5s %6s %9s %7s %11s %8s %8s %7s %7s %8s" % header)

    for sg in SG_RATIOS:
        b04s, c5s, bbs, b5s, c55s, branches, nbonds = [], [], [], [], [], [], []
        for rep in range(N_REPS):
            seed = BASE_SEED + rep + int(sg * 100)
            fracs, total, branch = run_one(sg, rates, seed)
            b04s.append(fracs[BO4])
            c5s.append(sum(fracs[k] for k in C5_LINKAGES))
            bbs.append(fracs[BB]); b5s.append(fracs[B5]); c55s.append(fracs[C5C5])
            branches.append(branch); nbonds.append(total)
        row = {
            "sg_ratio": sg, "pct_s": sg / (1 + sg),
            "bo4_mean": np.mean(b04s), "bo4_sd": np.std(b04s),
            "c5link_mean": np.mean(c5s), "bb_mean": np.mean(bbs),
            "b5_mean": np.mean(b5s), "c55_mean": np.mean(c55s),
            "branch_mean": np.mean(branches), "n_bonds_mean": np.mean(nbonds),
        }
        rows.append(row)
        print("%5.2f %6.2f %8.1f%% %6.1f %10.1f%% %7.1f%% %7.1f%% %6.1f%% %7.3f %8.1f"
              % (sg, row["pct_s"], 100 * row["bo4_mean"], 100 * row["bo4_sd"],
                 100 * row["c5link_mean"], 100 * row["bb_mean"], 100 * row["b5_mean"],
                 100 * row["c55_mean"], row["branch_mean"], row["n_bonds_mean"]))

    # write CSV
    cols = ["sg_ratio", "pct_s", "bo4_mean", "bo4_sd", "c5link_mean",
            "bb_mean", "b5_mean", "c55_mean", "branch_mean", "n_bonds_mean"]
    with open(out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join("%.5f" % r[c] for c in cols) + "\n")
    print("\nWrote %s" % out_csv)


if __name__ == "__main__":
    main()
