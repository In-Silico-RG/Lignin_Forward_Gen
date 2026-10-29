#!/usr/bin/env python
"""
calibrate_kmc_ref.py

The fair comparison: calibrate Lignin-KMC (Orella et al. 2019, ligninkmc 0.2.2)
against the SAME reference set, objective, DP, seeds, replicate counts and
optimiser as our generator (code/calibrate_ref.py), with the SAME number of free
parameters.

Lignin-KMC's rates are Eyring rates from the Gani et al. DFT barriers. We keep
those barriers and fit one dimensionless MULTIPLIER per bond-forming reaction
class -- beta-5, beta-1, beta-beta, 5-5, 4-O-5 (ligninkmc's '5o4') and alpha-O-4
-- with beta-O-4 pinned at 1.0 as the scale, plus the monomer-addition rate
(GROW), which sets the supply regime (endwise vs bulk). Seven free parameters,
as in our model; oxidation and hydration rates stay at their defaults. Each
multiplier scales every substrate pair and size class of its reaction alike, so
the DFT pair specificity is preserved and only the class-to-class balance is
refit.

Lignin-KMC has no H monomer, so a sample's H fraction is ignored (noted in the
log). Observables are the same ratios to beta-O-4 as reference_set.py, with
DC = (2 n_55 + n_5o4 + n_b5) / n_units.

Written for the `lignin_kmc` environment (python 3.6): no f-strings.

    python code/calibrate_kmc_ref.py [--source-weight] [--hold-out ID ...] [--procs 24]
"""
from __future__ import print_function
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

from ligninkmc.create_lignin import (calc_rates, create_initial_monomers,
                                     create_initial_events, create_initial_state)
from ligninkmc.kmc_functions import run_kmc, analyze_adj_matrix
from ligninkmc.kmc_common import (Event, GROW, ADJ_MATRIX, BONDS, MONO_LIST,
                                  DEF_E_BARRIER_KCAL_MOL, BO4, BB, B5, B1, C5O4, AO4, C5C5)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reference_set as rs                    # noqa: E402

OUT = os.path.join(HERE, os.pardir, "results")
TEMP = 298.15
N = 120                       # DP, as in calibrate_ref / calibrate_full
INI_MONOS = 2
SIM_TIME = 1.0e9              # growth stops on DP, not on time
SEED_BASE = 7000
FIT_REPS = 8

# free parameters: multipliers on the DFT rates of six reaction classes (beta-O-4
# pinned at 1) and the monomer-addition rate. Bounds in log space are generous.
BOUNDS = [("m_b5", 1e-3, 1e3), ("m_b1", 1e-3, 1e3), ("m_bb", 1e-3, 1e3),
          ("m_55", 1e-3, 1e3), ("m_5o4", 1e-3, 1e3), ("m_ao4", 1e-3, 1e3),
          ("add_rate", 1e-2, 1e4)]
X0 = {"m_b5": 1.0, "m_b1": 1.0, "m_bb": 1.0, "m_55": 1.0, "m_5o4": 1.0, "m_ao4": 1.0,
      "add_rate": 1.0}
MULT_KEY = {B5: "m_b5", B1: "m_b1", BB: "m_bb", C5C5: "m_55", C5O4: "m_5o4", AO4: "m_ao4"}

BASE_RATES = calc_rates(TEMP, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)


def scaled_rates(params):
    rates = {}
    for rxn, subs in BASE_RATES.items():
        m = params.get(MULT_KEY.get(rxn, ""), 1.0)
        rates[rxn] = {sub: {sz: v * m for sz, v in sizes.items()} for sub, sizes in subs.items()}
    return rates


def unpack(x):
    params = {}
    for (name, lo, hi), v in zip(BOUNDS, np.exp(x)):
        params[name] = float(np.clip(v, lo, hi))
    return params


def run_one(args):
    """One polymerization; returns the bond counts and the unit count."""
    sg, pct_s, params, seed = args
    rates = scaled_rates(params)
    np.random.seed(seed)
    draw = np.random.rand(INI_MONOS)
    mons = create_initial_monomers(pct_s, draw)
    events = create_initial_events(mons, rates)
    state = create_initial_state(events, mons)
    events.append(Event(GROW, [], rate=params["add_rate"]))
    res = run_kmc(rates, state, events, n_max=N, t_max=SIM_TIME, random_seed=seed, sg_ratio=sg)
    bonds = analyze_adj_matrix(res[ADJ_MATRIX])[BONDS]
    return dict(bonds), len(res[MONO_LIST])


def observables(results):
    n_units = sum(n for _, n in results)
    b = {}
    for bonds, _ in results:
        for k, v in bonds.items():
            b[k] = b.get(k, 0) + v
    per = {"bo4": b.get(BO4, 0), "b5": b.get(B5, 0), "bb": b.get(BB, 0), "b1": b.get(B1, 0),
           "c55": b.get(C5C5, 0), "c4o5": b.get(C5O4, 0),
           "dc": 2 * b.get(C5C5, 0) + b.get(C5O4, 0) + b.get(B5, 0),
           "ao4": b.get(AO4, 0)}
    return {k: 100.0 * v / n_units for k, v in per.items()}


def seed_offset(row):
    return int(round(100 * row["sg_ratio"]))


def evaluate(params, rows, reps, pool):
    tasks, index = [], []
    for row in rows:
        pct_s = row["sg_ratio"] / (1.0 + row["sg_ratio"])
        for r in range(reps):
            tasks.append((row["sg_ratio"], pct_s, params, SEED_BASE + r + seed_offset(row)))
            index.append(row["sample_id"])
    outs = pool.map(run_one, tasks, chunksize=4) if pool else [run_one(t) for t in tasks]
    by = {}
    for sid, o in zip(index, outs):
        by.setdefault(sid, []).append(o)
    return {sid: observables(v) for sid, v in by.items()}


def total_chi2(obs_by, rows, weights):
    tot, per = 0.0, {}
    for row, w in zip(rows, weights):
        c, terms = rs.chi2_sample(row, obs_by[row["sample_id"]])
        per[row["sample_id"]] = (c, terms, obs_by[row["sample_id"]])
        tot += w * c
    return tot, per


def source_weights(rows, by_source):
    if not by_source:
        return [1.0] * len(rows)
    n = {}
    for r in rows:
        n[r["source"]] = n.get(r["source"], 0) + 1
    return [1.0 / n[r["source"]] for r in rows]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", type=int, default=FIT_REPS)
    ap.add_argument("--eval-reps", type=int, default=30)
    ap.add_argument("--source-weight", action="store_true")
    ap.add_argument("--hold-out", nargs="*", default=[])
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    rows_all = rs.load_reference_set()
    rows = [r for r in rows_all if r["sample_id"] not in args.hold_out]
    held = [r for r in rows_all if r["sample_id"] in args.hold_out]
    weights = source_weights(rows, args.source_weight)
    tag = args.tag or ("_sw" if args.source_weight else "") + \
        ("_ho-" + "-".join(args.hold_out) if args.hold_out else "")
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    log_path = os.path.join(OUT, "calibration_kmc_ref%s.log" % tag)
    lines = []

    def log(s=""):
        lines.append(s); print(s); sys.stdout.flush()

    log("Lignin-KMC calibration on the expanded reference set (ratios to beta-O-4)")
    log("ligninkmc rates from DEF_E_BARRIER_KCAL_MOL at %.2f K; free: six class multipliers "
        "(beta-O-4 pinned at 1) + add_rate" % TEMP)
    log("samples in fit: %d  held out: %s  reps: %d  N: %d  seed base: %d  source-weighted: %s  procs: %d"
        % (len(rows), args.hold_out or "-", args.reps, N, SEED_BASE, args.source_weight, args.procs))
    log("NOTE: Lignin-KMC has no H monomer; samples' h_frac is ignored (max in set: %.2f)"
        % max([r["h_frac"] or 0.0 for r in rows_all]))

    pool = Pool(args.procs) if args.procs > 1 else None
    state = {"n": 0, "best": np.inf, "t0": time.time()}

    def loss(x):
        obs_by = evaluate(unpack(x), rows, args.reps, pool)
        L, _ = total_chi2(obs_by, rows, weights)
        state["n"] += 1
        state["best"] = min(state["best"], L)
        if state["n"] % 10 == 0:
            log("  eval %3d  chi2=%.3f  best=%.3f  (%.0f s)" % (state["n"], L, state["best"], time.time() - state["t0"]))
        return L

    x0 = np.log([X0[name] for name, _, _ in BOUNDS])
    simplex = [x0] + [x0 + np.eye(len(x0))[i] * 0.7 for i in range(len(x0))]
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"maxiter": args.maxiter, "xatol": 1e-2, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    best = unpack(res.x)
    log("\nfitted params: %s" % {k: round(v, 4) for k, v in best.items()})
    log("minimised chi2 = %.3f (%d evals, %.0f s)" % (res.fun, state["n"], time.time() - state["t0"]))

    obs_by = evaluate(best, rows + held, args.eval_reps, pool)
    chi2, per = total_chi2(obs_by, rows, weights)
    log("\nchi2 at %d reps: %.3f over %d samples, %d terms"
        % (args.eval_reps, chi2, len(rows), sum(len(t) for _, t, _ in per.values())))
    log("%-34s %5s %6s  %s" % ("sample", "S/G", "chi2", "z per ratio (model - measured)/sigma"))
    for row in rows:
        c, terms, _ = per[row["sample_id"]]
        log("%-34s %5.2f %6.2f  %s" % (row["sample_id"], row["sg_ratio"], c,
                                      " ".join("%s:%+.2f" % (k, z) for k, z in terms.items())))
    if held:
        log("\nheld-out samples (predicted with the fitted multipliers):")
        for row in held:
            c, terms = rs.chi2_sample(row, obs_by[row["sample_id"]])
            log("%-34s %5.2f %6.2f  %s" % (row["sample_id"], row["sg_ratio"], c,
                                          " ".join("%s:%+.2f" % (k, z) for k, z in terms.items())))
    log("\nmodel per 100 units at the fitted multipliers (%d reps):" % args.eval_reps)
    log("%-34s %5s %6s %5s %5s %5s %5s %5s %5s %5s" % ("sample", "S/G", "bo4", "b5", "bb", "b1", "5-5", "4O5", "aO4", "DC"))
    for row in rows + held:
        o = obs_by[row["sample_id"]]
        log("%-34s %5.2f %6.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f" % (
            row["sample_id"], row["sg_ratio"], o["bo4"], o["b5"], o["bb"], o["b1"], o["c55"], o["c4o5"], o["ao4"], o["dc"]))
    if pool:
        pool.close()

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, "fitted_params_kmc_ref%s.json" % tag), "w") as fh:
        d = {"_note": "Lignin-KMC multipliers fitted on data/reference_inventories.csv; see calibration_kmc_ref%s.log" % tag,
             "_held_out": args.hold_out, "_source_weighted": args.source_weight, "_chi2_eval": chi2}
        d.update(best)
        json.dump(d, fh, indent=1)
    print("Wrote", os.path.relpath(log_path))


if __name__ == "__main__":
    main()
