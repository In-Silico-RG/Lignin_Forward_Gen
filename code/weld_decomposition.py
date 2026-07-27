"""
weld_decomposition.py

How much of the model's condensed (5-5 / 4-O-5) linkage is structurally forced?

The generator emits ONE molecule. Only 5-5 and 4-O-5 join two fragments that are
both larger than a monomer, so welding the n_nuclei separately grown chains into a
single macromolecule costs exactly n_nuclei - 1 condensed bonds, whatever the
coupling constants say. Those bonds are not a fitted preference for condensation;
they are the price of the single-molecule construction.

Every other condensed bond attaches a FREE MONOMER through its C5. That one was
elective: a beta-O-4 chain extension could have consumed the same monomer instead,
so its frequency is set by k_55 / k_4o5 against the beta constants.

This script separates the two. At each 5-5 / 4-O-5 event it records the sizes of
the two fragments being joined:

    both > 1 unit  -> WELD      (structurally required; only this reaction can do it)
    either == 1    -> ELECTIVE  (a monomer attachment the beta reactions could have made)

and reports the elective fraction against the reference inventories. If the
single-molecule constraint is what inflates the model's condensed linkages, the
ELECTIVE fraction -- not the total -- is the one that should match the references.

Usage:
    python code/weld_decomposition.py [--reps 24] [--n 120]
Writes results/weld_decomposition.csv and prints a report.
"""
import argparse
import csv
import json
import os
import sys

import mini_forward_generator as mfg  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import rel

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")

# References as % of the six modelled linkages (Table S3 / main-text Table 1).
REF = {
    "softwood": {"55": 6.0, "4o5": 1.5},
    "poplar":   {"55": 3.8, "4o5": 1.9},
    "high-S":   {"55": 1.1, "4o5": 0.6},
}
WOODS = [(0.0, "softwood"), (1.0, "poplar"), (3.0, "high-S")]


def load_params():
    with open(os.path.join(RESULTS, "fitted_params.json")) as fh:
        p = json.load(fh)
    return {k: v for k, v in p.items() if not k.startswith("_")}


def run_one(sg, seed, params, n):
    """Run the generator once, classifying every condensed bond as weld or elective."""
    g = mfg.MiniForwardKMC(sg, n_target=n, params=dict(params), seed=seed)
    rec = {("55", "weld"): 0, ("55", "elect"): 0,
           ("4o5", "weld"): 0, ("4o5", "elect"): 0}

    orig_fire, orig_union = g._fire_coupling, g.dsu.union
    state = {"sizes": None}

    def union_with_sizes(a, b):
        ra, rb = g.dsu.find(a), g.dsu.find(b)
        sizes = {}
        for uid in g.units:
            r = g.dsu.find(uid)
            sizes[r] = sizes.get(r, 0) + 1
        state["sizes"] = (sizes.get(ra, 1), sizes.get(rb, 1))
        return orig_union(a, b)

    def fire_and_classify(kind, *a, **k):
        before = {x: g.bonds.get(x, 0) for x in ("55", "4o5")}
        state["sizes"] = None
        out = orig_fire(kind, *a, **k)
        for x in ("55", "4o5"):
            if g.bonds.get(x, 0) > before[x] and state["sizes"]:
                sa, sb = state["sizes"]
                rec[(x, "weld" if (sa > 1 and sb > 1) else "elect")] += 1
        return out

    g.dsu.union = union_with_sizes
    g._fire_coupling = fire_and_classify
    g.run()
    return rec, g.summary()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=24)
    ap.add_argument("--n", type=int, default=120)
    args = ap.parse_args()

    params = load_params()
    n_bonds = args.n - 1
    required = round(params["n_nuclei"]) - 1

    print("weld decomposition of the condensed linkages")
    print("  n_target = %d  (%d bonds)   reps = %d   n_nuclei = %.3f -> %d"
          % (args.n, n_bonds, args.reps, params["n_nuclei"],
             round(params["n_nuclei"])))
    print("  welds structurally required to reach one molecule = n_nuclei - 1 = %d"
          "  (%.2f%% of %d bonds)\n" % (required, required / n_bonds * 100, n_bonds))

    hdr = ("%-9s %-6s %8s %8s %9s %8s %9s %11s"
           % ("wood", "bond", "model %", "weld %", "elective %", "ref %",
              "model-ref", "elective-ref"))
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for sg, name in WOODS:
        agg = {}
        for seed in range(args.reps):
            rec, _ = run_one(sg, seed, params, args.n)
            for k, v in rec.items():
                agg[k] = agg.get(k, 0) + v
        for bond in ("55", "4o5"):
            weld = agg[(bond, "weld")] / args.reps
            elect = agg[(bond, "elect")] / args.reps
            model_pct = (weld + elect) / n_bonds * 100
            weld_pct = weld / n_bonds * 100
            elect_pct = elect / n_bonds * 100
            ref = REF[name][bond]
            print("%-9s %-6s %8.2f %8.2f %9.2f %8.2f %+9.2f %+11.2f"
                  % (name, bond, model_pct, weld_pct, elect_pct, ref,
                     model_pct - ref, elect_pct - ref))
            rows.append({"wood": name, "sg": sg, "bond": bond,
                         "model_pct": round(model_pct, 3),
                         "weld_pct": round(weld_pct, 3),
                         "elective_pct": round(elect_pct, 3),
                         "ref_pct": ref,
                         "model_minus_ref": round(model_pct - ref, 3),
                         "elective_minus_ref": round(elect_pct - ref, 3)})

    print("\nREADING. The welds are the same %d bonds at every wood, because the"
          % required)
    print("construction always closes to one molecule. Removing them is what the")
    print("last column shows: the condensed linkage the model forms BY CHOICE, against")
    print("the measured value. If that column is near zero while 'model-ref' is not,")
    print("the model's condensed excess is the single-molecule constraint, not the fit.")

    out = os.path.join(RESULTS, "weld_decomposition.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote %s" % rel(out))


if __name__ == "__main__":
    main()
