#!/usr/bin/env python
"""
2026-09-14: the authors' unit defaults and the provisional vetting, applied to
data/reference_inventories.csv. RUN ONCE (commits bd7f381); kept for provenance.

Defaults (Aldo, chat 2026-09-14, "safe defaults" accepted):
  * Capanema 2005 E. grandis 5-5 printed 3 per 100 Ar with units unstated -> treated as
    C9 units in biphenyl -> 1.5 bonds (consistent with the paper's resinol convention).
  * Peng 2021 beta-beta printed 13.02 / 10.76 per 100 Ar -> units in resinol -> 6.51 / 5.38 bonds.
  * The three legacy anchors get distinct `source` strings so source weighting treats
    them as three sources.
Provisional vetting (Aldo delegated, "take over"):
  * wood samples with quantitative or relative HSQC / 13C -> use_fit = 1
  * grasses, bark, cork, phloem, nut shell -> use_fit = 0 (validation-only), because
    gamma-acylation, tricin and stilbene chemistry are not in the model
  * every touched row: vetted = provisional-claude
"""
import csv
import os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "reference_inventories.csv")
VAL_ONLY = {"lourenco2016_qsuber_phloem_ml", "lourenco2016_qsuber_cork_ml", "rosado2021_rice_husk_mwl",
            "rosado2021_rice_straw_mwl", "wangdeuss2021_walnut_mwl", "wangdeuss2021_reed_mwl",
            "delrio2012_elephantgrass_cortex_mwl", "delrio2012_elephantgrass_pith_mwl"}


def apply(rows):
    for r in rows:
        if r["sample_id"] == "capanema2005_egrandis":
            r["c55"] = "1.5"
            r["notes"] = r["notes"].replace(
                "5-5 printed 3 -- UNITS OR BONDS NOT STATED, recorded as printed (FLAG)",
                "5-5 printed 3, taken as C9 units in biphenyl -> 1.5 bonds (authors' default 2026-09-14, consistent with the resinol convention)")
        if r["sample_id"].startswith("peng2021"):
            v = float(r["bb"]); r["bb"] = "%.2f" % (v / 2)
            r["notes"] = r["notes"].replace(
                "bb likely units-in-resinol convention -- CHECK (13.02 per 100 Ar is high for bonds)",
                "bb printed %.2f taken as units in resinol -> %.2f bonds (authors' default 2026-09-14)" % (v, v / 2)
            ).replace("as above", "bb printed %.2f -> %.2f bonds (units in resinol, authors' default)" % (v, v / 2))
        if r["sample_id"].startswith("anchor_"):
            r["source"] = "legacy anchor %s (calibrate_full.RAW_REFERENCE)" % r["sample_id"].split("_", 1)[1]
        if r["vetted"] == "pending":
            r["vetted"] = "provisional-claude"
            r["use_fit"] = "0" if r["sample_id"] in VAL_ONLY else "1"
    return rows


if __name__ == "__main__":
    raise SystemExit("Already applied on 2026-09-14 (commit bd7f381). Kept for provenance; do not re-run.")
    rows = list(csv.DictReader(open(P)))
    hdr = list(rows[0].keys())
    with open(P, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr); w.writeheader(); w.writerows(apply(rows))
