#!/usr/bin/env python
"""
2026-09-15: halve the beta-beta values of the HSQC sources that report the resinol
Calpha-Halpha volume without halving it, i.e. count C9 UNITS in resinol (two per bond),
the same convention Balakshin 2016 prints explicitly and we already halve.

Evidence (docs/15 section 5, 2026-09-15):
  * Rosado 2023 (del Rio group) Methods: "the integration volumes were divided by two"
    for S2,6 and H2,6 -- stated for the aromatics, NOT stated for resinol, which is
    "estimated from the Calpha/Halpha correlations".
  * Wang & Deuss 2021, Table 1 footnote a: molar % = alpha-signal integral / total of
    S (1/2 S2,6) + S' + G + G' + H (1/2 H2,6): aromatics halved, resinol not.
  * Lourenco 2016, Rencoret 2019: "from the Calpha-Halpha correlations", nothing else.
  * Happs 2021: refers to the Mansfield 2012 protocol, no detail; its beta-beta/beta-O-4
    (0.18-0.40) is 3-5x the 13C bonds ratio for poplar (0.05-0.09).
  * Rencoret 2011 E. globulus resinol 13-14 per 100 Ar vs Balakshin 2016 E. globulus
    10 UNITS (5 bonds); Wang & Deuss birch 11 vs Balakshin birch 9 units.
  * Ralph-group poplars (Stewart 2009, Ralph 2012, Dou 2018) already sit at 0.07-0.11,
    the bonds level, so they are NOT halved. Sette QQ-HSQC compares its resinol count
    directly with Capanema's 13C structures (Table 1), so it counts bonds: NOT halved.
    Peng 2021 was halved on 2026-09-14.

Delegated decision (Aldo, 2026-09-14: "take over"); reversible by editing `bb` and the
note. RUN ONCE; kept for provenance.
"""
import csv
import os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "reference_inventories.csv")
HALVE_PREFIXES = ("rencoret2011_", "rencoret2019_", "rencoret2009_", "rosado2023_", "lourenco2016_",
                  "wangdeuss2021_", "happs2021_", "rosado2021_", "delrio2012_")
NOTE = "; bb printed %s taken as C9 units in resinol (Calpha-Halpha volume, not halved by the source) -> %s bonds (2026-09-15 rule, see data/staging)"


def apply(rows):
    n = 0
    for r in rows:
        if r["sample_id"].startswith(HALVE_PREFIXES) and r["bb"]:
            v = float(r["bb"])
            r["bb"] = ("%.2f" % (v / 2)).rstrip("0").rstrip(".")
            r["notes"] = r["notes"] + NOTE % (("%g" % v), r["bb"])
            n += 1
    return rows, n


if __name__ == "__main__":
    rows = list(csv.DictReader(open(P)))
    hdr = list(rows[0].keys())
    if any("2026-09-15 rule" in r["notes"] for r in rows):
        raise SystemExit("Already applied; refusing to halve twice.")
    rows, n = apply(rows)
    with open(P, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr); w.writeheader(); w.writerows(rows)
    print("halved beta-beta in %d rows" % n)
