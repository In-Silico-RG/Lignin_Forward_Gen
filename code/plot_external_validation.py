#!/usr/bin/env python
"""
plot_external_validation.py -- reviewers' Major #1 ("n=3 / across feedstocks")

Independent-literature check on beta-O-4(S/G). We do NOT plot the model's normalised
(% of six linkages) output against mol% literature -- that would collide with the
connectivity-gap denominator (the single molecule over-bonds, so its beta-O-4 per 100
units ~67 is not comparable to a measured mol% ~46). The honest, denominator-clean
comparison is: do the manuscript's REFERENCE beta-O-4 anchors (per 100 C9 ~ mol%, which
the model reproduces to z<1) sit within the independent multi-feedstock scatter?

Points from data/independent_feedstock_linkages.csv (bo4_basis == mol%):
  Bourmaud 2024 (whole-cell-wall HSQC0) and Laskar 2013 (isolated MWL). S/G values are
  type-typical placeholders (softwood 0, poplar ~1, birch ~2.3, eucalyptus ~2.7,
  high-S ~3, grasses ~0.8) -- replace with each paper's own HSQC S/G to finalise.

The figure makes three honest points at once:
  1. Across independent woods, beta-O-4 rises with S/G, and the paper's anchors sit
     centrally in that scatter.
  2. The field's own numbers scatter ~15 points by METHOD (e.g. birch 60 MWL vs 72
     HSQC0) -- context for "your n=3": no three lignins pin beta-O-4 tightly.
  3. Grasses are method-split outliers (WCW-HSQC0 ~30 vs MWL ~75); they cannot arbitrate
     the model without resolving the method gap.

Run:    conda run -n base python code/plot_external_validation.py
Output: manuscript/submission/figures/fig_external_validation.pdf (+ .png)
"""
import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
CSV = os.path.join(ROOT, "data", "independent_feedstock_linkages.csv")
OUT = os.path.join(ROOT, "manuscript", "submission", "figures", "fig_external_validation")

BLUE, GREEN, ORANGE, INK, MUTED, GRID = (
    "#2a78d6", "#008300", "#eb6834", "#0b0b0b", "#898781", "#e1e0d9")

# manuscript reference anchors (RAW per-100-C9 beta-O-4, which the model reproduces)
REF = [(0.0, 46.0), (1.0, 62.0), (3.0, 78.0)]

METHOD_STYLE = {           # marker + colour by measurement method
    "WCW-HSQC0": ("o", BLUE, "Bourmaud 2024 (whole cell wall, HSQC0)"),
    "MWL":       ("s", GREEN, "Laskar 2013 (isolated MWL)"),
}


def load():
    rows = []
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            if r["bo4_basis"] != "mol%":
                continue
            rows.append({
                "sg": float(r["sg_ratio"]), "bo4": float(r["bo4"]),
                "err": float(r["bo4_err"]) if r["bo4_err"] else None,
                "method": r["method"], "type": r["type"], "name": r["feedstock"],
            })
    return rows


def main():
    rows = load()
    fig, ax = plt.subplots(figsize=(5.2, 3.6))

    # independent points, styled by method; grasses ringed to flag them as outliers
    for r in rows:
        mk, col, _ = METHOD_STYLE[r["method"]]
        is_grass = r["type"] == "grass"
        ax.errorbar(r["sg"], r["bo4"], yerr=r["err"], fmt=mk, color=col,
                    ms=8 if not is_grass else 8, mfc=col if not is_grass else "white",
                    mec=col, mew=1.4, ecolor=col, elinewidth=0.8, capsize=2.5, zorder=3)
        dy = 2.4 if r["bo4"] > 45 else -3.2
        ax.annotate(r["name"], (r["sg"], r["bo4"]), xytext=(4, dy),
                    textcoords="offset points", fontsize=6.2, color=MUTED)

    # manuscript reference anchors + connecting trend
    rx = [p[0] for p in REF]; ry = [p[1] for p in REF]
    ax.plot(rx, ry, "-", color=ORANGE, lw=1.6, zorder=4)
    ax.plot(rx, ry, "D", color=ORANGE, ms=9, mec="white", mew=1.2, zorder=5)

    ax.set_xlabel("S/G ratio  (type-typical; verify per source)")
    ax.set_ylabel("β-O-4  (mol% / per 100 aromatic units)")
    ax.set_xlim(-0.3, 3.4)
    ax.set_ylim(20, 90)
    ax.set_title("β-O-4(S/G): the paper's anchors sit within the independent scatter",
                 fontsize=9.5, loc="left", color=INK)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Line2D([0], [0], marker="D", color=ORANGE, lw=1.6, ms=9, mec="white",
               label="this work, reference anchors (model reproduces, z<1)"),
        Line2D([0], [0], marker="o", ls="none", color=BLUE, ms=8,
               label=METHOD_STYLE["WCW-HSQC0"][2]),
        Line2D([0], [0], marker="s", ls="none", color=GREEN, ms=8,
               label=METHOD_STYLE["MWL"][2]),
        Line2D([0], [0], marker="o", ls="none", color=MUTED, ms=8, mfc="white",
               mec=MUTED, label="grass (method-split outlier; hollow)"),
    ]
    ax.legend(handles=handles, fontsize=6.6, frameon=False, loc="upper left",
              bbox_to_anchor=(0.0, 0.98))

    fig.tight_layout()
    fig.savefig(OUT + ".pdf", bbox_inches="tight")
    fig.savefig(OUT + ".png", dpi=200, bbox_inches="tight")
    print("wrote", os.path.relpath(OUT + ".pdf", ROOT))


if __name__ == "__main__":
    main()
