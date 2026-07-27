#!/usr/bin/env python
"""
plot_intermediate_sg.py  --  figure for the intermediate-S/G experiment

Reads results/intermediate_sg_experiment.csv (written by intermediate_sg_experiment.py)
and draws the predicted linkage inventory vs S/G, with:
  * the three CALIBRATION anchors (S/G 0, 1, 3) as filled markers,
  * the interpolated predictions as the smooth curve between them,
  * the S/G > 3 EXTRAPOLATION region shaded (mechanism-guided but unvalidated),
  * +/- 1 std bands over seeds.

Two panels because beta-O-4 (68-91%) and the five minor linkages (<17%) do not
share a sensible y-scale.

Only matplotlib + numpy needed -- run in an env that has them (the `lignin` env's
matplotlib is broken; `base` works):
    conda run -n base python code/plot_intermediate_sg.py
Output: manuscript/submission/figures/fig_intermediate_sg.pdf (+ .png)
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
CSV = os.path.join(ROOT, "results", "intermediate_sg_experiment.csv")
OUT = os.path.join(ROOT, "manuscript", "submission", "figures", "fig_intermediate_sg")

KEYS = ["bo4", "b5", "bb", "b1", "55", "4o5"]
LABELS = {"bo4": r"$\beta$-O-4", "b5": r"$\beta$-5", "bb": r"$\beta$-$\beta$",
          "b1": r"$\beta$-1", "55": "5-5", "4o5": "4-O-5"}
# Okabe-Ito colourblind-safe palette
COLOR = {"bo4": "#0072B2", "b5": "#D55E00", "bb": "#009E73",
         "b1": "#CC79A7", "55": "#E69F00", "4o5": "#56B4E9"}
MINOR = ["b5", "bb", "b1", "55", "4o5"]
EXTRAP_START = 3.0


def load():
    rows = []
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    sg = np.array([float(r["sg_ratio"]) for r in rows])
    role = [r["role"] for r in rows]
    mean = {k: np.array([100 * float(r[k]) for r in rows]) for k in KEYS}
    std = {k: np.array([100 * float(r[k + "_std"]) for r in rows]) for k in KEYS}
    return sg, role, mean, std


def draw_panel(ax, sg, role, mean, std, keys, label_dy=None):
    xmax = 6.0  # calibrated window [0,3] + extrapolation tail (S/G 4, 6); 10/49 off-panel
    label_dy = label_dy or {}
    anch = np.array([r == "anchor" for r in role])
    m = sg <= xmax
    ax.axvspan(EXTRAP_START, xmax + 0.9, color="0.92", zorder=0)
    for k in keys:
        ax.fill_between(sg[m], (mean[k] - std[k])[m], (mean[k] + std[k])[m],
                        color=COLOR[k], alpha=0.15, lw=0)
        ax.plot(sg[m], mean[k][m], "-", color=COLOR[k], lw=1.8, zorder=3)
        # anchors as filled diamonds, predictions implicit in the line
        am = m & anch
        ax.plot(sg[am], mean[k][am], "D", color=COLOR[k], ms=7,
                mec="white", mew=1.0, zorder=4)
        # label at the right edge, with optional manual de-collision offset
        xr = sg[m][-1]
        ax.annotate(LABELS[k], (xr, mean[k][m][-1]), xytext=(5, label_dy.get(k, 0)),
                    textcoords="offset points", va="center", fontsize=9,
                    color=COLOR[k])
    ax.set_xlim(-0.1, xmax + 0.9)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    sg, role, mean, std = load()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.35]})

    draw_panel(ax1, sg, role, mean, std, ["bo4"])
    ax1.set_ylabel("% of modelled linkages")
    ax1.set_title("Predicted linkage inventory vs monomer composition (S/G)",
                  fontsize=11, loc="left")

    # nudge the crowded minor-panel labels apart (they converge near 2-4%)
    # nudge the crowded minor-panel labels apart (at S/G 6 they stack 4.0/3.1/2.6/2.6/1.4)
    draw_panel(ax2, sg, role, mean, std, MINOR,
               label_dy={"bb": 6, "4o5": 2, "b1": -3, "b5": -12, "55": -6})
    ax2.set_ylabel("% of modelled linkages")
    ax2.set_xlabel("S/G ratio")
    ax2.set_ylim(0, None)

    # shared legend explaining the encoding
    handles = [
        Line2D([0], [0], marker="D", color="0.3", ls="none", ms=7, mec="white",
               label="calibration anchor (S/G 0, 1, 3)"),
        Line2D([0], [0], color="0.3", lw=1.8, label="model prediction (interpolation)"),
        Patch(facecolor="0.92", label="extrapolation (S/G > 3, unvalidated)"),
    ]
    ax1.legend(handles=handles, fontsize=8, frameon=False, loc="center right")

    fig.tight_layout()
    fig.savefig(OUT + ".pdf", bbox_inches="tight")
    fig.savefig(OUT + ".png", dpi=200, bbox_inches="tight")
    print("wrote", os.path.relpath(OUT + ".pdf", ROOT))
    print("wrote", os.path.relpath(OUT + ".png", ROOT))


if __name__ == "__main__":
    main()
