#!/usr/bin/env python
"""
feedstock_yield_screen.py -- theoretical aromatic yields of real feedstocks, computed on
molecules built by the forward generator at the adopted constants (v1.1.1, _v11e source-
weighted fit, results/fitted_params_ref_v11f_sw.json).

For each feedstock (S/G and H fraction from data/reference_inventories.csv) the generator
builds REPS molecules of N = 120 units. Every bond is recorded with its two partners
(donor first, as _fire_coupling / _nucleate_one pass them to _bond), and the molecule is
then "processed" as a graph:

  ether cleavage (reductive catalytic fractionation, hydrogenolysis, acidolysis):
      cut every beta-O-4; a unit is a MONOMER if it ends in a fragment of one unit,
      part of a DIMER if in a fragment of two. Theoretical ceiling (100 % cleavage of
      beta-O-4, no other bond cleaved, no repolymerization).
  side-chain oxidation (nitrobenzene / alkaline oxidation, aldehydes):
      every side chain is cleaved; a unit's ring is released as a monomer (vanillin,
      syringaldehyde, p-hydroxybenzaldehyde) unless the ring itself is bonded to another
      unit: acceptor C5 of beta-5, acceptor C1 of beta-1, both rings of 5-5, both rings of
      4-O-5 (diaryl ether). Same idealization.
  condensation: degree of condensation as in the paper, substituted aromatic carbons per
      100 units, (2 n_5-5 + n_4-O-5 + n_beta-5) / units x 100; and the share of rings
      that are condensed (the complement of the oxidation yield).

Outputs: results/feedstock_yield_screen.csv, results/figures/feedstock_yield_screen.{png,pdf}

    conda run -n base python code/feedstock_yield_screen.py
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mini_forward_generator as mfg  # noqa: E402

N, REPS, SEED = 120, 30, 91000
# Share of a lignin's inter-unit bonds that are of the six modelled types: median of the
# six-linkage sums of the 25 quantitative per-100-aromatic-unit NMR inventories in the fit
# (13C-quant + HSQC-quant, anchors and semi-quantitative HSQC excluded): 67.5 (range 59.8-84.5).
# The rest are linkages the model does not build; ether cleavage does not cut them.
KEEP6 = 0.675
PARAMS = json.load(open(os.path.join(REPO, "results", "fitted_params_ref_v11f_sw.json")))
PARAMS = {k: v for k, v in PARAMS.items() if not k.startswith("_")}

# (sample_id in the reference set, short label, group). Anchors excluded (no single source).
FEEDSTOCKS = [
    ("sette2011_norwayspruce_mwl", "Norway spruce", "softwood"),
    ("wangdeuss2021_pine_mwl", "Pine", "softwood"),
    ("rosado2021_rice_straw_mwl", "Rice straw", "grass"),
    ("bauer2012_miscanthus_eol95_4h", "Miscanthus", "grass"),
    ("rencoret2009_paulownia_mwl", "Paulownia", "hardwood"),
    ("rencoret2019_olive_mwl", "Olive pruning", "hardwood"),
    ("wangdeuss2021_walnut_mwl", "Walnut shell", "other"),
    ("balakshin2016_acacia", "Acacia mangium", "hardwood"),
    ("balakshin2016_redmaple", "Red maple", "hardwood"),
    ("delrio2012_elephantgrass_cortex_mwl", "Elephant grass", "grass"),
    ("sette2011_beech_mwl", "Beech", "hardwood"),
    ("peng2021_poplar_ck_del", "Poplar", "hardwood"),
    ("balakshin2016_egrandis", "E. grandis", "hardwood"),
    ("balakshin2016_redoak", "Red oak", "hardwood"),
    ("balakshin2016_eglobulus", "E. globulus", "hardwood"),
    ("balakshin2016_birch", "Birch", "hardwood"),
    ("ralph2012_poplar_c3h14_cel", "Poplar C3H-down", "transgenic"),
    ("stewart2009_poplar_f5h_acml", "Poplar F5H-up", "transgenic"),
]


class RecordingKMC(mfg.MiniForwardKMC):
    """The generator unchanged, plus a list of (kind, donor uid, acceptor uid)."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.edges = []

    def _bond(self, kind, a, b):
        self.edges.append((kind, a.uid, b.uid))
        super()._bond(kind, a, b)


def components(n_ids, edges):
    parent = {u: u for u in n_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    size = {}
    for u in n_ids:
        r = find(u)
        size[r] = size.get(r, 0) + 1
    return {u: size[find(u)] for u in n_ids}


def process(g, keep6=None, rng=None):
    """keep6: fraction of bonds that are of the six modelled types in the real lignin (the
    sample's measured six-linkage sum per 100 aromatic units / 100). If given, each bond is
    turned, with probability 1 - keep6, into an unmodelled bond that ether cleavage does not
    cut; returned as a 5th value (monomers after that substitution)."""
    ids = list(g.units)
    n = len(ids)
    kept = [(a, b) for k, a, b in g.edges if k != "bo4"]
    comp = components(ids, kept)
    mono = sum(1 for u in ids if comp[u] == 1) / n
    di = sum(1 for u in ids if comp[u] == 2) / n
    ring = set()
    for k, a, b in g.edges:
        if k in ("b5", "b1"):
            ring.add(b)
        elif k in ("55", "4o5"):
            ring.update((a, b))
    ox = 1 - len(ring) / n
    c = g.bonds
    dc = 100 * (2 * c["55"] + c["4o5"] + c["b5"]) / n
    mono_x = np.nan
    if keep6 is not None:
        kept_x = [(a, b) for k, a, b in g.edges if k != "bo4" or rng.rand() > keep6]
        comp_x = components(ids, kept_x)
        mono_x = sum(1 for u in ids if comp_x[u] == 1) / n
        di_x = sum(1 for u in ids if comp_x[u] == 2) / n
        return mono, di, ox, dc, mono_x, di_x
    return mono, di, ox, dc, mono_x, np.nan


def main():
    rows = {r["sample_id"]: r for r in csv.DictReader(open(os.path.join(REPO, "data", "reference_inventories.csv")))}
    out = []
    for sid, label, group in FEEDSTOCKS:
        r = rows[sid]
        sg = float(r["sg_ratio"])
        h = float(r["h_frac"]) if r["h_frac"].strip() else 0.0
        p = dict(PARAMS)
        p["pct_h"] = h
        six = [r[k] for k in ("bo4", "b5", "bb", "b1", "c55", "c4o5")]
        absolute = r["denominator"] == "per100Ar"
        keep6 = KEEP6   # same for every feedstock: the predictor gets S/G and H only
        nmr = (float(r["bo4"]) / 100) ** 2 * 100 if absolute else np.nan
        rng = np.random.RandomState(SEED)
        res = []
        for rep in range(REPS):
            g = RecordingKMC(sg, n_target=N, params=p, seed=SEED + rep)
            g.run()
            res.append(process(g, keep6, rng))
        a = np.array(res)
        m, s = a.mean(0), a.std(0)
        out.append(dict(sample_id=sid, label=label, group=group, sg=sg, h=h,
                        rcf_mono=100 * m[0], rcf_mono_sd=100 * s[0],
                        rcf_mono_di=100 * (m[0] + m[1]), rcf_mono_di_sd=100 * (a[:, 0] + a[:, 1]).std(),
                        ox_mono=100 * m[2], ox_mono_sd=100 * s[2], dc=m[3], dc_sd=s[3],
                        six_sum_measured=sum(float(v) for v in six if v.strip()) if absolute else np.nan,
                        rcf_mono_pred=100 * m[4], rcf_mono_pred_sd=100 * s[4],
                        rcf_mono_di_pred=100 * (m[4] + m[5]), rcf_mono_di_pred_sd=100 * (a[:, 4] + a[:, 5]).std(),
                        bo4_measured=float(r["bo4"]) if absolute else np.nan,
                        rcf_ceiling_nmr=nmr, nmr_method=r["method"] if absolute else ""))
        print("%-16s S/G %5.2f H %3.0f%%  RCF mono %5.1f  pred %5.1f  NMR %5.1f  oxid %5.1f  DC %5.1f"
              % (label, sg, 100 * h, out[-1]["rcf_mono"], out[-1]["rcf_mono_pred"],
                 nmr, out[-1]["ox_mono"], out[-1]["dc"]))
    with open(os.path.join(REPO, "results", "feedstock_yield_screen.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    plot(out)


def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    out = sorted(out, key=lambda o: (o["group"] == "transgenic", o["sg"]))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    INK, MUTED = "#1a1a1a", "#6b6a66"
    C_MONO, C_DI, C_DC = "#2a78d6", "#9cc3ef", "#b3261e"
    GROUP_BG = {"softwood": "#eef3ea", "grass": "#fbf3e2", "other": "#f1eff7",
                "hardwood": "#ffffff", "transgenic": "#f6ecec"}
    x = np.arange(len(out))
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 7.4), sharex=True,
                                  gridspec_kw={"height_ratios": [1.6, 1], "hspace": 0.12})
    for i, o in enumerate(out):
        for a_ in (ax, ax2):
            a_.axvspan(i - 0.5, i + 0.5, color=GROUP_BG[o["group"]], zorder=0, lw=0)
    mono = [o["rcf_mono_pred"] for o in out]
    di = [o["rcf_mono_di_pred"] - o["rcf_mono_pred"] for o in out]
    ax.bar(x, mono, 0.6, color=C_MONO, zorder=3, yerr=[o["rcf_mono_pred_sd"] for o in out],
           error_kw=dict(lw=0.6, capsize=1.5, ecolor=MUTED), label="monomers")
    ax.bar(x, di, 0.6, bottom=mono, color=C_DI, zorder=3, label="dimers")
    ax.legend(loc="upper left", frameon=False, fontsize=8, ncol=2)
    ax.set_ylabel("reductive depolymerization (RCF)\npredicted yield ceiling,\n% of aromatic units")
    ax.set_ylim(0, 75)
    ax.set_title("What the generator predicts from monomer composition alone (S/G and H)",
                 fontsize=11, loc="left", color=INK, fontweight="bold")
    ax2.bar(x, [o["dc"] for o in out], 0.6, color=C_DC, yerr=[o["dc_sd"] for o in out],
            error_kw=dict(lw=0.6, capsize=1.5, ecolor=MUTED), zorder=3)
    ax2.set_ylabel("processing challenge:\ncondensed ring carbons\nper 100 units")
    for a_ in (ax, ax2):
        a_.spines[["top", "right"]].set_visible(False)
        a_.grid(axis="y", color="#e1e0d9", lw=0.6, zorder=1)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["%s%s\nS/G %s · H %s%%" % (o["label"], "†" if o["group"] in ("grass", "other") else "",
                                                    ("%.2g" % o["sg"]), ("%.0f" % (100 * o["h"]))) for o in out],
                        rotation=55, ha="right", fontsize=7.5)
    fig.text(0.01, 0.005,
             "Generator at the adopted constants: 30 molecules of 120 units per feedstock, input S/G and H only. Yield ceiling: every "
             "β-O-4 cleaved, nothing else, no repolymerization;\n"
             "two-thirds of the bonds are of the six modelled linkage types (67.5 %, median of 25 quantitative NMR inventories), "
             "the rest are not cleaved by RCF.\n"
             "Background: green softwood, white hardwood, beige grass, lilac nut shell, red transgenic poplar. "
             "† not in the calibration (γ-acylation, tricin not modelled).",
             fontsize=7, color=MUTED, va="bottom")
    fig.subplots_adjust(left=0.10, right=0.99, top=0.94, bottom=0.25)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(REPO, "results", "figures", "feedstock_yield_screen." + ext), dpi=300)


if __name__ == "__main__":
    main()
