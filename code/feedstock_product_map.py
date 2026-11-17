#!/usr/bin/env python
"""
feedstock_product_map.py -- what monomer composition does to the products of reductive
depolymerization, from molecules built by the forward generator (adopted constants,
results/fitted_params_ref_v11f_sw.json). Educational companion to feedstock_yield_screen.py.

Panels
  a, b  composition map: S/G (x) and H fraction (y) on a grid; colour = predicted monomer
        ceiling (a) and degree of condensation (b). Real feedstocks of the reference set are
        placed on the map at their measured S/G and H.
  c     where the aromatic units of six representative feedstocks end up after every
        beta-O-4 is cleaved: G, S and H monomers; dimers named by the bond that holds them
        (resinol beta-beta, phenylcoumaran beta-5, spirodienone beta-1, biphenyl 5-5,
        diaryl ether 4-O-5, or a linkage the model does not build); larger fragments.

Assumptions (same as feedstock_yield_screen.py): 120-unit molecules; every beta-O-4 and
nothing else cleaved, no repolymerization; 67.5 % of bonds are of the six modelled types
(median of the 25 quantitative NMR inventories in the fit), the rest are not cleaved.

    conda run -n base python code/feedstock_product_map.py
    conda run -n base python code/feedstock_product_map.py --si   # redraw from the saved JSON
        without title and footnote -> manuscript/figures/fig_product_map_si.png
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

N, SEED, KEEP6 = 120, 92000, 0.675
PARAMS = {k: v for k, v in json.load(open(os.path.join(REPO, "results", "fitted_params_ref_v11f_sw.json"))).items()
          if not k.startswith("_")}
SG_GRID = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4])
H_GRID = np.array([0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
GRID_REPS, BAR_REPS = 60, 60

MAP_POINTS = [  # sample_id, label (placed on the composition map)
    ("sette2011_norwayspruce_mwl", "spruce"), ("wangdeuss2021_pine_mwl", "pine"),
    ("rosado2021_rice_straw_mwl", "rice straw†"), ("bauer2012_miscanthus_eol95_4h", "Miscanthus†"),
    ("rencoret2009_paulownia_mwl", "Paulownia"), ("wangdeuss2021_walnut_mwl", "walnut shell†"),
    ("balakshin2016_acacia", "acacia"), ("sette2011_beech_mwl", "beech"),
    ("balakshin2016_egrandis", "E. grandis"), ("balakshin2016_eglobulus", "E. globulus"),
    ("balakshin2016_birch", "birch"), ("peng2021_poplar_ck_del", "poplar"),
    ("ralph2012_poplar_c3h14_cel", "poplar C3H-down"),
]
BAR_FEEDS = [
    ("sette2011_norwayspruce_mwl", "Spruce\n(softwood)"),
    ("bauer2012_miscanthus_eol95_4h", "Miscanthus†\n(grass)"),
    ("peng2021_poplar_ck_del", "Poplar\n(hardwood)"),
    ("balakshin2016_birch", "Birch\n(hardwood)"),
    ("ralph2012_poplar_c3h14_cel", "Poplar C3H-down\n(31 % H)"),
    ("stewart2009_poplar_f5h_acml", "Poplar F5H-up\n(S/G 38)"),
]
DIMER_NAME = {"bb": "resinol dimer (β-β)", "b5": "phenylcoumaran dimer (β-5)",
              "b1": "spirodienone dimer (β-1)", "55": "biphenyl dimer (5-5)",
              "4o5": "diaryl-ether dimer (4-O-5)", "x": "dimer by an unmodelled linkage"}


class RecordingKMC(mfg.MiniForwardKMC):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.edges = []

    def _bond(self, kind, a, b):
        self.edges.append((kind, a.uid, b.uid))
        super()._bond(kind, a, b)


def fragments(ids, edges):
    parent = {u: u for u in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for _, a, b in edges:
        parent[find(a)] = find(b)
    groups = {}
    for u in ids:
        groups.setdefault(find(u), []).append(u)
    return list(groups.values())


def products(g, rng):
    """Fractions of aromatic units in each product class after beta-O-4 cleavage."""
    ids = list(g.units)
    n = len(ids)
    # a beta-O-4 survives only if it was in truth an unmodelled bond (prob 1 - KEEP6);
    # every other bond survives; surviving ex-beta-O-4 bonds are labelled "x"
    kept = []
    for k, a, b in g.edges:
        if k == "bo4":
            if rng.rand() > KEEP6:
                kept.append(("x", a, b))
        else:
            kept.append((k, a, b))
    out = {}
    for frag in fragments(ids, kept):
        if len(frag) == 1:
            key = "%s monomer" % g.units[frag[0]].mtype
        elif len(frag) == 2:
            s = set(frag)
            kind = [k for k, a, b in kept if a in s and b in s][0]
            key = DIMER_NAME[kind] if kind != "x" else "residue"
        else:
            key = "residue"
        out[key] = out.get(key, 0) + len(frag)
    c = g.bonds
    dc = 100 * (2 * c["55"] + c["4o5"] + c["b5"]) / n
    return {k: v / n for k, v in out.items()}, dc


def run(sg, h, reps, seed):
    p = dict(PARAMS)
    p["pct_h"] = h
    rng = np.random.RandomState(seed)
    acc, dcs = {}, []
    for r in range(reps):
        g = RecordingKMC(sg, n_target=N, params=p, seed=seed + r)
        g.run()
        fr, dc = products(g, rng)
        for k, v in fr.items():
            acc[k] = acc.get(k, 0) + v / reps
        dcs.append(dc)
    return acc, float(np.mean(dcs))


def mono(acc):
    return 100 * sum(v for k, v in acc.items() if k.endswith("monomer"))


def main():
    rows = {r["sample_id"]: r for r in csv.DictReader(open(os.path.join(REPO, "data", "reference_inventories.csv")))}
    comp = lambda sid: (float(rows[sid]["sg_ratio"]), float(rows[sid]["h_frac"] or 0))  # noqa: E731
    Y = np.zeros((len(H_GRID), len(SG_GRID)))
    D = np.zeros_like(Y)
    for i, h in enumerate(H_GRID):
        for j, sg in enumerate(SG_GRID):
            acc, dc = run(sg, h, GRID_REPS, SEED + 1000 * i + 37 * j)
            Y[i, j], D[i, j] = mono(acc), dc
    bars = []
    for sid, lab in BAR_FEEDS:
        sg, h = comp(sid)
        acc, dc = run(sg, h, BAR_REPS, SEED + 7)
        bars.append((lab, sg, h, acc, dc))
        print("%-28s S/G %5.2f H %3.0f%%  monomers %5.1f  DC %5.1f" % (lab.replace("\n", " "), sg, 100 * h, mono(acc), dc))
    with open(os.path.join(REPO, "results", "feedstock_product_map.json"), "w") as fh:
        json.dump({"sg_grid": SG_GRID.tolist(), "h_grid": H_GRID.tolist(), "monomer_ceiling": Y.tolist(),
                   "condensation": D.tolist(),
                   "bars": [{"label": l, "sg": s, "h": h, "products": a, "dc": d} for l, s, h, a, d in bars]},
                  fh, indent=1)
    plot(Y, D, bars, [(lab, *comp(sid)) for sid, lab in MAP_POINTS])


def plot(Y, D, bars, pts, si=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    INK, MUTED = "#1a1a1a", "#6b6a66"
    fig = plt.figure(figsize=(12, 12.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.05, 0.55], hspace=0.48, wspace=0.18,
                          left=0.07, right=0.98, top=0.97 if si else 0.95, bottom=0.04 if si else 0.08)
    for col, (Z, title, cmap, lab) in enumerate([
            (Y, "a  Monomers released by β-O-4 cleavage", "Blues", "monomer ceiling, % of aromatic units"),
            (D, "b  Condensation: bonds that resist cleavage", "Reds", "condensed ring carbons per 100 units")]):
        ax = fig.add_subplot(gs[0, col])
        im = ax.contourf(SG_GRID, 100 * H_GRID, Z, levels=12, cmap=cmap)
        cs = ax.contour(SG_GRID, 100 * H_GRID, Z, levels=6, colors="white", linewidths=0.6)
        ax.clabel(cs, fmt="%.0f", fontsize=7, colors="white")
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(lab, fontsize=8)
        for name, sg, h in pts:
            ax.scatter(min(sg, 4), 100 * h, s=18, facecolor="white", edgecolor=INK, lw=0.8, zorder=5)
            right = min(sg, 4) > 3.2
            ax.annotate(name, (min(sg, 4), 100 * h), xytext=(-4 if right else 4, 3), textcoords="offset points",
                        fontsize=6.5, color=INK, zorder=6, ha="right" if right else "left")
        ax.set_xlabel("S/G  (more syringyl →)")
        ax.set_ylabel("H units, %")
        ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    # c: product slate
    ax = fig.add_subplot(gs[1, :])
    order = ["G monomer", "S monomer", "H monomer", "resinol dimer (β-β)", "phenylcoumaran dimer (β-5)",
             "spirodienone dimer (β-1)", "biphenyl dimer (5-5)", "diaryl-ether dimer (4-O-5)",
             "residue"]
    legend = {"G monomer": "G monomers (4-propylguaiacol type)", "S monomer": "S monomers (4-propylsyringol type)",
              "H monomer": "H monomers (4-propylphenol type)",
              "residue": "residue: oligomers, and fragments held by\nlinkages other than the six (not cleaved)"}
    colors = ["#1f5fa8", "#5b9bd5", "#a9cdee", "#2e8b57", "#7fbf7f", "#c4e3b5", "#d55e00", "#f0a35e",
              "#8a8a8a"]
    x = np.arange(len(bars))
    bottom = np.zeros(len(bars))
    for k, c in zip(order, colors):
        v = np.array([100 * b[3].get(k, 0) for b in bars])
        ax.bar(x, v, 0.62, bottom=bottom, color=c, label=legend.get(k, k), zorder=3, edgecolor="white", lw=0.4)
        for i, (vi, bi) in enumerate(zip(v, bottom)):
            if vi >= 4:
                ax.text(i, bi + vi / 2, "%.0f" % vi, ha="center", va="center", fontsize=7,
                        color="white" if c in ("#1f5fa8", "#2e8b57", "#d55e00", "#8a8a8a") else INK, zorder=4)
        bottom += v
    ax.set_xticks(x)
    ax.set_xticklabels(["%s\nS/G %.3g · H %.0f%%" % (b[0], b[1], 100 * b[2]) for b in bars], fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of the lignin's aromatic units")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=7.5)
    ax.set_title("c  Where each aromatic unit ends up after reductive depolymerization",
                 loc="left", fontsize=10, fontweight="bold", color=INK)
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width * 0.74, pos.height])
    # d: application streams (alternative uses of the same product slate)
    axd = fig.add_subplot(gs[2, :])
    mono_ = np.array([100 * sum(v for k, v in b[3].items() if k.endswith("monomer")) for b in bars])
    dim_ = np.array([100 * sum(v for k, v in b[3].items() if "dimer" in k) for b in bars])
    res_ = np.array([100 * b[3].get("residue", 0) for b in bars])
    M = np.vstack([mono_, mono_ + dim_, res_])
    rows_lab = ["phenolic chemicals and polymer\nbuilding blocks  ← monomers",
                "aromatic oil for fuels after\nhydrodeoxygenation  ← monomers + dimers",
                "lignin materials (resins, polyols,\ncarbon) or energy  ← residue"]
    im = axd.imshow(M, aspect="auto", cmap="YlGn", vmin=0, vmax=80)
    for r in range(M.shape[0]):
        for c in range(M.shape[1]):
            axd.text(c, r, "%.0f %%" % M[r, c], ha="center", va="center", fontsize=8.5,
                     color="white" if M[r, c] > 50 else INK)
    axd.set_yticks(range(3))
    axd.set_yticklabels(rows_lab, fontsize=8)
    axd.yaxis.tick_right()
    axd.set_xticks(range(len(bars)))
    axd.set_xticklabels([b[0].split("\n")[0] for b in bars], fontsize=8)
    axd.tick_params(length=0)
    for sp in axd.spines.values():
        sp.set_visible(False)
    axd.set_title("d  Possible uses: share of each lignin's aromatic units available to each stream "
                  "(alternatives, not additive)", loc="left", fontsize=10, fontweight="bold", color=INK)
    pd_ = axd.get_position()
    axd.set_position([pos.x0, pd_.y0, pos.width * 0.74, pd_.height])
    ax.set_xlim(-0.5, len(bars) - 0.5)
    if si:
        fig.savefig(os.path.join(REPO, "manuscript", "figures", "fig_product_map_si.png"), dpi=300)
        return
    fig.suptitle("Monomer composition decides the products: predictions of the calibrated generator",
                 x=0.07, ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.text(0.07, 0.005,
             "120-unit molecules at the adopted constants (60 per grid point and per bar). Every β-O-4 cleaved, nothing else, "
             "no repolymerization.\n67.5 % of bonds are of the six modelled types (median of 25 quantitative NMR inventories); "
             "the rest are not cleaved.\n† grass / shell: γ-acylation and tricin not modelled. Points at S/G > 4 are drawn at the edge.",
             fontsize=7, color=MUTED, va="bottom")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(REPO, "results", "figures", "feedstock_product_map." + ext), dpi=300)


def replot_si():
    d = json.load(open(os.path.join(REPO, "results", "feedstock_product_map.json")))
    rows = {r["sample_id"]: r for r in csv.DictReader(open(os.path.join(REPO, "data", "reference_inventories.csv")))}
    pts = [(lab, float(rows[sid]["sg_ratio"]), float(rows[sid]["h_frac"] or 0)) for sid, lab in MAP_POINTS]
    bars = [(b["label"], b["sg"], b["h"], b["products"], b["dc"]) for b in d["bars"]]
    plot(np.array(d["monomer_ceiling"]), np.array(d["condensation"]), bars, pts, si=True)


if __name__ == "__main__":
    if "--si" in sys.argv:
        replot_si()
    else:
        main()
