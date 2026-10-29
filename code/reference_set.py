#!/usr/bin/env python
"""
reference_set.py

The expanded reference set (data/reference_inventories.csv) and the objective
that compares a simulated lignin with it. Design: docs/14_reference_set_design.md.

One row per measured sample. A sample reports whichever linkages its method
resolved; blanks mean "not measured", never zero. Every observable is a RATIO TO
beta-O-4 (r_k = x_k / x_bo4), which is independent of the sample's denominator
(per 100 Ar, % of side chains, % of the six) and of the model's over-bonding, and
is defined for partial inventories. The degree of condensation (DC, quantitative
13C) is mapped to the model as the number of substituted aromatic carbons per unit:
    DC_model = (2 n_55 + n_4O5 + n_b5) / n_units * 100
and enters as DC / bo4 like everything else.

Usage from a calibration script:
    refs = load_reference_set()                 # rows with use_fit == 1
    obs  = model_observables(gens)              # from a list of finished runs
    chi2, terms = chi2_sample(row, obs)
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, os.pardir, "data", "reference_inventories.csv")

# measured columns that map onto a model quantity
LINKAGE_COLS = ["bo4", "b5", "bb", "b1", "c55", "c4o5", "dc"]
RATIO_KEYS = ["b5", "bb", "b1", "c55", "c4o5", "dc"]       # each divided by bo4

# per-method relative precision (docs/14 section 2); FLOOR is in the sample's own
# units (1 per 100 Ar or 1 % of linkages), the same floor calibrate_full.py used.
REL_SIGMA = {"13C-quant": 0.05, "HSQC-quant": 0.05, "HSQC-rel": 0.10,
             "degradative": 0.15, "legacy": 0.05}
FLOOR_SIGMA = 1.0
# The degree of condensation is a difference of large integrals and, in the
# words of the authors who report it (Balakshin 2016, Holzforschung 70:95), its
# quantification carries an error "in the range of 10-20 %". It also contains
# C6-condensed and S-condensed units the model cannot form, so the model is
# expected to sit below it (docs/14 section 3). It gets its own relative sigma.
REL_SIGMA_DC = 0.20


def _num(s):
    s = (s or "").strip()
    return float(s) if s else None


def load_reference_set(path=CSV_PATH, fit_only=True):
    """Rows of the reference CSV as dicts with numeric fields parsed.
    Blank measurements become None."""
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["sample_id"].startswith("#"):
                continue
            if fit_only and r.get("use_fit", "0").strip() != "1":
                continue
            for k in LINKAGE_COLS + ["sg_ratio", "h_frac"]:
                r[k] = _num(r.get(k))
            if r["bo4"] is None:
                continue                      # no ratio can be formed
            rows.append(r)
    return rows


def sigma_abs(row, key):
    """Absolute uncertainty of one measured quantity, in the row's own units."""
    rel = REL_SIGMA_DC if key == "dc" else REL_SIGMA.get(row["method"], 0.05)
    return max(rel * row[key], FLOOR_SIGMA)


def measured_ratios(row):
    """{key: (ratio, sigma_ratio)} for every quantity the sample reports."""
    x0, s0 = row["bo4"], sigma_abs(row, "bo4")
    out = {}
    for k in RATIO_KEYS:
        x = row[k]
        if x is None:
            continue
        s = sigma_abs(row, k)
        r = x / x0
        # first-order propagation, bo4 and x_k independent
        sr = r * ((s / x) ** 2 + (s0 / x0) ** 2) ** 0.5 if x > 0 else s / x0
        out[k] = (r, sr)
    return out


def model_observables(gens):
    """Model quantities per 100 units, pooled over a list of finished generator
    runs: the six bond counts and the mapped degree of condensation."""
    n_units = sum(g.summary()["n_units"] for g in gens)
    b = {k: sum(g.bonds[k] for g in gens) for k in ("bo4", "b5", "bb", "b1", "55", "4o5")}
    per100 = {"bo4": b["bo4"], "b5": b["b5"], "bb": b["bb"], "b1": b["b1"],
              "c55": b["55"], "c4o5": b["4o5"],
              "dc": 2 * b["55"] + b["4o5"] + b["b5"]}
    return {k: 100.0 * v / n_units for k, v in per100.items()}


def model_ratios(obs):
    return {k: obs[k] / obs["bo4"] for k in RATIO_KEYS} if obs["bo4"] > 0 else \
        {k: float("inf") for k in RATIO_KEYS}


def chi2_sample(row, obs):
    """chi-square of one sample against pooled model observables; returns
    (chi2, {key: z}) over the quantities the sample reports."""
    mr = model_ratios(obs)
    terms = {}
    for k, (r, sr) in measured_ratios(row).items():
        terms[k] = (mr[k] - r) / sr
    return sum(z * z for z in terms.values()), terms


def sim_inputs(row):
    """What the generator needs for this sample: S/G and the H fraction."""
    return {"sg": row["sg_ratio"], "pct_h": row["h_frac"] or 0.0}


if __name__ == "__main__":
    rows = load_reference_set(fit_only=False)
    print("%-28s %5s %5s %-11s %-12s  measured" % ("sample", "S/G", "fit", "method", "denominator"))
    for r in rows:
        meas = [k for k in LINKAGE_COLS if r[k] is not None]
        print("%-28s %5.2f %5s %-11s %-12s  %s" % (r["sample_id"], r["sg_ratio"], r["use_fit"],
                                                  r["method"], r["denominator"], " ".join(meas)))
    print(len(rows), "rows;", sum(1 for r in rows if r["use_fit"] == "1"), "in the fit")
