#!/usr/bin/env python
"""
mini_forward_generator.py

A fresh, minimal, fully-owned FORWARD lignin structure generator.

Motivation (see docs/06): Lignin-KMC demonstrated that a forward generator makes
the S/G->beta-O-4 relationship and chain linearity EMERGE, but its published Gani
barriers over-favour beta-5 for guaiacyl (softwood irreproducible) and pure
kinetic control cannot give high beta-O-4 + linearity + modest condensation at
once. This generator addresses both by owning its parameters:

  (P) coupling propensities are OUR tunable rate constants (recalibrate softwood);
  (M) supply is endwise (slow monomer feed) with an independent, low-rate
      oligomer-oligomer channel -> modest condensation without collapsing beta-O-4.

Chemistry (minimal but mechanistic):
  * Monomers G, S, H. Reactive positions per unit: O4 (phenolic O), C5, C1, and a
    one-shot beta (C8) radical donor. S is methoxylated at C5 -> NO free C5
    (the S-blocks-C5 rule), so S cannot form beta-5 / 5-5 / 4-O-5. H and G have C5.
  * Bonds form by radical cross-coupling between units on DIFFERENT fragments:
      beta-O-4 : donor.beta + acceptor.O4   (chain extension; donor.O4 = new end)
      beta-5   : donor.beta + acceptor.C5   (acceptor must have free C5)
      beta-1   : donor.beta + acceptor.C1
      beta-beta: donor.beta + donor2.beta   (monomer-monomer)
      5-5      : unitA.C5   + unitB.C5       (oligomer-oligomer branching)
      4-O-5    : unitA.O4   + unitB.C5       (oligomer-oligomer branching)
    Only free monomers carry a live beta donor, so beta-coupling is intrinsically
    endwise; 5-5 / 4-O-5 are the minority oligomer-oligomer channel.

This is a research prototype: propensities use aggregate counts and participants
are drawn at firing time with a different-fragment constraint (documented
simplification). Kinetic exactness is Lignin-KMC's job; this tool's job is a
transparent, controllable forward model.

Usage:
    python code/mini_forward_generator.py            # S/G sweep vs NMR
    python code/mini_forward_generator.py --calibrate  # (placeholder hook)
"""
import argparse
import os
import numpy as np

# ---- default OWNED parameters (tunable; this is the point) ------------------
DEFAULT_PARAMS = {
    # relative coupling rate constants (arbitrary units; only ratios matter)
    "k_bo4": 8.0,    # beta-O-4  (endwise chain extension)
    "k_b5":  2.0,    # beta-5    (needs acceptor free C5)  -> softwood knob (P)
    "k_b1":  0.4,    # beta-1
    "k_bb":  3.0,    # beta-beta (monomer-monomer)
    "k_55":  2.5,    # 5-5    (oligomer-oligomer, needs both C5) -> softwood condensation
    "k_4o5": 0.4,    # 4-O-5  (oligomer-oligomer, needs C5)
    "supply_rate": 0.4,   # (deprecated; supply is now deterministic)
    # DECOUPLED NUCLEATION (M): beta-beta is an EARLY, self-limiting event --
    # monomers dimerize at the start (nucleating chains), then growth is endwise.
    # n_nuclei sets how many dimers nucleate (=> beta-beta level + #chains);
    # pool_target stays ~1 for the endwise phase so chains grow LINEARLY and are
    # not forced to coalesce via condensed C-C bonds. This separates beta-beta
    # nucleation from ongoing condensation (which conflated under one big pool).
    "n_nuclei":  6.0,     # number of nucleation dimers formed up front
    "bb_s_boost": 4.0,    # syringyl radicals favour beta-beta (syringaresinol);
                          # boosts beta-beta weight when both nucleating units are S
                          # -> beta-beta rises with S/G, as observed in real lignin
    "pool_target": 1.0,   # standing free-monomer pool during endwise growth
    "oo_scale": 1.0,      # global multiplier on oligomer-oligomer channel (M)
    # ---- grass/palm (commelinid) extension; all default OFF (0) ----
    "pct_h":  0.0,        # fraction of p-hydroxyphenyl (H) monomers
    "frac_pca": 0.0,      # fraction of S units gamma-p-coumaroylated (BAHD PMT)
    "frac_tricin": 0.0,   # fraction of nucleation events initiated by tricin
}

MONO_C5_T = False         # tricin: flavone initiator, no free C5

MONO_C5 = {"G": True, "S": False, "H": True, "T": False}   # S/T have no free C5


class Unit:
    __slots__ = ("uid", "mtype", "o4", "c5", "c1", "beta")

    def __init__(self, uid, mtype):
        self.uid = uid
        self.mtype = mtype
        self.o4 = True                 # phenolic O available
        self.c5 = MONO_C5[mtype]       # C5 available (G/H only)
        self.c1 = True
        self.beta = True               # one-shot beta radical donor


class DSU:
    """union-find for fragment tracking"""
    def __init__(self):
        self.parent = {}

    def add(self, x):
        self.parent.setdefault(x, x)

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        self.parent[self.find(a)] = self.find(b)


class MiniForwardKMC:
    def __init__(self, sg_ratio, n_target=100, params=None, seed=0):
        self.sg = sg_ratio
        self.pct_s = sg_ratio / (1.0 + sg_ratio)
        self.n_target = n_target
        self.p = dict(DEFAULT_PARAMS)
        if params:
            self.p.update(params)
        self.rng = np.random.RandomState(seed)
        self.units = {}          # uid -> Unit
        self.dsu = DSU()
        self.n_introduced = 0
        self.bonds = {"bo4": 0, "b5": 0, "b1": 0, "bb": 0, "55": 0, "4o5": 0}
        self.branch_points = 0
        self.fired = {}       # debug: productive events by kind
        self.wasted = {}      # debug: selected-but-rejected events by kind
        # grass/palm extension bookkeeping (0 unless the extension params are set)
        self.comp = {"G": 0, "S": 0, "H": 0}
        self.n_pca = 0        # S units gamma-acylated by p-coumarate
        self.n_tricin = 0     # tricin chain initiators

    def _new_monomer(self):
        # composition: optional H fraction (pct_h); the rest split G/S by S/G
        r = self.rng.rand()
        pct_h = self.p.get("pct_h", 0.0)
        if r < pct_h:
            mtype = "H"
        else:
            mtype = "S" if self.rng.rand() < self.pct_s else "G"
        u = Unit(self.n_introduced, mtype)
        self.units[u.uid] = u
        self.dsu.add(u.uid)
        self.n_introduced += 1
        self.comp[mtype] += 1
        # grass/palm gamma-p-coumaroylation, predominantly on S units (BAHD PMT)
        if mtype == "S" and self.rng.rand() < self.p.get("frac_pca", 0.0):
            self.n_pca += 1
        return u

    # ---- site lists (recomputed each step; fine for minimal n) --------------
    def _sites(self):
        beta = [u for u in self.units.values() if u.beta]          # donors (monomers)
        o4 = [u for u in self.units.values() if u.o4]
        c5 = [u for u in self.units.values() if u.c5]
        c1 = [u for u in self.units.values() if u.c1]
        return beta, o4, c5, c1

    def _diff_frag_pairs(self, listA, listB):
        """count ordered/unordered pairs on different fragments (approx via product)."""
        # approximate: product minus same-fragment overlaps is expensive; use product
        # and rely on rejection sampling at firing. Return counts for propensity.
        return len(listA) * len(listB)

    def _count_free_monomers(self):
        """free monomer = beta-carrying unit that is still its own size-1 fragment."""
        size = {}
        for uid in self.units:
            r = self.dsu.find(uid); size[r] = size.get(r, 0) + 1
        return sum(1 for u in self.units.values()
                   if u.beta and size[self.dsu.find(u.uid)] == 1)

    def step(self):
        p = self.p
        # deterministic standing-pool refill: keep ~pool_target free monomers
        # available so early monomer-monomer beta-beta dimerization can occur.
        # (A slow KMC supply event cannot maintain a pool -- coupling is far
        # faster, so the pool would always collapse to ~0 and beta-beta never form.)
        if self.n_introduced < self.n_target:
            deficit = int(round(p["pool_target"])) - self._count_free_monomers()
            deficit = min(deficit, self.n_target - self.n_introduced)
            for _ in range(max(0, deficit)):
                self._new_monomer()

        beta, o4, c5, c1 = self._sites()
        n_beta = len(beta)

        # oligomer-oligomer (5-5, 4-O-5) require two units on DIFFERENT fragments;
        # count only inter-fragment pairs so single-fragment states give zero rate
        # (otherwise the KMC selects impossible events and spins on rejection).
        # fragment sizes: a beta radical is only a valid donor while the unit is a
        # FREE MONOMER (its own size-1 fragment). Once incorporated (even via a
        # C5/O4 acceptor coupling), its beta is no longer an available radical --
        # enforcing this is both correct endwise chemistry and prevents "stuck"
        # beta donors trapped inside the main fragment (which cannot satisfy the
        # different-fragment constraint and would spin forever).
        size = {}
        for uid in self.units:
            r = self.dsu.find(uid); size[r] = size.get(r, 0) + 1
        beta = [u for u in beta if size[self.dsu.find(u.uid)] == 1]
        n_beta = len(beta)

        o4_by, c5_by = {}, {}
        for u in o4:
            r = self.dsu.find(u.uid); o4_by[r] = o4_by.get(r, 0) + 1
        for u in c5:
            r = self.dsu.find(u.uid); c5_by[r] = c5_by.get(r, 0) + 1
        t_o4, t_c5 = len(o4), len(c5)
        inter_c5_pairs = (t_c5 * t_c5 - sum(v * v for v in c5_by.values())) / 2.0
        inter_o4_c5 = t_o4 * t_c5 - sum(o4_by.get(r, 0) * c5_by.get(r, 0)
                                        for r in set(o4_by) | set(c5_by))

        # Endwise chemistry: an incoming monomer beta-radical couples to a
        # phenolic END of the polymer, offering that end unit's O4 (->beta-O-4),
        # C5 (->beta-5, G/H only) or C1 (->beta-1). So beta-5 / beta-1 acceptors
        # are restricted to END units (open O4), NOT arbitrary internal C5/C1 --
        # otherwise the huge internal-C5 pool inflates beta-5 regardless of rate
        # constants. Internal C5 still participates in oligomer-oligomer 5-5/4-O-5.
        ends = o4                              # units with an open phenolic O4
        b5_acc = [u for u in ends if u.c5]
        b1_acc = [u for u in ends if u.c1]

        # propensities (aggregate). beta-couplings: donor(free monomer) x end.
        # (Supply is handled deterministically above, not as a KMC event.)
        props = {}
        props["bo4"] = p["k_bo4"] * n_beta * len(ends)
        props["b5"] = p["k_b5"] * n_beta * len(b5_acc)
        props["b1"] = p["k_b1"] * n_beta * len(b1_acc)
        props["bb"] = p["k_bb"] * n_beta * max(n_beta - 1, 0) / 2.0
        props["55"] = p["k_55"] * p["oo_scale"] * inter_c5_pairs
        props["4o5"] = p["k_4o5"] * p["oo_scale"] * max(inter_o4_c5, 0.0)

        total = sum(props.values())
        if total <= 0:
            return False
        r = self.rng.rand() * total
        acc = 0.0
        choice = None
        for k, v in props.items():
            acc += v
            if r <= acc:
                choice = k
                break

        before = sum(self.bonds.values())
        res = self._fire_coupling(choice, beta, o4, c5, c1, b5_acc, b1_acc)
        if sum(self.bonds.values()) > before:
            self.fired[choice] = self.fired.get(choice, 0) + 1
        else:
            self.wasted[choice] = self.wasted.get(choice, 0) + 1
        return res

    def _pick_diff(self, listA, listB, same_ok=False, tries=12):
        """pick a on A, b on B, on different fragments (and different units)."""
        for _ in range(tries):
            a = listA[self.rng.randint(len(listA))]
            b = listB[self.rng.randint(len(listB))]
            if a.uid == b.uid:
                continue
            if self.dsu.find(a.uid) == self.dsu.find(b.uid):
                continue
            return a, b
        return None, None

    def _fire_coupling(self, kind, beta, o4, c5, c1, b5_acc, b1_acc):
        if kind == "bo4":
            d, a = self._pick_diff(beta, o4)
            if d is None:
                return True
            d.beta = False
            a.o4 = False              # ether formed on acceptor's O4
            self.dsu.union(d.uid, a.uid)
            self.bonds["bo4"] += 1
        elif kind == "b5":
            d, a = self._pick_diff(beta, b5_acc)
            if d is None:
                return True
            d.beta = False
            a.c5 = False
            self.dsu.union(d.uid, a.uid)
            self.bonds["b5"] += 1
        elif kind == "b1":
            d, a = self._pick_diff(beta, b1_acc)
            if d is None:
                return True
            d.beta = False
            a.c1 = False
            self.dsu.union(d.uid, a.uid)
            self.bonds["b1"] += 1
        elif kind == "bb":
            d, d2 = self._pick_diff(beta, beta)
            if d is None:
                return True
            d.beta = False
            d2.beta = False
            self.dsu.union(d.uid, d2.uid)
            self.bonds["bb"] += 1
        elif kind == "55":
            a, b = self._pick_diff(c5, c5)
            if a is None:
                return True
            a.c5 = False
            b.c5 = False
            self.dsu.union(a.uid, b.uid)
            self.bonds["55"] += 1
            self.branch_points += 1
        elif kind == "4o5":
            a, b = self._pick_diff(o4, c5)
            if a is None:
                return True
            a.o4 = False
            b.c5 = False
            self.dsu.union(a.uid, b.uid)
            self.bonds["4o5"] += 1
            self.branch_points += 1
        return True

    def _nucleate_one(self):
        """Form one nucleation dimer from two fresh monomers (donor beta + acceptor
        position), with the coupling type drawn from the monomer-monomer rate
        constants. This is where beta-beta predominantly comes from."""
        if self.n_introduced + 2 > self.n_target:
            return False
        p = self.p
        # tricin chain initiation (grass/palm): a monolignol beta couples to
        # tricin's 4'-O (beta-O-4-like); tricin caps the chain start (no further
        # reactivity). A grass-specific nucleation mechanism.
        if self.rng.rand() < p.get("frac_tricin", 0.0):
            d = self._new_monomer()
            t = Unit(self.n_introduced, "T")
            t.c1 = False; t.beta = False          # tricin: only its 4'-O couples
            self.units[t.uid] = t
            self.dsu.add(t.uid)
            self.n_introduced += 1
            self.n_tricin += 1
            d.beta = False
            t.o4 = False
            self.dsu.union(d.uid, t.uid)
            self.bonds["bo4"] += 1
            return True
        d = self._new_monomer()
        a = self._new_monomer()
        bb_w = p["k_bb"]
        if d.mtype == "S" and a.mtype == "S":
            bb_w *= p.get("bb_s_boost", 1.0)      # syringaresinol favoured
        opts = [("bb", bb_w), ("bo4", p["k_bo4"]), ("b1", p["k_b1"])]
        if a.c5:
            opts.append(("b5", p["k_b5"]))
        tot = sum(w for _, w in opts)
        r = self.rng.rand() * tot
        acc, kind = 0.0, opts[-1][0]
        for k, w in opts:
            acc += w
            if r <= acc:
                kind = k
                break
        d.beta = False
        if kind == "bb":
            a.beta = False
        elif kind == "bo4":
            a.o4 = False
        elif kind == "b5":
            a.c5 = False
        elif kind == "b1":
            a.c1 = False
        self.dsu.union(d.uid, a.uid)
        self.bonds[kind] += 1
        return True

    def run(self, max_steps=200000):
        # nucleation phase: form n_nuclei dimers up front (chain starts + beta-beta)
        n_nuc = max(1, int(round(self.p.get("n_nuclei", 1))))
        for _ in range(n_nuc):
            if not self._nucleate_one():
                break
        if self.n_introduced == 0:            # safety: ensure something exists
            self._new_monomer()
        # endwise growth phase: pool_target ~1 keeps chains linear
        steps = 0
        while steps < max_steps:
            if not self.step():
                break
            steps += 1
        return self.summary()

    def summary(self):
        total_bonds = sum(self.bonds.values())
        frac = {k: (v / total_bonds if total_bonds else 0.0)
                for k, v in self.bonds.items()}
        # fragment sizes
        sizes = {}
        for uid in self.units:
            root = self.dsu.find(uid)
            sizes[root] = sizes.get(root, 0) + 1
        n_units = len(self.units)
        branch_coeff = self.branch_points / n_units if n_units else 0.0
        comp = {k: (v / n_units if n_units else 0.0) for k, v in self.comp.items()}
        return {"frac": frac, "n_bonds": total_bonds, "n_units": n_units,
                "n_frags": len(sizes), "max_frag": max(sizes.values()) if sizes else 0,
                "branch_coeff": branch_coeff, "comp": comp,
                "pca_frac": self.n_pca / n_units if n_units else 0.0,
                "n_tricin": self.n_tricin}


# ---- driver: S/G sweep, compared to NMR anchors -----------------------------
SG_RATIOS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
NMR = {  # rough whole-cell-wall anchors for context (beta-O-4 %)
    0.0: 46, 1.0: 60, 2.0: 72, 3.0: 79,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=120, help="monomers per polymer")
    ap.add_argument("--reps", type=int, default=25)
    ap.add_argument("--seed", type=int, default=20260714)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, os.pardir, "results")
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "mini_forward_sweep.csv")

    print("Mini forward generator — S/G sweep (endwise + minority oligomer-oligomer)")
    print("params:", DEFAULT_PARAMS)
    print("n=%d monomers, reps=%d\n" % (args.n, args.reps))
    print("%5s %8s %7s %6s %6s %6s %6s %7s %6s   %s" %
          ("S/G", "bO4", "b5", "bb", "55", "4o5", "b1", "branch", "NMR", "note"))

    rows = []
    for sg in SG_RATIOS:
        acc = {k: [] for k in ["bo4", "b5", "b1", "bb", "55", "4o5"]}
        branches = []
        for rep in range(args.reps):
            gen = MiniForwardKMC(sg, n_target=args.n, seed=args.seed + rep + int(sg * 100))
            s = gen.run()
            for k in acc:
                acc[k].append(s["frac"][k])
            branches.append(s["branch_coeff"])
        m = {k: 100 * np.mean(acc[k]) for k in acc}
        br = np.mean(branches)
        nmr = NMR.get(sg, None)
        note = "SW anchor" if sg == 0 else ("" if nmr is None else "HW anchor")
        rows.append((sg, m, br, nmr))
        print("%5.2f %7.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %7.3f %5s   %s" %
              (sg, m["bo4"], m["b5"], m["bb"], m["55"], m["4o5"], m["b1"], br,
               ("%d" % nmr) if nmr else "-", note))

    with open(out_csv, "w") as fh:
        fh.write("sg_ratio,bo4,b5,bb,55,4o5,b1,branch_coeff,nmr_bo4\n")
        for sg, m, br, nmr in rows:
            fh.write("%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s\n" %
                     (sg, m["bo4"], m["b5"], m["bb"], m["55"], m["4o5"], m["b1"], br,
                      nmr if nmr else ""))
    print("\nWrote %s" % out_csv)


if __name__ == "__main__":
    main()
