#!/usr/bin/env python
"""
mini_forward_generator.py

A fresh, minimal, fully-owned FORWARD lignin structure generator.

Lignin, the final frontier. These are the runs of the forward generator. Its
continuing mission: to explore strange new feedstocks, to seek out new linkage
distributions and new condensations, to boldly polymerize where no model has
polymerized before.

Motivation (see docs/06): Lignin-KMC demonstrated that a forward generator makes
the S/G->beta-O-4 relationship and chain linearity EMERGE, but its published Gani
barriers over-favour beta-5 for guaiacyl (softwood irreproducible) and pure
kinetic control cannot give high beta-O-4 + linearity + modest condensation at
once. This generator addresses both by owning its parameters:

  (P) coupling propensities are OUR tunable rate constants (recalibrate softwood);
  (M) supply is endwise (slow monomer feed) with an independent, low-rate
      INTER-FRAGMENT channel -> modest condensation without collapsing beta-O-4.

NOTE on (M)'s name: this channel requires only that its two partners lie on
DIFFERENT fragments. A free monomer is itself a fragment, so ~66% of 5-5 / 4-O-5
events actually couple a free monomer through C5 rather than two oligomers (100%
when n_nuclei = 1, where no second oligomer exists). In standard nomenclature that
is monomer-monomer (MM) and monomer-polymer (MP) coupling, not polymer-polymer
(PP), so the channel is called INTER-FRAGMENT. The name was wrong; the chemistry is
not. Thioacidolysis of pine and spruce reports 5-5 as formable by MM, MP and PP
coupling alike -- unlike phenolic beta-5 and beta-1, which are MP-only (Saito et al.
2021, ChemSusChem 14:2554) -- and the 5-5 dimer of two coniferyl alcohols is a
competent precursor that polymerizes from both side chains (Shigetomi et al. 2025).

oo_true = 1 restricts the channel to PP only. It is a DIAGNOSTIC, NOT A FIX: it
improves poplar and high-S and BREAKS softwood, because PP bonds each merge two
fragments and so cap at n_nuclei - 1. That failure is the point -- it is the
biochemistry companion's central result (softwood cannot be condensing via PP,
because every chain a PP bond consumes was paid for in a resinol softwood does not
have). Do not enable it to make the model "more correct"; it makes it less so. See
../../Lignin_BioChem/docs/01_findings.md and
../docs/07_correction_condensed_channel.md.

Chemistry (minimal but mechanistic):
  * Monomers G, S, H. Reactive positions per unit: O4 (phenolic O), C5, C1, and a
    one-shot beta (C8) radical donor. S is methoxylated at C5 -> NO free C5
    (the S-blocks-C5 rule), so S cannot form beta-5 / 5-5 / 4-O-5. H and G have C5.
  * Bonds form by radical cross-coupling between units on DIFFERENT fragments:
      beta-O-4 : donor.beta + acceptor.O4   (chain extension; donor.O4 = new end)
      beta-5   : donor.beta + acceptor.C5   (acceptor must have free C5)
      beta-1   : donor.beta + acceptor.C1
      beta-beta: donor.beta + donor2.beta   (monomer-monomer)
      5-5      : unitA.C5   + unitB.C5       (inter-fragment branching)
      4-O-5    : unitA.O4   + unitB.C5       (inter-fragment branching)
    Only free monomers carry a live beta donor, so beta-coupling is intrinsically
    endwise; 5-5 / 4-O-5 are the minority inter-fragment channel (see NOTE above).

This is a research prototype. There is no time variable: at each step one
reaction is drawn with probability proportional to k times the number of pairs
that could form it right now, and the pair is then drawn uniformly from that
same list. Only ratios of the constants matter. Kinetic exactness is Lignin-KMC's
job; this tool's job is a transparent, controllable forward model.

v1.1 (2026-09, docs/12_code_audit.md) corrected three things the v1.0.1 code did
that its own chemistry forbids:
  A. beta-5 and beta-1 now consume the acceptor's phenolic 4-OH (phenylcoumaran
     ring closure / spirodienone), so a beta-5 acceptor cannot later accept a
     beta-O-4 on the same oxygen. v1.0.1 left it free, doubling the phenolic-end
     count.
  B. 5-5 and 4-O-5 now require a PHENOLIC C5 partner (a phenoxy radical needs a
     free 4-OH). v1.0.1 let etherified internal units condense.
  C. partners are drawn exactly from the valid inter-fragment pairs. v1.0.1 drew
     from the full site lists and discarded the event after 12 same-fragment
     draws, which scaled the condensed constants by ~0.5 at softwood and ~0.95
     at high-S.
At the v1.0.1 constants the mean error is unchanged (1.62 -> 1.60 points) but the
mechanism is not: phenolic ends drop from ~35 to ~12 per 120 units at softwood,
and at high-S the corrected chemistry no longer always welds the chains into one
molecule (S ends have no C5). See docs/12_code_audit.md section 3.

Usage:
    python code/mini_forward_generator.py            # S/G sweep vs NMR
    python code/mini_forward_generator.py --calibrate  # (placeholder hook)
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import rel

# ---- default OWNED parameters (tunable; this is the point) ------------------
DEFAULT_PARAMS = {
    # relative coupling rate constants (arbitrary units; only ratios matter)
    "k_bo4": 8.0,    # beta-O-4  (endwise chain extension)
    "k_b5":  2.0,    # beta-5    (needs acceptor free C5)  -> softwood knob (P)
    "k_b1":  0.4,    # beta-1
    "k_bb":  3.0,    # beta-beta (monomer-monomer)
    "k_55":  2.5,    # 5-5    (inter-fragment, needs both C5) -> softwood condensation
    "k_4o5": 0.4,    # 4-O-5  (inter-fragment, needs C5)
    "supply_rate": 0.4,   # (deprecated; supply is now deterministic)
    # DECOUPLED NUCLEATION (M): beta-beta is an EARLY, self-limiting event --
    # monomers dimerize at the start (nucleating chains), then growth is endwise.
    # n_nuclei sets how many dimers nucleate, which fixes BOTH the beta-beta level
    # and the number of fragments available to condense -- the two are coupled, and
    # sweeping n_nuclei moves them together (r = +0.96..+0.99), not in opposition.
    # pool_target stays ~1 for the endwise phase so chains grow LINEARLY during
    # growth. They are NOT left separate, though: beta-coupling can never join two
    # chains (only free monomers carry a live beta donor), so the inter-fragment
    # channel welds every nucleus into ONE macromolecule -- n_frags = 1 at any
    # n_nuclei, and even at oo_scale = 0.01. At oo_scale = 0, n_frags = n_nuclei
    # exactly. The model's DP is therefore n_target, a simulation size: make no
    # chain-count or molar-mass claim from it.
    # This still separates beta-beta nucleation from ongoing condensation (which
    # conflated under one big pool), which is what the phase separation buys.
    "n_nuclei":  6.0,     # number of nucleation dimers formed up front
    "k_bo4_nuc": 1.0,     # beta-O-4 at the MONOMER-MONOMER (nucleation) step.
                          # Separate from k_bo4 because nucleation is a different
                          # regime: k_bo4 describes a monomer adding to a chain
                          # END, which is not the same reaction as two free
                          # monomers coupling. Gani (docs/04) puts mon-mon
                          # beta-O-4 at 6.3 kcal/mol against beta-beta's 5.2, i.e.
                          # beta-O-4 does NOT dominate monomer-monomer the way it
                          # dominates chain extension -- it is ~6x SLOWER than
                          # beta-beta there. Default ~1 reflects that; the fitted
                          # value is 1.09 (results/fitted_params.json). It was
                          # inheriting k_bo4 = 8 until 2026-07-15, which is what
                          # suppressed beta-beta. See _pair_mod and docs/06 Exp 6.
    "pool_target": 1.0,   # standing free-monomer pool during endwise growth
    "oo_scale": 1.0,      # global multiplier on the inter-fragment channel (M)
    # ---- mechanism switches (>=0.5 is ON). Both default ON = the calibrated
    # model; turning them OFF is how the two load-bearing mechanistic choices are
    # ablated rather than merely asserted. See Lignin_BioChem/code/ablate_mechanisms.py
    "end_restricted": 1.0,   # beta-5/beta-1 acceptors limited to phenolic END units
    "nucleation_on": 1.0,    # run the up-front monomer-monomer nucleation phase
    # oo_true = 1 restricts 5-5 / 4-O-5 to PP (both partners in fragments of size
    # > 1). Default 0 = the calibrated inter-fragment channel, in which ~66% of
    # condensed events couple a free monomer through C5 (see the module NOTE).
    #
    # This is a DIAGNOSTIC, not a correction. PP-only is not the truer chemistry:
    # Saito 2021 finds no 5-5 in pine or spruce carrying a free phenolic OH, so both
    # sides are already etherified -- the biphenyl formed early and the polymer grew
    # outward from it, which is a monomer-level origin, not a late join between two
    # long chains. Enabling this breaks softwood (RMSE 2.25 -> 7.39) precisely
    # because PP bonds cap at n_nuclei - 1; that failure is a result, not a bug.
    "oo_true": 0.0,
    # ---- route resolution for 5-5 (see docs/07 option (d)) --------------------
    # The inter-fragment channel lumps three chemically distinct reactions that
    # differ in whether a condensed bond CONSUMES a chain:
    #   MM (both partners free monomers)  -> CREATES a fragment (a bisCA dimer)
    #   MP (one monomer, one oligomer)    -> fragment count UNCHANGED
    #   PP (both oligomers)               -> CONSUMES a fragment; capped at N-1
    # Only PP is bounded, which is why oo_true (= PP-only) breaks softwood. Setting
    # oo_route resolves 5-5 into one route so the three can be compared as matched
    # hypotheses rather than lumped:
    #   0 = any inter-fragment (default; the calibrated channel, an MM/MP/PP mix)
    #   1 = PP only (identical to oo_true)
    #   2 = MP only
    # NOTE: with pool_target = 1 there is at most ONE free monomer during growth, so
    # MM 5-5 cannot occur in the growth phase at all -- a bisCA dimer can only form
    # at nucleation. That is what k_55_nuc is for.
    # 4-O-5 is deliberately NOT split: Saito et al. 2021 put it at <=1% of lignin
    # units and find no evidence it is further etherified, i.e. it does not form the
    # branching structure usually envisioned, so the MM/MP/PP question is about 5-5.
    "oo_route": 0.0,
    # 5-5 formed at NUCLEATION from two free monomers -- the bisCA route. Creates a
    # fragment carrying a condensed bond, so unlike PP it costs no resinol. Default
    # 0 = off (the published model cannot nucleate a biphenyl at all).
    # Shigetomi et al. 2025 synthesized this dimer and showed it polymerizes with
    # chain growth from BOTH side chains.
    "k_55_nuc": 0.0,
    # ---- CONTINUOUS nucleation -----------------------------------------------
    # n_nuclei front-loads every nucleation event into a burst before growth
    # starts. That was a modelling convenience, and it has a consequence: the
    # number of chains -- and hence every quantity that scales with it -- is capped
    # at n_nuclei, decided before the first monomer is added.
    #
    # Real lignification has no burst. Monolignols arrive at the wall continuously,
    # and occasionally two meet and dimerize (starting a NEW chain) instead of
    # adding to an existing one. Chang et al. 2020 infer exactly this for biphenyls,
    # which are "formed while monolignols are being transported to the lignifying
    # zone" -- i.e. throughout the supply stream, not in one burst.
    #
    # k_nuc is the propensity of a nucleation event DURING growth: two fresh
    # monomers are introduced and coupled by the monomer-monomer constants (so it
    # can produce beta-beta, beta-5, beta-O-4, beta-1, or -- with k_55_nuc -- a
    # bisCA biphenyl). Rate is taken constant because the standing monomer pool is
    # held constant, so [M]^2 is too. Default 0 = the published front-loaded model.
    "k_nuc": 0.0,
    # ---- oligomer beta-beta: re-radicalised chain-end resinols ----------------
    # The beta-beta ceiling may be self-inflicted rather than chemical: the model
    # lets only FREE MONOMERS carry a live beta radical (so beta-beta is monomer-
    # monomer at nucleation only -> beta-beta% = p_bb * n_nuclei / (N-1), a hard cap),
    # whereas peroxidase/laccase continuously RE-RADICALISE oligomer units, enabling
    # monomer-oligomer and oligomer-oligomer beta-beta that would lift the cap.
    # oligo_bb (>0) opens exactly that channel: phenolic END units of oligomers
    # (accessible, unlike buried internal units whose beta the growth mechanism has
    # already spent) are re-radicalised at beta and may couple beta-beta, JOINING two
    # fragments. Rate = k_bb * oligo_bb * (n_free_monomer x n_ends + inter-fragment
    # end pairs). Default 0 = the published monomer-only model. Judged like k_ao4 /
    # k_nuc: whether a refit accepts it and lowers held-out error, or collapses it to
    # the floor (which would show the ceiling is chemistry, not an artifact).
    "oligo_bb": 0.0,
    # ---- Poisson nucleus count ------------------------------------------------
    # n_nuclei is normally rounded to a fixed integer, identical in every molecule.
    # That makes the beta-beta ceiling a hard line: beta-beta% = p_bb * N / (N_tot-1)
    # with N fixed, so poplar (5.13%) and high-S (4.57%) are unreachable at N = 5.
    # With this switch each MOLECULE draws its own N ~ Poisson(n_nuclei), so the
    # nucleus count -- and hence the ceiling -- varies rep to rep and the line
    # becomes a BAND (also what SI Fig. 6 wants). This is a MEAN-preserving change:
    # E[N] = n_nuclei exactly, so it tests whether integer rounding, rather than the
    # mean nucleus count, was what pinned beta-beta below the two hardwood anchors.
    # Default 0 = the published fixed-integer model.
    "nuclei_poisson": 0.0,
    # ---- alpha-O-4: the ETHER route to joining two chains ---------------------
    # Every beta-O-4 coupling leaves a quinone methide at the donor's alpha
    # carbon. Its fate is a nucleophilic addition, NOT a radical coupling -- and at
    # physiological pH the preferred nucleophile is a PHENOL, not water:
    #   "under neutral conditions (pH 6.5) the predominant reaction is addition of
    #    phenolic hydroxyl to form benzyl non-cyclic aryl ethers"  (Brunow 1989)
    #   "...the quinone methide ... rapidly led to the formation of significant
    #    amounts of the non-cyclic alpha-O-4 benzyl aryl ethers ... via phenol
    #    addition. The addition of water became predominant only under acid
    #    catalysis."  (Quideau & Ralph 1994, whose adducts include alpha-pinoresinol
    #    and alpha-dehydrodiconiferyl ethers -- i.e. the phenol donor is an OLIGOMER)
    # Chang 2020 counts noncyclic alpha-O-4' at ~4/100 C9 in softwood, a major
    # branch point; Balakshin 2020 finds a large share of spruce MWL's crosslinks
    # sit in the SIDE CHAINS rather than on aromatic rings.
    #
    # Why it matters here: alpha-O-4 JOINS TWO FRAGMENTS WITHOUT SPENDING A
    # CONDENSED C-C BOND. Every argument in this project that treats 5-5 / 4-O-5 as
    # the currency of chain-joining assumes this channel does not exist. It does.
    # Because it is not a radical coupling, the combinatorial account -- and
    # Lignin-KMC, and this generator until now -- omits it by construction.
    # Default 0 = off, preserving the published calibration.
    "k_ao4": 0.0,
    # ---- grass/palm (commelinid) extension; all default OFF (0) ----
    "pct_h":  0.0,        # fraction of p-hydroxyphenyl (H) monomers
    "frac_pca": 0.0,      # fraction of S units gamma-p-coumaroylated (BAHD PMT)
    "frac_tricin": 0.0,   # fraction of nucleation events initiated by tricin
    # ---- beta-5 donor specificity (tested 2026-09-15 as an eighth parameter) -------
    # The corrected-set fit over-predicts beta-5 for high-S woods and under-predicts it
    # for low-S woods (docs/15 section 4.6). The acceptor side is already in the model
    # (an S end has no C5); what is missing is the DONOR side: a sinapyl beta radical
    # may add to a C5 less readily than a coniferyl one. b5_s_donor multiplies the
    # chain-extension beta-5 propensity when the free monomer is S. 1.0 = no effect
    # (the published model). Nucleation is untouched: Gani's monomer-monomer beta-5
    # barriers are equal for (G,G) and (S,G) donors, so the DFT gives no donor effect
    # there.
    "b5_s_donor": 1.0,
}

MONO_C5_T = False         # tricin: flavone initiator, no free C5

MONO_C5 = {"G": True, "S": False, "H": True, "T": False}   # S/T have no free C5

# ---- Gani pair-specificity at the MONOMER-MONOMER (nucleation) step ----------
# Nucleation couples two FREE MONOMERS. That is a different reaction from a
# monomer adding to a chain end, and Gani's DFT barriers (docs/04, the mon,mon
# column) say so plainly: beta-beta is the FASTEST monomer-monomer channel
# (5.2 kcal/mol, G-G and S-S alike), while beta-O-4 is middling for G-G (6.3) and
# COLLAPSES for syringyl (9.8). Reusing the chain-extension constants here, as
# this generator used to, silently asserts that beta-O-4 dominates nucleation the
# way it dominates chain growth. It does not, and that was the beta-beta deficit.
#
# We take ONLY THE PAIR-SPECIFICITY from Gani -- each channel's barrier relative
# to its own (G,G) value -- and leave every absolute scale to our own fitted
# constants. This is deliberate. docs/06 Exp 1-2 established that Gani's ABSOLUTE
# monomer-oligomer beta-5/beta-O-4 balance is wrong (gap P), so we do not import
# it; normalising to (G,G) means g == 1 there and the fitted constants keep their
# meaning. What we do import is the pattern -- syringyl kills beta-O-4 and (having
# no free C5) beta-5, leaving beta-beta to win by default -- which is the same
# S-blocks-C5 chemistry the rest of the model rests on, and is why syringaresinol
# is abundant in high-S lignin while pinoresinol is minor in softwood.
#
# This replaces `bb_s_boost`, a flat x4 fudge on S-S beta-beta: the S/G trend in
# beta-beta is now a consequence of the competing channels failing, not a knob.
RT_298 = 0.001987 * 298.15        # kcal/mol at 298.15 K (the barriers' reference)

_GANI_MM = {   # kcal/mol, monomer-monomer; docs/04 DEF_E_BARRIER_KCAL_MOL
    "bb":  {("G", "G"): 5.2, ("S", "S"): 5.2, ("S", "G"): 6.5, ("G", "S"): 6.5},
    "bo4": {("G", "G"): 6.3, ("S", "S"): 9.8, ("S", "G"): 9.1, ("G", "S"): 8.9},
    "b5":  {("G", "G"): 5.5, ("S", "G"): 5.5},   # acceptor needs a free C5
    "b1":  {("G", "G"): 9.6, ("S", "S"): 11.9, ("S", "G"): 11.7, ("G", "S"): 10.7},
}


def _pair_mod(kind, donor_type, acceptor_type):
    """Gani monomer-monomer pair modifier for `kind`, normalised to the (G,G)
    pair so that g(G,G) == 1 and the fitted constants are unchanged there.

    Returns 1.0 for pairs Gani does not tabulate (H, tricin) -- those are the
    commelinid extension's business and have no DFT entry to borrow.
    """
    tbl = _GANI_MM[kind]
    barrier = tbl.get((donor_type, acceptor_type))
    if barrier is None:
        return 1.0
    return float(np.exp(-(barrier - tbl[("G", "G")]) / RT_298))


class Unit:
    __slots__ = ("uid", "mtype", "o4", "c5", "c1", "beta", "alpha", "deg")

    def __init__(self, uid, mtype):
        self.uid = uid
        self.mtype = mtype
        self.deg = 0                   # bonds this unit carries (>=3 = branch point)
        self.o4 = True                 # phenolic O available
        self.c5 = MONO_C5[mtype]       # C5 available (G/H only)
        # beta-1 acceptor: a phenolic unit whose own side chain is already in a
        # beta-O-4 linkage -- its alpha-OH traps the quinone methide and closes the
        # spirodienone (Zhang & Gellerstedt 2001; Zhang, Gellerstedt, Ralph & Lu 2006). A fresh monomer is not one;
        # a unit becomes one when it donates a beta-O-4 (fix D, 2026-09-16,
        # Lignin_Structure_Toolkit/docs/notes_for_paper1_beta1.md).
        self.c1 = False
        self.beta = True               # one-shot beta radical donor
        # alpha: a quinone methide, created when THIS unit donates a beta-O-4.
        # Capturable by a phenol on another fragment -> alpha-O-4 (see k_ao4).
        # Modelled as persistent; real QMs are transient (t_1/2 ~ 13-15 min for
        # G-type at pH 7, Zhu 2019), so this overstates the window. A documented
        # simplification: the rate constant absorbs the capture probability.
        self.alpha = False


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
        self.bonds = {"bo4": 0, "b5": 0, "b1": 0, "bb": 0, "55": 0, "4o5": 0,
                      "ao4": 0}
        self.branch_points = 0
        self.n_nuclei_formed = 0   # nucleation dimers made in the up-front burst
        self.n_nucleation_events = 0  # ALL nucleation events (burst + continuous k_nuc)
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
        alpha = [u for u in self.units.values() if u.alpha]        # quinone methides
        c5 = [u for u in self.units.values() if u.c5]
        c1 = [u for u in self.units.values() if u.c1]
        return beta, o4, c5, c1, alpha

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
            target = int(round(p["pool_target"]))
            # A lone free monomer with no oligomer to add onto cannot react -- its
            # beta radical needs a DIFFERENT fragment, and there is none -- so the
            # pool would spin on rejection. This only arises in degenerate starts
            # where no dimer seeded growth (Poisson nucleation drawing 0 nuclei, or
            # nucleation_on = 0). Guarantee >= 2 free monomers then, so a chain can
            # nucleate from the pool. With the default burst (>= 1 dimer) an oligomer
            # always exists here, so target stays pool_target and the calibrated path
            # is untouched (verified bit-identical).
            sizes = {}
            for uid in self.units:
                r = self.dsu.find(uid); sizes[r] = sizes.get(r, 0) + 1
            if not any(v > 1 for v in sizes.values()):
                target = max(target, 2)
            deficit = target - self._count_free_monomers()
            deficit = min(deficit, self.n_target - self.n_introduced)
            for _ in range(max(0, deficit)):
                self._new_monomer()

        beta, o4, c5, c1, alpha = self._sites()

        # A beta radical is a valid donor only while its unit is a FREE MONOMER
        # (its own size-1 fragment). Once incorporated -- even through a C5/O4
        # coupling -- the side chain is spent. That is what makes beta-coupling
        # endwise, and it is why a monomer that attached through C5 keeps an
        # unreacted cinnamyl-alcohol side chain (reported by summary()).
        size = {}
        for uid in self.units:
            r = self.dsu.find(uid); size[r] = size.get(r, 0) + 1
        self._size = size
        beta = [u for u in beta if size[self.dsu.find(u.uid)] == 1]

        # Coupling at C5 needs a phenoxy radical, i.e. a unit whose 4-OH is still
        # free. Internal etherified units keep a chemically free C5 but cannot be
        # oxidised, so they do not condense (fix B; the same rule that restricts
        # beta-5 to phenolic ends below).
        c5_phen = [u for u in c5 if u.o4]

        # Candidate partners for the condensed channels. Under oo_true / oo_route=1
        # these exclude free monomers, so 5-5 / 4-O-5 couple two genuine oligomers.
        # Otherwise (default) any phenolic unit qualifies and the standing free
        # monomer participates. See DEFAULT_PARAMS for why the routes differ.
        route = int(round(p.get("oo_route", 0.0)))
        if p.get("oo_true", 0.0) >= 0.5:
            route = 1
        if route == 1:
            o4_oo = [u for u in o4 if size[self.dsu.find(u.uid)] > 1]
            c5_oo = [u for u in c5_phen if size[self.dsu.find(u.uid)] > 1]
        else:
            o4_oo, c5_oo = o4, c5_phen
        c5_mono = [u for u in c5_phen if size[self.dsu.find(u.uid)] == 1]
        c5_olig = [u for u in c5_phen if size[self.dsu.find(u.uid)] > 1]

        # Endwise chemistry: an incoming monomer beta-radical couples to a
        # phenolic END of the polymer, offering that end unit's O4 (->beta-O-4),
        # C5 (->beta-5, G/H only) or C1 (->beta-1). So beta-5 / beta-1 acceptors
        # are restricted to END units (open O4), NOT arbitrary internal C5/C1 --
        # otherwise the huge internal-C5 pool inflates beta-5 regardless of rate
        # constants.
        ends = o4                              # units with an open phenolic O4
        if p.get("end_restricted", 1.0) >= 0.5:
            b5_acc = [u for u in ends if u.c5]
            b1_acc = [u for u in ends if u.c1]
        else:
            # ABLATION: every free C5/C1 in the polymer accepts, whether or not it
            # sits at a phenolic end. This is the "all sites open" chemistry the
            # end-restriction rule denies; the internal-C5 pool then dwarfs the end
            # pool and beta-5 inflates regardless of k_b5.
            b5_acc = c5
            b1_acc = c1

        # Each reaction is the list of pairs that could form it right now: distinct
        # units on different fragments. The propensity is k times the length of
        # that list, and the event, once chosen, is drawn uniformly from the same
        # list, so the count and the draw agree exactly (fix C).
        pairs = {}
        pairs["bo4"] = self._pairs(beta, ends)
        pairs["b5"] = self._pairs(beta, b5_acc)
        pairs["b1"] = self._pairs(beta, b1_acc)
        pairs["bb"] = self._pairs(beta, beta, unordered=True)
        # oligomer beta-beta (re-radicalised chain ends). Off by default; the test
        # for whether the beta-beta ceiling is an artifact. Phenolic END units of
        # oligomers re-radicalise at beta and couple beta-beta with a free monomer
        # or with another oligomer end, joining two fragments.
        if p.get("oligo_bb", 0.0) > 0.0:
            pool = beta + [u for u in o4 if size[self.dsu.find(u.uid)] > 1]
            pairs["bb_olig"] = self._pairs(pool, pool, unordered=True)
        if route == 2:      # MP only: one free monomer x one oligomer
            pairs["55"] = self._pairs(c5_mono, c5_olig)
        else:
            pairs["55"] = self._pairs(c5_oo, c5_oo, unordered=True)
        pairs["4o5"] = self._pairs(o4_oo, c5_oo)
        # alpha-O-4: a quinone methide captures a phenol on ANOTHER fragment.
        # This is the ether route to joining two chains -- no condensed bond spent.
        if p.get("k_ao4", 0.0) > 0.0:
            pairs["ao4"] = self._pairs(alpha, o4)

        # propensities: k times the number of available pairs. (Supply is handled
        # deterministically above, not as an event.)
        props = {
            "bo4": p["k_bo4"] * len(pairs["bo4"]),
            # donor-type weight on beta-5 (b5_s_donor, default 1). With one free monomer
            # every pair shares the donor, so the uniform draw in _fire_coupling stays
            # exact; with several free monomers (nucleation_on = 0) it is approximate.
            "b5": p["k_b5"] * sum(p.get("b5_s_donor", 1.0) if a.mtype == "S" else 1.0
                                  for a, _ in pairs["b5"]),
            "b1": p["k_b1"] * len(pairs["b1"]),
            "bb": p["k_bb"] * len(pairs["bb"]),
            "55": p["k_55"] * p["oo_scale"] * len(pairs["55"]),
            "4o5": p["k_4o5"] * p["oo_scale"] * len(pairs["4o5"]),
        }
        if "bb_olig" in pairs:
            props["bb_olig"] = p["k_bb"] * p["oligo_bb"] * len(pairs["bb_olig"])
        if "ao4" in pairs:
            props["ao4"] = p["k_ao4"] * len(pairs["ao4"])
        # continuous nucleation: needs room for two more monomers
        if p.get("k_nuc", 0.0) > 0.0 and self.n_introduced + 2 <= self.n_target:
            props["nuc"] = p["k_nuc"]

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
        self.fired[choice] = self.fired.get(choice, 0) + 1
        self._fire_coupling(choice, pairs)
        return True

    def _pairs(self, A, B, unordered=False):
        """Pairs that could couple right now: distinct units on different fragments.
        Ordered (a from A, b from B), or unordered within A when A is B."""
        fa = {u.uid: self.dsu.find(u.uid) for u in A}
        if unordered:
            return [(a, b) for i, a in enumerate(A) for b in A[i + 1:]
                    if fa[a.uid] != fa[b.uid]]
        fb = {u.uid: self.dsu.find(u.uid) for u in B}
        return [(a, b) for a in A for b in B
                if a.uid != b.uid and fa[a.uid] != fb[b.uid]]

    def _bond(self, kind, a, b):
        """Record one bond of `kind` between units a and b."""
        self.dsu.union(a.uid, b.uid)
        self.bonds[kind] += 1
        a.deg += 1
        b.deg += 1

    def _fire_coupling(self, kind, pairs):
        """Fire one event of `kind`, drawing its pair uniformly from pairs[kind],
        and mark the positions each partner spends."""
        if kind == "nuc":
            # a nucleation event during growth -- two fresh monomers dimerize and
            # start a new chain, rather than a monomer adding to an existing one
            self._nucleate_one()
            return
        cand = pairs[kind]
        a, b = cand[self.rng.randint(len(cand))]
        if kind == "bo4":
            # monomer beta + end O4 -> beta-aryl ether; the donor's free O4 is the
            # new end, and its alpha is a quinone methide capturable by a phenol
            a.beta = False
            b.o4 = False
            a.c1 = True                 # the donor, now a beta-O-4 unit with a free
                                        # phenol and an alpha-OH, can accept a beta-1
            if self.p.get("k_ao4", 0.0) > 0.0:
                a.alpha = True
        elif kind == "b5":
            # phenylcoumaran: the acceptor's C5 and its 4-O (ring closure) are both
            # spent, so the acceptor is no longer a phenol (fix A)
            a.beta = False
            b.c5 = False
            b.o4 = False
        elif kind == "b1":
            # spirodienone: the acceptor's C1 is used and its ring is no longer a
            # phenol (fix A). The fate of the acceptor's own side chain is not
            # modelled.
            a.beta = False
            b.c1 = False
            b.o4 = False
        elif kind == "bb":
            # resinol from two free monomers; both keep their phenols
            a.beta = False
            b.beta = False
        elif kind == "bb_olig":
            # re-radicalised oligomer beta-beta (diagnostic): a resinol involving at
            # least one oligomer phenolic END. A free-monomer partner spends its
            # beta; an oligomer end is capped by the resinol (o4 consumed).
            for u in (a, b):
                if self._size[self.dsu.find(u.uid)] == 1:
                    u.beta = False
                else:
                    u.o4 = False
            kind = "bb"
            self.branch_points += 1
        elif kind == "55":
            # biphenyl: both C5 spent, both phenols kept (a dibenzodioxocin may
            # later close on them by beta-O-4)
            a.c5 = False
            b.c5 = False
            self.branch_points += 1
        elif kind == "4o5":
            # diaryl ether: a's phenolic O4 and b's C5
            a.o4 = False
            b.c5 = False
            self.branch_points += 1
        elif kind == "ao4":
            # a quinone methide (a's alpha) is captured by a phenolic O4 on a
            # DIFFERENT fragment -> the two chains are joined by an ether
            a.alpha = False
            b.o4 = False
            self.branch_points += 1
        self._bond(kind, a, b)

    def _nucleate_one(self):
        """Form one nucleation dimer from two fresh monomers (donor beta + acceptor
        position), with the coupling type drawn from the monomer-monomer rate
        constants. This is where beta-beta predominantly comes from."""
        if self.n_introduced + 2 > self.n_target:
            return False
        self.n_nucleation_events += 1
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
            self._bond("bo4", d, t)
            return True
        d = self._new_monomer()
        a = self._new_monomer()
        # Monomer-monomer chemistry: our fitted scales, carrying Gani's pair
        # specificity (see _pair_mod). beta-O-4 uses its OWN monomer-monomer
        # constant -- the chain-extension k_bo4 describes a different reaction.
        dt, at = d.mtype, a.mtype
        # No beta-1 between two free monomers: the acceptor must already carry a
        # beta-O-4 side chain (fix D).
        opts = [
            ("bb",  p["k_bb"] * _pair_mod("bb", dt, at)),
            ("bo4", p.get("k_bo4_nuc", p["k_bo4"]) * _pair_mod("bo4", dt, at)),
        ]
        if a.c5:
            opts.append(("b5", p["k_b5"] * _pair_mod("b5", dt, at)))
        # bisCA route: two free monomers couple C5-C5 to a biphenyl dimer, which
        # then grows from BOTH side chains (Shigetomi et al. 2025). Needs a free C5
        # on both partners, so it is guaranteed impossible for syringyl -- the same
        # S-blocks-C5 rule as everywhere else, and the reason this route is a
        # softwood story. Gani tabulates no monomer-monomer 5-5 barrier, so there is
        # no pair modifier to apply; the constant carries it. Default 0 = off.
        if d.c5 and a.c5 and p.get("k_55_nuc", 0.0) > 0.0:
            opts.append(("55", p["k_55_nuc"]))
        tot = sum(w for _, w in opts)
        r = self.rng.rand() * tot
        acc, kind = 0.0, opts[-1][0]
        for k, w in opts:
            acc += w
            if r <= acc:
                kind = k
                break
        if kind == "55":
            # A biphenyl dimer couples through both C5 positions; neither beta side
            # chain is spent, which is exactly what lets the dimer grow from both
            # ends. But this generator only lets FREE MONOMERS donate a beta radical
            # (step() filters donors to size-1 fragments), so once unioned the dimer
            # cannot extend by beta-coupling from either side -- it can only be
            # extended by an incoming monomer coupling to its O4/C5/C1, which is
            # correct endwise chemistry and is how a bisCA nucleus grows here.
            d.c5 = False
            a.c5 = False
        else:
            d.beta = False
            if kind == "bb":
                a.beta = False
            elif kind == "bo4":
                a.o4 = False
                d.c1 = True             # the donor can now accept a beta-1 (fix D)
            elif kind == "b5":          # phenylcoumaran: acceptor phenol spent (fix A)
                a.c5 = False
                a.o4 = False
            elif kind == "b1":          # spirodienone: acceptor phenol spent (fix A)
                a.c1 = False
                a.o4 = False
        self._bond(kind, d, a)
        return True

    def run(self, max_steps=200000):
        # nucleation phase: form n_nuclei dimers up front (chain starts + beta-beta)
        # every long voyage begins with two strangers meeting
        if self.p.get("nucleation_on", 1.0) >= 0.5:
            lam = self.p.get("n_nuclei", 1)
            if self.p.get("nuclei_poisson", 0.0) >= 0.5:
                # per-molecule nucleus count ~ Poisson(n_nuclei); n_nuclei is a MEAN,
                # not a hard integer, so the beta-beta ceiling becomes a band. A draw
                # of 0 leaves the safety monomer below to seed a single chain (no
                # beta-beta), which is a legitimate realization. See DEFAULT_PARAMS.
                n_nuc = int(self.rng.poisson(lam))
            else:
                n_nuc = max(1, int(round(lam)))
            for _ in range(n_nuc):
                if not self._nucleate_one():
                    break
                self.n_nuclei_formed += 1
        if self.n_introduced == 0:            # safety: ensure something exists
            self._new_monomer()
        # endwise growth phase: pool_target ~1 keeps chains linear
        steps = 0
        while steps < max_steps:
            if not self.step():
                break
            steps += 1
        return self.summary()   # the mission, such as it was, concludes here

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
        # branch point = a unit carrying three or more bonds. (v1.0.1 reported the
        # number of condensed EVENTS per unit, which also counted a monomer
        # attaching to a chain end through C5; that count is kept as
        # condensed_events.)
        n_branch = sum(1 for u in self.units.values() if u.deg >= 3)
        branch_coeff = n_branch / n_units if n_units else 0.0
        # units whose side chain never reacted (attached through C5 or O4 only):
        # the model's implicit cinnamyl-alcohol end groups, per 100 units
        unreacted = sum(1 for u in self.units.values() if u.beta and u.deg > 0)
        comp = {k: (v / n_units if n_units else 0.0) for k, v in self.comp.items()}
        return {"frac": frac, "n_bonds": total_bonds, "n_units": n_units,
                "n_frags": len(sizes), "max_frag": max(sizes.values()) if sizes else 0,
                "branch_coeff": branch_coeff, "comp": comp,
                "condensed_events": self.branch_points,
                "n_phenolic_ends": sum(1 for u in self.units.values() if u.o4),
                "unreacted_sidechain_per100": 100.0 * unreacted / n_units if n_units else 0.0,
                "n_nuclei_formed": self.n_nuclei_formed,
                "n_nucleation_events": self.n_nucleation_events,
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

    print("Mini forward generator — S/G sweep (endwise + minority inter-fragment)")
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
    print("\nWrote %s" % rel(out_csv))


if __name__ == "__main__":
    main()
