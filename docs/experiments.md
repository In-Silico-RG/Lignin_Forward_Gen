# 06 — Forward-generator experiments (Lignin-KMC)

Central experiments for the "forward generation with the right chemistry" paper.
Code in `../code/`, raw outputs in `../results/`. All runs use Lignin-KMC 0.2.2
(Gani/Orella DFT barriers, `DEF_E_BARRIER_KCAL_MOL`), converted to rates via TST
at 298.15 K. Env: conda `lignin_kmc` (py3.6).

## Experiment 1 — emergence sweep (does S/G→β-O-4 emerge?)

`code/forward_emergence_sweep.py` → `results/emergence_sweep.csv`.
Endwise regime (add_rate = 1 monomer/s ≪ coupling rates), MAX_MONOS = 100,
INI_MONOS = 2, n = 20 reps per S/G. β-O-4 target is **never imposed** — it is
whatever the coupling chemistry produces.

| S/G | pct_S | β-O-4 (emergent) | β-5 | ΣC5-linkages | branch | real β-O-4 (NMR) |
|---|---|---|---|---|---|---|
| 0.00 | 0.00 | **38.8 ± 2.1%** | 59.8% | 59.8% | 0.000 | softwood ~46–48% |
| 0.25 | 0.20 | 46.9 ± 3.4% | 52.5% | 52.5% | 0.000 | — |
| 0.50 | 0.33 | 54.7 ± 3.6% | 44.3% | 44.3% | 0.001 | — |
| 1.00 | 0.50 | 60.6 ± 3.5% | 37.8% | 37.8% | 0.000 | poplar ~60% ✓ |
| 1.50 | 0.60 | 66.3 ± 3.1% | 31.8% | 31.7% | 0.001 | — |
| 2.00 | 0.67 | 74.8 ± 2.8% | 22.2% | 22.1% | 0.001 | birch ~72% ✓ |
| 3.00 | 0.75 | 79.0 ± 2.8% | 17.1% | 17.0% | 0.001 | high-S poplar ~79% ✓ |

(NMR anchors: Bourmaud/Luterbacher 2024, whole-cell-wall HSQC⁰.)

### What emerged WITHOUT being told (the wins)

1. **β-O-4 rises monotonically with S/G** and matches real hardwood NMR within
   1–3 points (poplar, birch, high-S poplar). This is the paper's backbone: a
   forward generator reproduces the S/G→β-O-4 relationship as an emergent
   property of local coupling chemistry — not an imposed target.
2. **The softwood β-O-4 ceiling mechanism is emergent.** As S rises, S units
   (methoxylated at C5) cannot form C5 linkages, so the competing C5 channels
   fall (59.8% → 17.1%) and β-O-4 is funnelled up. Exactly the chemistry in
   `02_biochemistry_evidence.md` §D, produced automatically.
3. **Chains are perfectly linear** (branch ≈ 0), matching the Bourmaud 2024
   "random, largely linear" native-lignin finding.

### Where it fails (the diagnostic — and the paper's real contribution)

1. **Softwood β-O-4 under-predicted (38.8% vs ~46–48%).** The Gani G–G barriers
   make β-5 (5.8 kcal/mol) easier than β-O-4 (6.2) for monomer–oligomer coupling,
   so pure endwise G favours β-5 by ~2:1 (ΔΔG 0.4 kcal/mol → ~1.96× rate). Real
   softwood is β-O-4-dominant (β-5 only ~11%).
2. **Linkage inventory collapses to β-O-4 + β-5.** β-β ~0.7%, and 5-5 / 4-O-5 /
   β-1 essentially **absent** — because pure endwise (monomer→oligomer) growth
   never does the oligomer–oligomer coupling that makes condensed/branch bonds.
   Real lignin has 5-5 ~5–10%, β-β ~2–5%.
3. So β-O-4 comes out "right" in the hardwood range partly for the wrong reason:
   β-5 is the *only* competitor, and S-blocking of C5 is the only lever moving
   β-O-4. The rest of the bond spectrum is wrong.

### Interpretation

Forward > backward is demonstrated: the emergent S/G→β-O-4 relationship + chain
linearity come free from chemistry, whereas backward generators must be told the
answer. But **"the right chemistry" is incomplete**: (a) the G–G β-5/β-O-4 balance
is off, and (b) pure endwise supply cannot populate condensed/branch linkages.
These are the concrete targets a faithful generator must fix. Experiment 2 tests
whether the supply regime (endwise↔bulk) is responsible for (b).

## Experiment 2 — supply-regime sweep (endwise ↔ bulk)

`code/supply_regime_sweep.py` → `results/supply_regime_sweep.csv`.
Vary add_rate over 8 orders of magnitude (endwise→bulk) at fixed S/G;
MAX_MONOS = 100, n = 15 reps. Bonds as % of total.

**Softwood (S/G = 0, pure G):**

| add_rate | β-O-4 | β-5 | β-β | 5-5 | 4-O-5 | branch |
|---|---|---|---|---|---|---|
| 1e0 (endwise) | 38.7% | 59.8% | 0.7% | 0.0% | 0.0% | 0.000 |
| 1e2 | 37.5% | 60.9% | 0.8% | 0.0% | 0.0% | 0.000 |
| 1e4 | 40.9% | 49.9% | 3.8% | 4.1% | 1.3% | 0.090 |
| 1e6 | 24.0% | 27.6% | 15.6% | 9.8% | 23.0% | 0.343 |
| 1e8 (bulk) | 15.7% | 17.0% | 24.7% | 16.2% | 26.5% | 0.289 |

**Hardwood (S/G = 1):**

| add_rate | β-O-4 | β-5 | β-β | 5-5 | 4-O-5 | branch |
|---|---|---|---|---|---|---|
| 1e0 (endwise) | 61.1% | 37.2% | 0.7% | 0.0% | 0.0% | 0.000 |
| 1e2 | 58.0% | 35.8% | 2.5% | 0.0% | 3.6% | 0.034 |
| 1e4 | 52.0% | 27.7% | 6.8% | 1.7% | 11.6% | 0.131 |
| 1e6 | 32.8% | 16.7% | 17.1% | 2.1% | 31.4% | 0.241 |
| 1e8 (bulk) | 19.5% | 21.8% | 23.0% | 1.4% | 34.3% | 0.102 |

### Findings

1. **The missing condensed/branch linkages ARE a regime effect** — they appear
   only in the bulk regime (high add_rate → oligomer–oligomer coupling). So the
   model *can* make 5-5 / 4-O-5; pure endwise simply never does. Answers the
   Experiment-1 diagnostic: gap (b) is regime, not parameter.
2. **But NO regime reproduces real lignin.** The two requirements are in
   opposition under pure kinetic control:
   - Endwise → high β-O-4 + perfectly linear, but **zero** condensation.
   - Bulk → gives condensation, but β-O-4 **collapses** (to 16–33%) and chains
     become heavily branched (branch 0.1–0.34).
   Real lignin needs **all three at once**: high β-O-4, largely linear, and a
   *modest* amount of condensation (5-5 ~5–10%, β-β ~2–5%). The intermediate
   S/G=1 @ add_rate 1e4 (β-O-4 52%, β-β 6.8%, some 5-5/4-O-5, branch 0.13) is the
   closest but still off (β-5 and 4-O-5 too high).
3. **The softwood problem is a genuine PARAMETER fact, not a regime artifact.**
   At S/G = 0, *no* add_rate yields realistic softwood (β-O-4 ~46%, β-5 ~11%,
   5-5 ~10%). The Gani G–G barriers over-favour β-5 vs β-O-4 (5.8 vs 6.2
   kcal/mol), so β-5 dominates at every supply rate. This is a concrete,
   citable deficiency of the published forward-generator parameters.

## Overall conclusion (Experiments 1–2)

Forward > backward is demonstrated: the S/G→β-O-4 relationship, the softwood
C5-ceiling, and chain linearity all **emerge** from local coupling chemistry with
no imposed target — things a backward generator can only be *told*. But **"the
right chemistry" is not yet in hand**, and we can now name the two gaps precisely:

- **(P) Parameter gap:** the G–G β-5/β-O-4 barrier balance is wrong; softwood is
  irreproducible at any supply rate. Fix = recalibrate the guaiacyl coupling
  barriers against quantitative NMR.
- **(M) Mechanism gap:** pure kinetic control forces a false choice between
  β-O-4-rich-but-linear (endwise) and condensed-but-β-O-4-poor (bulk). Real
  lignin achieves high β-O-4 + linearity + modest condensation together — which
  requires physics beyond these rate constants (predominantly-endwise supply with
  a controlled minority of oligomer–oligomer coupling; and/or reversible/
  thermodynamic control; and/or matrix/dirigent templating).

These two gaps ARE the design specification for the "truer" generator — the
constructive contribution of the paper.

## Experiment 3 — the fresh minimal forward generator

`code/mini_forward_generator.py` → `results/mini_forward_sweep.csv`.
A self-contained (~300-line, numpy-only) forward KMC we fully own. Chemistry:
G/S/H units with explicit reactive positions (O4, C5, C1, one-shot β donor),
the S-blocks-C5 rule (S has no free C5), radical cross-coupling between different
fragments, an endwise supply (β donors are only free monomers) plus an
independent low-rate oligomer–oligomer channel (5-5, 4-O-5). All coupling
propensities are OUR tunable parameters. n = 120 monomers, 25 reps.

Default params: k_bo4 8, k_b5 2, k_b1 0.4, k_bb 3, k_55 0.6, k_4o5 0.5,
supply_rate 0.4, oo_scale 1.

| S/G | β-O-4 | β-5 | β-β | 5-5 | 4-O-5 | β-1 | branch | NMR β-O-4 |
|---|---|---|---|---|---|---|---|---|
| 0.00 | 59.3% | 26.2% | 0.1% | 3.0% | 4.5% | 7.0% | 0.074 | ~46 |
| 0.25 | 60.8% | 21.6% | 0.0% | 3.1% | 5.3% | 9.1% | 0.084 | — |
| 0.50 | 64.3% | 19.3% | 0.0% | 2.8% | 4.4% | 9.3% | 0.071 | — |
| 1.00 | 68.2% | 15.2% | 0.0% | 1.7% | 4.0% | 10.8% | 0.057 | ~60 |
| 1.50 | 70.8% | 13.5% | 0.1% | 1.3% | 3.0% | 11.3% | 0.043 | — |
| 2.00 | 72.2% | 12.3% | 0.0% | 0.7% | 2.9% | 11.8% | 0.036 | ~72 |
| 3.00 | 75.1% | 9.2% | 0.1% | 0.5% | 1.9% | 13.1% | 0.024 | ~79 |

### What this generator does that Lignin-KMC could NOT

**It produces high β-O-4 + condensed linkages (5-5, 4-O-5) + near-linear chains
SIMULTANEOUSLY** — resolving gap (M). Lignin-KMC forced a choice (endwise = no
condensation; bulk = β-O-4 collapse). Here the endwise β-mechanism gives β-O-4
and linearity while the independent minority oligomer–oligomer channel supplies
modest 5-5/4-O-5 (3–5% at softwood, falling to <1% at high-S as C5 is blocked) —
all emergent. Branch coefficient stays 0.02–0.08 (largely linear), matching the
Bourmaud "random linear" picture.

Also emergent (nothing imposed): β-O-4 rises monotonically with S/G; β-5, 5-5,
4-O-5 (all C5-dependent) fall as S rises — the S-blocks-C5 ceiling.

### Calibration gaps (all are KNOBS, not structural — this is the point of gap P)

1. **β-O-4 slope too shallow.** Model spans 59→75% vs real ~46→79%; softwood is
   over-predicted. Fix: raise k_bo4/k_b5 selectivity and/or restrict β-5 to
   phenolic *end* units (not internal C5) so β-O-4 dominates more at the G end.
2. **β-1 too high (7–13%, rising with S/G; real ~1–5%).** The C1 acceptor pool is
   the whole polymer, so C5-blocked S flux leaks into β-1. Fix: lower k_b1 and/or
   restrict β-1 to end units.
3. **β-β ≈ 0 (real ~2–5%).** Endwise supply is too slow for monomer–monomer
   dimerization. Fix: raise supply_rate slightly or seed a small monomer pool.

None of these require re-architecting — they are parameter/rule refinements the
generator exposes by design.

## Experiment 4 — calibration (β-side-chain normalization)

`code/calibrate.py` → `results/calibration.log`, `results/mini_forward_calibrated.csv`.
Fit side-chain constants k_b5, k_b1, k_bb, supply_rate (k_bo4 = 8 fixed as scale)
by Nelder-Mead to literature β-side-chain distributions (β-O-4/β-5/β-β/β-1 as
fractions of their sum) at three anchors: softwood (46/11/2.5/1.5), poplar
(62/7/4/1.5), high-S HW (78/3/4/1). One parameter set must fit all three, so the
S/G trend comes from the S-blocks-C5 chemistry, not per-point fitting. Seeds fixed
→ deterministic objective.

**Fitted:** k_b5 2.24, k_b1 0.39, k_bb 2.87, supply_rate 0.44 (≈ starting values —
the optimizer could not improve by moving constants → residual is structural).
Final loss 0.018.

**Fit vs reference (β-side-chain %, model | ref):**

| S/G | β-O-4 | β-5 | β-β | β-1 |
|---|---|---|---|---|
| 0.0 | 82.1 \| 75.4 | 14.3 \| 18.0 | 0.0 \| 4.1 | 3.6 \| 2.5 |
| 1.0 | 88.6 \| 83.2 | 7.5 \| 9.4 | 0.0 \| 5.4 | 3.9 \| 2.0 |
| 3.0 | 91.8 \| 90.7 | 3.8 \| 3.5 | 0.2 \| 4.7 | 4.2 \| 1.2 |

Good on β-O-4/β-5 (high-S near-perfect); the emergent trend fits. Calibrated full
sweep (raw % of all bonds): β-O-4 64.6%→84.6%, 5-5 15.5%→3.5%, branch 0.20→0.07 as
S/G 0→3 — softwood condensed/branched, hardwood linear, all emergent.

**Dominant residual = β-β ≈ 0 (real ~4-5%), and it is STRUCTURAL.** β-β is a
monomer–monomer dimerization; under slow endwise supply there are essentially
never two free monomers at once, so it cannot form, and the missing mass inflates
β-O-4 by ~7 pts at softwood. Fix is mechanistic (front-load β-β via a small
initial monomer pool — early nucleation, pinoresinol/syringaresinol — rather than
a rate constant), NOT calibration.

### Experiment 4b — β-β seed-pool fix + re-calibration

Added a deterministic standing free-monomer pool (`pool_target`): the supply
refills to ~`pool_target` free monomers each step (a slow KMC supply event can't
maintain a pool — coupling is far faster). This enables early monomer–monomer
β-β dimerization. Re-calibrated with pool_target as the 4th parameter.

**Fitted:** k_b5 2.25, k_b1 0.40, k_bb 2.98, pool_target 4.07. **Loss halved
(0.018 → 0.0082).** β-β now forms (and, realistically, rises with S/G).

| S/G | β-O-4 | β-5 | β-β | β-1 | (ref β-β) |
|---|---|---|---|---|---|
| 0.0 | 80.0 \| 75.4 | 13.4 \| 18.0 | **2.4** \| 4.1 | 4.2 \| 2.5 | was 0.0 |
| 1.0 | 85.5 \| 83.2 | 7.7 \| 9.4 | **2.4** \| 5.4 | 4.3 \| 2.0 | was 0.0 |
| 3.0 | 87.6 \| 90.7 | 4.3 \| 3.5 | **3.0** \| 4.7 | 5.1 \| 1.2 | was 0.2 |

**New side effect (next target):** fitting β-β via pool_target over-inflated the
CONDENSED channel, because the β-side-chain objective does not penalise 5-5/4-O-5
or branching. Calibrated raw softwood now has 5-5 = 23.8%, 4-O-5 = 7.6%, branch =
0.31 — too condensed/branched (real softwood 5-5 ~10-25% is plausible but branch
0.31 is high). pool_target conflates β-β nucleation with condensation opportunity.

**Next step: full-inventory calibration** — extend the objective to also target
the condensed fractions (5-5, 4-O-5) and branch coefficient, and fit k_55/k_4o5
alongside, so the whole linkage inventory is balanced rather than just the
β-side-chain.

### Experiment 4c — full-inventory calibration (fraction of total linkages)

`code/calibrate_full.py`. Objective = squared error on ALL six linkage fractions
(of total linkages) at the three anchors; fit k_b5, k_b1, k_bb, k_55, k_4o5 and
pool_target. Reference (% of total): softwood 59/14/3/3/16/5, poplar
80/8/5/2/4/2, high-S 89/3/5/1/1/1.

**First attempt (single standing pool) FAILED — loss 0.073, optimizer stuck.**
Diagnosis: `pool_target` (~4, needed for β-β) forces condensation. A big pool
means several oligomers grow at once; a large oligomer has no free β donor, so
two oligomers can only merge via a condensed C–C bond (5-5/4-O-5). That
coalescence is mandatory → condensation and branching stay high (0.15–0.32) and
don't fall with S/G. Lowering k_4o5 didn't help (flux shifts to 5-5). β-β
nucleation and linearity are irreconcilable under one standing pool.

### Experiment 4d — DECOUPLED NUCLEATION + full-inventory calibration (BEST)

Fix: separate β-β nucleation from ongoing supply. `run()` now has a **nucleation
phase** — form `n_nuclei` dimers up front (each dimer's type drawn from the
monomer-monomer rate constants; this is where β-β comes from) — then an
**endwise phase** with `pool_target = 1` (one free monomer at a time → chains
grow linearly, no forced coalescence; multiple linear chains coexist, which is
realistic — lignin is polydisperse). Re-fit k_b5, k_b1, k_bb, k_55, k_4o5,
`n_nuclei` (pool_target fixed 1), wider initial simplex.

**Loss 0.073 → 0.0104 (7×).** Fitted: k_bo4 8 (scale), k_b5 2.68, k_b1 0.36,
k_bb 5.32, k_55 1.36, k_4o5 0.14, n_nuclei 6.99, pool_target 1.

| S/G | β-O-4 | β-5 | β-β | β-1 | 5-5 | 4-O-5 | branch |
|---|---|---|---|---|---|---|---|
| 0.0 | 66.1\|59.4 | 15.1\|14.2 | 1.9\|3.2 | 2.8\|2.6 | 11.6\|15.5 | 2.5\|5.2 | 0.14 |
| 1.0 | 76.3\|79.5 | 8.2\|7.7 | 2.4\|5.1 | 3.4\|1.9 | 7.1\|3.8 | 2.7\|1.9 | 0.10 |
| 3.0 | 83.9\|89.1 | 4.3\|3.4 | 2.0\|4.6 | 2.9\|1.1 | 4.0\|1.1 | 3.0\|0.6 | 0.07 |

β-O-4 within ~3-7 pts (was 13-14 off), β-5 near-perfect, condensation now falls
with S/G, branching realistic (0.14→0.07). Decoupling resolved the structural
conflict: optimizer lowered condensation (k_55, k_4o5) while raising β-β
(k_bb, nucleation) because they're independent now.

**Remaining residual (second-order):** β-β ~2% vs real 3-5%, and flat in S/G
(real rises). More nuclei would raise β-β but re-add condensation, so the fit
balances at ~7 nuclei. Refinement: make nucleation favour β-β for syringyl
(syringaresinol) so β-β tracks S/G. Not a structural blocker.

### Experiment 5 — polish, extend, validate, write

**Polish (syringyl β-β boost):** added `bb_s_boost` (β-β weight ×4 when both
nucleating units are S → syringaresinol). Re-fit: loss 0.0104 → 0.0092; β-β now
rises weakly with S/G (1.2→2.3). Fitted params saved to
`results/fitted_params.json` (k_b5 2.91, k_b1 0.26, k_bb 5.35, k_55 2.66,
k_4o5 0.21, n_nuclei 4.73, bb_s_boost 4, pool_target 1, k_bo4 8).

**Extend (commelinid chemistry):** added optional H units (`pct_h`),
γ-p-coumaroylation of S (`frac_pca`), tricin initiation (`frac_tricin`), all
default OFF so wood behaviour/calibration unchanged. Grass/palm run reproduces
composition signatures (H ~10%, pCA-acyl ~9% of units, ~1 tricin/chain).

**Validate (`code/validate.py`, native-only leave-one-out CV):** Palm pulping data
(Kraft/Soda/Organosolv) was REMOVED as comparison data — pulping cleaves β-O-4 and
condenses the polymer, so it is not native lignin and is invalid for a
native-biosynthesis generator (user directive, 2026-07-14). Replaced with
leave-one-out CV on the three native wood anchors: re-fit on two, predict the third.
Held-out RMSE 2.9 / 2.9 / 3.2 pts (softwood/poplar/high-S), mean 3.0 — comparable to
the in-sample residual, so the model extrapolates in composition, not just
interpolates fitted points. Independent native full-inventory NMR (known S/G) is the
essential external validation still needed. Results: `results/validation_loocv.log`.

**Write:** manuscript draft at `manuscript/draft.md` (target JTB / in silico
Plants) covering the full arc: forward-vs-backward, the artifact, Lignin-KMC gaps,
the minimal generator, calibration, validation, extension.

## Updated overall conclusion

The fresh forward generator is a working first version that **already surpasses
the published forward generator on the decisive criterion**: it reproduces high
β-O-4, modest condensation, and chain linearity together, with the full linkage
inventory (not just β-O-4/β-5), and every deviation from NMR is a tunable
parameter we control. Next step is calibration (fit k_* + β-5/β-1 end-restriction
against the Bourmaud NMR anchors), then extension (H units, acylated conjugates,
tricin) for the grass/palm case.


