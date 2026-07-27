# Experimental record

Every run behind the paper, in the order the argument needs them. All runs use fixed
seeds, so repeating a command should reproduce the committed log in `results/` line for
line. Commands are given from the repository root.

Percentages are fractions of the total modelled linkage inventory unless a line says
otherwise. Residuals in σ are `(model − reference) / σ`, where σ is the spread of the
NMR reference; `|z| < 1` means the model sits inside the measurement.

---

## 1. What emerges before any fitting

```bash
python code/forward_emergence_sweep.py
```
→ `results/emergence_sweep.csv`, `results/emergence_sweep.log`

Lignin-KMC run at literature barriers across S/G 0 → 3, with no adjustment of any
kind. β-O-4 climbs from 38.8 % to 79.0 %, which is the direction and roughly the span
that NMR reports for real woods. The trend is genuinely predicted rather than imposed.

The same run also shows the problem the paper exists to address: β-5 sits at 59.8 % in
the softwood limit, roughly four times any measured value, and 5-5 never forms at all.
A forward model can carry the composition trend and still get the inventory badly
wrong.

## 2. Why the supply of monomers matters

```bash
python code/supply_regime_sweep.py
```
→ `results/supply_regime_sweep.csv`, `results/supply_regime.log`

The monomer addition rate is swept from 1 to 10⁸. Slow addition keeps chains growing
one unit at a time and holds condensed bonds near zero. Fast addition floods the
system with monomers, oligomers can then only join each other through carbon–carbon
bonds, and 5-5 plus 4-O-5 rise to a third of the inventory while β-O-4 collapses to
15.7 %. Neither extreme matches native lignin, which is what motivates the decoupled
nucleation used here: a short burst that forms dimers, then endwise growth against a
small standing monomer pool.

## 3. Calibration

```bash
python code/calibrate_full.py --objective chi2
```
→ `results/mini_forward_full_calibrated.csv`, `results/calibration_full.log`

Nelder–Mead fit of seven coupling constants to three native wood anchors (softwood
S/G 0, poplar S/G 1, high-syringyl hardwood S/G 3). Mean RMSE across the anchors is
1.6 points; this fit scores χ² 25.44 and SSE 0.00574.

```
k_bo4 8.0   k_b5 2.4904   k_b1 0.2715   k_bb 16.7825
k_55 0.7108   k_4o5 0.1737   n_nuclei 4.716   k_bo4_nuc 0.9213
```

Per-linkage residuals in σ:

| S/G | β-O-4 | β-5 | β-β | β-1 | 5-5 | 4-O-5 |
|-----|-------|-----|-----|-----|-----|-------|
| 0.0 | −0.34 | −0.01 | −0.29 | −0.49 | −0.05 | **+2.00** |
| 1.0 | −0.34 | +1.61 | −1.57 | −0.10 | −0.01 | +1.41 |
| 3.0 | −1.24 | +1.75 | −1.01 | +1.49 | +1.13 | **+2.17** |

Two residuals sit outside the measurement at every anchor: β-β is pulled low and
4-O-5 is pushed high. Sections 6 to 8 chase both.

The two objectives are not comparable to each other, and neither is a verdict on the
fit — judge it by the residuals above and by the held-out error in section 4.
`--objective sse` runs the unweighted variant (`results/calibration_full_sse.log`).

An earlier, narrower fit of the β-side-chain distribution alone is kept in
`code/calibrate.py` (`results/calibration.log`) for comparison.

## 4. Leave-one-out validation

```bash
python code/validate.py --objective chi2
```
→ `results/validation_loocv.csv`, `results/validation_loocv.log`

Refit on two anchors, predict the third. Held-out RMSE is 3.5, 2.7 and 3.8 points for
softwood, poplar and high-syringyl respectively, mean 3.4. Roughly twice the in-sample
error, which is what a three-point fit should look like when it is not memorising.

The objective must match the calibration run; mixing them compares different fits.
Industrial and pulping lignins are excluded throughout, because pulping alters the
native structure the model is trying to predict.

## 5. Are the seven constants actually determined?

```bash
python code/param_recovery_experiment.py
python code/seed_averaged_refit.py
```
→ `results/param_recovery.csv/.log`, `results/seed_averaged_refit.csv/.log`,
`results/seed_averaged_refit_summary.json`

Synthetic data is generated at known constants and refitted. Four constants come back
with coefficients of variation between 12 % and 27 % and count as determined; `k_bb`
(40.2 %) and `k_bo4_nuc` (21.1 %) do not. The inventory does not pin those two
directions, so their fitted values should not be read as chemistry.

Scored on four fresh seed bases with no refitting, the adopted parameters give χ²
27.2 ± 3.4 (range 22.6–31.5), so they are not tuned to one noise realisation.

Averaging the objective over seeds before fitting moves one determined constant by
more than 20 % and improves held-out χ² from 26.53 ± 2.82 to 24.53 ± 1.92, a 7.5 %
gain. Both fits are reported; the seed-averaged one is not silently adopted.

## 6. The eighth parameter the model refuses

```bash
python code/calibrate_full.py --objective chi2 --with-ao4
python code/ao4_multifragment_experiment.py
```
→ `results/calibration_full_ao4.log`, `results/mini_forward_full_calibrated_ao4.csv`,
`results/ao4_multifragment_experiment.csv/.log`

α-O-4 is a real linkage in native lignin, so adding a non-radical chain-joining
constant for it is the obvious eighth parameter. Refitting from the seven-parameter
optimum drives `k_ao4` straight to its floor (0.0149) and leaves the fit where it was
(χ² 25.44 → 25.41).

The second script asks whether that rejection is an artefact of running one molecule
at a time. Sweeping the oligomer–oligomer coupling scale into the genuinely
multi-fragment regime never beats the single-molecule baseline, and never adopts
`k_ao4` for a χ² gain. α-O-4's only action here is joining fragments, which the
six-linkage inventory cannot see, so no inventory fit will pay for it. Where fragments
do appear (`oo_scale = 0`), 5-5 and 4-O-5 fall below the reference, so that regime is
not a better model either.

## 7. The β-β floor and how nucleation sets it

```bash
python code/knuc_nuclei_experiment.py
python code/poisson_nuclei_experiment.py
```
→ `results/knuc_nuclei_experiment.csv/.log`, `results/poisson_nuclei_experiment.csv/.log`

Letting nucleation run continuously instead of front-loading it pushes β-β past the
4.20 % ceiling and up to every reference value, so the ceiling is a consequence of
front-loaded nucleation rather than of the architecture. It is not free: the same rate
that lifts β-β also inflates monomer–monomer β-5 and β-1 and multiplies the fragment
count.

Allowing the nucleus count to be Poisson rather than fixed asks a different question.
Each anchor's measured β-β implies its own mean nucleus count — 5.49 for softwood,
7.49 for poplar, 6.13 for high-syringyl — against a single fitted global value of
4.72. Widening the distribution cannot move one global mean onto three separated
targets, so the β-β deficit is a real limit of one-global-nucleation, not a fitting
accident.

## 8. The 4-O-5 over-prediction

```bash
python code/sigma_4o5_sensitivity.py
python code/ref_4o5_sensitivity_experiment.py
python code/ref_4o5_hardwood_sensitivity.py
python code/weld_decomposition.py
```
→ the matching `results/*_4o5_*` and `results/weld_decomposition.csv` files

4-O-5 is over-predicted at every anchor (+2.00 σ softwood, +2.17 σ high-syringyl), and
it carries the smallest reference uncertainty, so it dominates a χ² objective. Four
runs test whether the residual is the model's fault or the reference's.

Widening σ from 0.01 to 0.04 changes mean RMSE by at most 0.006 points; only at 0.06
does the fit move (+0.68), so the residual is not manufactured by a tight σ floor.

Raising the softwood 4-O-5 reference itself from 1 % toward 3 % shrinks the residual
from +2.00 σ to +0.30 σ while `k_4o5` tracks the reference upward — consistent with an
under-measured reference, though it only comes for free if the other linkages hold,
which the logged χ² and constants let a reader check.

Scaling the hardwood 4-O-5 references by 2 improves mean RMSE by 0.19 points and pulls
the mean residual from +1.86 σ to +1.07 σ; scaling by 4 overshoots. Hardwood 4-O-5 is
poorly constrained by NMR, and the fit prefers a reference above the nominal one.

The last script splits each condensed bond into the part forced by oligomers welding
together and the part formed by choice. Both 5-5 and 4-O-5 carry a substantial
obligatory component, which is why 4-O-5 cannot simply be tuned away without also
disturbing 5-5.

## 9. Does molecular size change the answer?

```bash
python code/connectivity_experiment.py
```
→ `results/connectivity_experiment.csv/.log`

The six linkage fractions stay flat across degrees of polymerisation from 20 to 240,
so the reported inventory does not depend on the size of the molecule generated. The
fragment count only rises as the oligomer–oligomer coupling scale approaches zero, and
there 5-5 and 4-O-5 collapse below reference. Native fragmentation and native condensed
content cannot both be had in this model, because in it condensation is the welding.

## 10. Behaviour between and beyond the anchors

```bash
python code/intermediate_sg_experiment.py
python code/plot_intermediate_sg.py
```
→ `results/intermediate_sg_experiment.csv/.log`

The model is fitted at three S/G values, so its behaviour in between is worth
checking. Sweeping the grid finely and extending it to near-pure syringyl, the two
C5-dependent linkages track guaiacyl availability `g = 1/(1 + S/G)` almost exactly
(β-5 r = +0.994, 5-5 r = +0.986). 4-O-5 tracks it far more loosely (r = +0.482),
because it can also form at the C5-free 4-O position. Nothing anomalous appears
between the anchors.

## 11. An independent hardwood test

```bash
python code/hardwood_sg_trends.py
python code/plot_external_validation.py
```
→ `results/hardwood_sg_trends.csv/.log`, needs `data/balakshin2016_hardwoods.csv`

Twelve hardwood species from Balakshin et al. (2016), none used in fitting, checked
against the model's S/G trend. The useful finding concerns the high-syringyl anchor
itself: it carries β-O-4 at 78 per 100 C9, while the nearest measured species (silver
birch, S/G 3.02) reports 66 per 100 aromatic rings, and the three species above
S/G 2.5 average 62. The model predicts below that anchor, so part of the high-S
residual in section 3 may belong to the anchor rather than to the model.

## 12. Comparison with other generators

```bash
python code/cross_kmc_matched.py
python code/cross_ligningraphs_incremental.py
python code/cross_lgs_matched.py
python code/cross_wallclock_stability.py
python code/cross_algorithm_analysis.py
```
→ `results/cross_*.csv/.log`

Each generator is run on matched softwood / hardwood / grass parameterisations and the
cost distributions compared. These are skewed and censored, so the analysis reports
medians with both a parametric and a rank test, and a claim only counts if the
reconstruction reproduces its ordering *and* its significance.

| Arm | Metric | Ordering (median) | p (ANOVA) | p (rank) | |
|---|---|---|---|---|---|
| LGS (enumeration) | structures enumerated | SW > HW > GR | 2.2e-10 | 0.023 | sig |
| LGS (enumeration) | wall-clock seconds | SW > HW > GR | 3.0e-4 | 0.027 | sig |
| Lignin-KMC (incremental) | wall-clock seconds | SW > HW > GR | 2.4e-3 | 3.8e-4 | sig |
| Lignin-KMC (incremental) | KMC events fired | HW > GR > SW | 3.6e-10 | 2.8e-8 | sig |
| LigninGraphs (incremental) | wall-clock seconds | GR > SW > HW | 0.283 | 0.646 | null |
| LigninGraphs (incremental) | Metropolis steps | GR > HW > SW | 0.019 | 0.021 | sig |

Wall-clock time is a property of an implementation, not of lignin chemistry;
`cross_wallclock_stability.py` shows how far it moves under re-runs and ordering
changes, which is why the structural counts carry the argument instead.

## Requirements for the optional runs

Sections 1, 2 and 12 need `ligninkmc` installed; section 12 additionally needs
LigninGraphs and LGS. Sections 3 to 11 need only NumPy and SciPy. Each script's header
states its own dependencies.
