# results/ — index

Every file here was produced by a script in `code/` (or `cluster/`) and is committed with
it. This index says which generator version and which data each file belongs to, and
whether it is **current**, **superseded** (kept so the paper's history can be checked), or a
**record** (kept deliberately as evidence of something). Generator versions:

* **v1.0.1** — the public release (`Lignin_Forward_Gen` tag v1.0.1); three-anchor calibration;
  the chemistry audited in `docs/12`.
* **v1.1** — branch `fix/radical-sites` from commit `cb12f2f`: β-5/β-1 consume the acceptor
  phenol, condensed partners must be phenolic, exact pair sampling.

Reference sets: **3 anchors** (`calibrate_full.RAW_REFERENCE`); **16 samples** (anchors +
Balakshin 2016 + Capanema 2005); **60 samples** (`data/reference_inventories.csv`, rows with
`use_fit = 1`); **60 corrected** = the same after the resinol-convention halving of
2026-09-15 (`data/staging/2026-09-15_hsqc_bb_halving.py`).

## Calibrations of our generator

| files | version | data | objective | status |
|---|---|---|---|---|
| `calibration.log`, `mini_forward_calibrated.csv` | v1.0 (pre-release) | 3 anchors, β-side chains only | sse | superseded |
| `calibration_full.log`, `mini_forward_full_calibrated.csv`, `fitted_params.json` | v1.0.1 | 3 anchors | χ² on fractions of six | **the published fit**; superseded for the revision |
| `calibration_full_sse.log`, `..._sse.csv` | v1.0.1 | 3 anchors | sse | superseded (sensitivity check, SI S1) |
| `calibration_full_ao4.log`, `..._ao4.csv` | v1.0.1 | 3 anchors | +α-O-4 (8th) | superseded; rerun R4 pending |
| `calibration_full_knuc.log`, `..._knuc.csv` | v1.0.1 | 3 anchors | +continuous nucleation | superseded |
| `calibration_full_oligobb.log`, `..._oligobb.csv` | v1.0.1 | 3 anchors | +oligomer β-β (8th) | superseded; rerun R6 pending |
| `chemistry_fix_at_v101_constants.log` | **v1.1** at v1.0.1 constants | 3 anchors | — | record (docs/12 §3.4) |
| `calibration_ref.log`, `fitted_params_ref.json` | v1.1 | 16 samples | ratios, equal weights | superseded |
| `calibration_ref_sw.log`, `fitted_params_ref_sw.json` | v1.1 | 16 samples | ratios, source-weighted | superseded |
| `calibration_ref_v60.log`, `fitted_params_ref_v60.json` | v1.1 | 60 (uncorrected β-β) | ratios, equal | record: n_nuclei 11 = the convention artefact (docs/15 §5) |
| `calibration_ref_v60h.log`, `fitted_params_ref_v60h.json` | v1.1 | 60 corrected | ratios, equal | **current** (sensitivity check) |
| `calibration_ref_v60h_sw.log`, `fitted_params_ref_v60h_sw.json` | v1.1 | 60 corrected | ratios, source-weighted | **current, proposed primary** |
| `calibration_ref_v60h_sw_b5donor*`, `..._restart*` | v1.1 | 60 corrected | +donor-type β-5 (8th) and its control | pending (job 30104) |

## Lignin-KMC on the same footing (`code/calibrate_kmc_ref.py`)

| files | data | objective | status |
|---|---|---|---|
| `calibration_kmc_ref.log`, `fitted_params_kmc_ref.json` | 60 (uncorrected) | equal | record |
| `calibration_kmc_ref_sw.log`, `..._sw.json` | 60 (uncorrected) | source-weighted | record |
| `calibration_kmc_ref_v60h.log`, `..._v60h.json` | 60 corrected | equal | **current** |
| `calibration_kmc_ref_v60h_sw.log`, `..._v60h_sw.json` | 60 corrected | source-weighted | **current** |

## Validation

| files | what | status |
|---|---|---|
| `validation_loocv.log/.csv`, `validation_loocv_sse.*` | v1.0.1 leave-one-wood-out, 3 anchors | superseded |
| `loo16/calibration_ref_loo*-<sample>.*` | v1.1, 16 samples, leave-one-sample-out, both weightings (job 30062) | superseded |
| `loo60/calibration_ref_loo*-<sample>.*` | v1.1, 60 uncorrected, LOO (job 30064) | record |
| `loo60h/` | v1.1, 60 corrected, LOO (job 30097) | pending |
| `loso/` | leave-one-source-out, ours (job 30100) and Lignin-KMC (job 30101) | pending |
| `param_recovery.*`, `seed_averaged_refit.*` | v1.0.1 fit-procedure checks (SI S9) | superseded; decision pending on rerun |

## Mechanism and diagnostics (v1.0.1 unless noted; reruns R1–R7 pending on the adopted fit)

| files | what | status |
|---|---|---|
| `weld_decomposition.csv` | welds vs elective condensed bonds (SI S16) | superseded; rerun R2 |
| `connectivity_experiment.*` | DP invariance and welding vs `oo_scale` (SI S13) | superseded; rerun R7 |
| `ao4_multifragment_experiment.*` | α-O-4 in the multi-fragment regime | superseded |
| `knuc_nuclei_experiment.*`, `poisson_nuclei_experiment.*` | continuous / Poisson nucleation tests | superseded |
| `intermediate_sg_experiment.*`, `hardwood_sg_trends.*` | S/G trends vs Balakshin 2016 (SI S17) | superseded: those samples are now fitted |
| `ref_4o5_sensitivity.*`, `ref_4o5_hardwood_sensitivity.*`, `sigma_4o5_sensitivity.*` | reference-value sensitivity (SI S12) | superseded |
| `mini_forward.log`, `mini_forward_sweep.csv` | first S/G sweep of the v1.0 generator | superseded |

## Lignin-KMC at its published rates and the cross-generator timing experiments

| files | what | status |
|---|---|---|
| `emergence_sweep.*`, `supply_regime.log`, `supply_regime_sweep.csv` | Lignin-KMC S/G and supply-rate sweeps at the DFT rates (main text §"Behaviour of the published parameter set") | current, but must be reported on one denominator (docs/17 M3) |
| `cross_kmc_matched.*`, `cross_lgs_matched.*`, `cross_ligningraphs_incremental.*`, `cross_wallclock_stability.*`, `cross_algorithm_analysis.log` | timing experiments of the abandoned complexity paper (docs/03) | record; not a linkage comparison |
