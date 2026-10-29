# Lignin_Forward_Gen

A minimal, transparent **forward** kinetic Monte Carlo generator that predicts the
native lignin inter-unit **linkage inventory** (β-O-4, β-5, β-β, β-1, 5-5, 4-O-5)
from monomer composition (S/G ratio and H fraction) alone. Unlike *backward*
generators, which are given target bond fractions and reproduce them, this generator
evolves structure from radical-coupling propensities, so the linkage distribution is
an emergent *output*.

This repository contains the generator, its calibration and validation scripts, the
cross-generator comparison harness, and the reference data and result logs
supporting the paper below. Its continuing mission: to explore strange new
feedstocks, to seek out new linkage distributions and new condensations, to
boldly calibrate where no model has calibrated before.

## Paper

> M. Y. Combariza and A. F. Combariza, *Predicting linkage-resolved lignin structure
> across feedstocks with a chemically calibrated forward generator* (submitted to
> *ACS Sustainable Chemistry & Engineering*).

If you use this software, please cite the paper and the archived release (see
[Citation](#citation)). `CITATION.cff` carries both in machine-readable form.

## Requirements

Two separate environments, because the two generators being compared don't share a
Python version:

- **`lignin`** — Python 3.9+, NumPy, SciPy. Runs the forward generator itself and
  every calibration, ablation, validation and identifiability script.
  ```bash
  conda create -n lignin python=3.9 numpy scipy
  ```
- **`lignin_kmc`** — Python 3.6, plus [`ligninkmc`](https://github.com/michaelorella/ligninkmc)
  (the Lignin-KMC reference model, Orella et al. 2019). Only needed for the
  `calibrate_kmc_ref.py` comparison and the `cross_*.py` scripts; not needed to run
  the forward generator or reproduce its own calibration.
  ```bash
  conda create -n lignin_kmc python=3.6
  conda run -n lignin_kmc pip install ligninkmc
  ```

The committed logs in `results/` were produced with NumPy 2.3.3 / SciPy 1.16.2 under
`lignin`. The generator is seeded and deterministic; the Nelder–Mead calibration can
land on a marginally different optimum under a different SciPy version.

## Repository layout

```
code/
  mini_forward_generator.py   The forward generator itself (814 lines, 18 functions) — the core contribution
  reference_set.py            The 60-sample reference set and the ratio-to-beta-O-4 objective
  calibrate_ref.py            Calibrate the seven constants (the current, v1.1 fit)
  calibrate_kmc_ref.py        Calibrate Lignin-KMC on identical footing (needs the lignin_kmc env)
  identifiability_ref.py      Which of the seven constants the data determine
  score_validation.py         Score the eleven lignins never used in the fit
  weld_decomposition.py       Split condensed bonds into structural welds vs. elective attachments
  connectivity_experiment.py  Dependence of the fit on simulation size (DP)
  structure_from_run.py       Real 2D/3D structures from a recorded run (used by the tutorial)
  *_experiment.py, *_sensitivity*.py
                               Individual ablations and sensitivity tests
  cross_*.py                  Matched cross-generator comparison harness (LigninGraphs, LGS)
  calibrate_full.py           The original three-anchor calibration (superseded by calibrate_ref.py;
                               kept for provenance, not what the current paper reports)
docs/experiments.md           What each script does and what its committed log showed
data/reference_inventories.csv
                               The 60 fitted + 11 validation-only lignins, with provenance per row
results/                      Calibration/validation logs and the adopted fitted constants
requirements.txt
CITATION.cff                  How to cite the software and the paper
LICENSE                       MIT
```

`code/` is 6,432 lines across 35 files. This repository is the archived snapshot of the
generator itself; the full dated running log, the manuscript and its SI, and the
interactive tutorial live in the sibling repository
[`Testing_Lignin_Structure_Generators`](https://github.com/In-Silico-RG/Testing_Lignin_Structure_Generators),
where the paper's own analysis is developed. The sections below cover the main line of
the paper only.

## Reproducing the paper

All commands run from the repository root. Runs are seeded and deterministic, so
your numbers should match the committed logs in `results/` — engage.

### 1. Calibration — fit the seven constants on the 60-sample reference set

```bash
conda run -n lignin python code/calibrate_ref.py --source-weight --reps 8 --tag _v11e_sw
```

Fits `k_β5, k_β1, k_ββ, k_55, k_4o5, n_nuclei, k_β-O-4,nuc` against the 60 native
lignins (`data/reference_inventories.csv`), source-weighted so no single laboratory
sets the fit, with `k_β-O-4 = 8` pinned. Writes `results/calibration_ref_v11e_sw.log`
and `results/fitted_params_ref_v11e_sw.json` — the exact names already committed
(χ² = 244.1), so this reruns the paper's adopted fit under its own name rather than
a new one. Leaving `--tag` off still works, but lands on `fitted_params_ref_sw.json`
instead (the plain default tag is just `_sw`, not `_v11e_sw` — that part names which
dated round of reruns it was, not something the script picks on its own).

### 2. Validation — leave-one-sample-out and leave-one-source-out

```bash
# one fold, repeat per sample_id (or script it over data/reference_inventories.csv)
conda run -n lignin python code/calibrate_ref.py --source-weight --hold-out anchor_softwood --tag _loo-anchor_softwood

# collect all folds against the full-fit log
python code/summarize_loo.py 'results/calibration_ref_loo-*.log' results/calibration_ref_v11e_sw.log
```

Reproduces the held-out/in-sample ratios in the paper (LOO 0.98, LOSO 0.99).
Leave-one-source-out holds out every sample from one laboratory at once instead of
one sample.

### 3. Ablations — which rule does the work

```bash
conda run -n lignin python code/calibrate_ref.py --source-weight --fix end_restricted=0 --tag _ablate_endres
conda run -n lignin python code/calibrate_ref.py --source-weight --fix nucleation_on=0  --tag _ablate_single
```

Removing the phenolic-end restriction costs 23% χ²; removing the two-stage
nucleation structure is a near-tie in the source-weighted total but not sample by
sample (the sibling `Testing_Lignin_Structure_Generators` repo's docs/15 has the full
breakdown — the discussion the data actually leave open).

### 4. Eighth constants — declined on purpose

```bash
conda run -n lignin python code/calibrate_ref.py --source-weight --with-ao4     --start results/fitted_params_ref_v11e_sw.json --tag _ao4
conda run -n lignin python code/calibrate_ref.py --source-weight --with-oligobb --start results/fitted_params_ref_v11e_sw.json --tag _oligobb
```

(explicit `--tag` on both — without one, the default is just `_sw` for any
source-weighted run, so these two would otherwise overwrite the same output file.)

Both candidate reactions (a chain-joining α-O-4 ether, a re-oxidised oligomer β-β)
are tested and declined — the fit does not need them, and the paper says why.

### 5. Identifiability — which of the seven the data actually pin down

```bash
conda run -n lignin python code/identifiability_ref.py --params results/fitted_params_ref_v11e_sw.json
```

Sweeps each constant over its full range with the other six fixed, compared against
the model's own run-to-run noise. Six of seven are determined; the seventh
(`k_β-O-4,nuc`) barely moves the predictions.

### 6. Validation-only samples — eleven lignins never in the fit

```bash
conda run -n lignin python code/score_validation.py --params results/fitted_params_ref_v11e_sw.json
```

Nothing is fitted here — this just scores the held-out samples (cork, grasses,
*Miscanthus* at three processing levels) against the adopted constants.

### 7. Welds and size-dependence

```bash
conda run -n lignin python code/weld_decomposition.py      --params results/fitted_params_ref_v11e_sw.json
conda run -n lignin python code/connectivity_experiment.py --params results/fitted_params_ref_v11e_sw.json
```

Splits condensed bonds into structurally forced welds vs. elective C5 attachments,
and checks how much the fit depends on the simulated system size (120 units) —
reported honestly as part of the calibrated model, not hidden.

### 8. Lignin-KMC on identical footing

```bash
conda run -n lignin_kmc python code/calibrate_kmc_ref.py --source-weight --procs 24
```

Calibrates Lignin-KMC's own multipliers on the same 60 samples, same objective,
same number of free parameters (seven). Needs the separate `lignin_kmc`
environment — see Requirements.

### 9. Optional — cross-generator comparison

The `code/cross_*.py` scripts run the forward generator against LigninGraphs and
LGS (backward generators) on matched parameterizations. These require the
respective generators installed; see the header of each script.

## Data

`data/reference_inventories.csv` is the full reference set: 60 native lignins used
in the fit and 11 kept for validation only, one row per sample, with species,
method, S/G, H fraction, every linkage it reports, and the resinol (β-β)
units-vs-bonds convention applied to that source. Provenance for every row is in the
sibling repo's docs/15 §5 and the paper's SI §S4.

## License

Released under the [MIT License](LICENSE).

## Citation

```bibtex
@software{lignin_forward_gen,
  author  = {Combariza, Marianny Y. and Combariza, Aldo F.},
  title   = {Lignin\_Forward\_Gen: a forward kinetic Monte Carlo generator for
             linkage-resolved lignin structure},
  year    = {2026},
  url     = {https://github.com/In-Silico-RG/Lignin_Forward_Gen},
  doi     = {10.5281/zenodo.21630116}
}
```

Each release is archived on Zenodo and gets its own DOI; `10.5281/zenodo.21630116`
above is the concept DOI, which always resolves to the newest version.

## Contact

Aldo F. Combariza — <aldo.combariza@unisucre.edu.co>

Live long, and reproduce.
