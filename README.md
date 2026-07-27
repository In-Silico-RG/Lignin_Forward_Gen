# Lignin_Forward_Gen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A minimal, transparent **forward** kinetic Monte Carlo generator that predicts the
native lignin inter-unit **linkage inventory** (β-O-4, β-5, β-β, β-1, 5-5, 4-O-5)
from monomer composition (S/G ratio) alone. Unlike *backward* generators, which are
given target bond fractions and reproduce them, this generator evolves structure from
radical-coupling propensities, so the linkage distribution is an emergent *output*.

This repository contains the generator, its calibration and validation scripts, the
cross-generator comparison harness, and the reference data and result logs supporting
the paper below.

## Paper

> M. Y. Combariza and A. F. Combariza, *Predicting linkage-resolved lignin structure
> across feedstocks with a chemically calibrated forward generator* (submitted to
> *ACS Sustainable Chemistry & Engineering*).

If you use this software, please cite the paper and the archived release (see
[Citation](#citation)). `CITATION.cff` carries both in machine-readable form.

## Requirements

- Python 3.9+
- NumPy, SciPy
- [`ligninkmc`](https://github.com/michaelorella/ligninkmc) (the Lignin-KMC forward
  reference model; used in the Section 3.1 comparison and the cross-generator scripts)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

The core generator, calibration, and validation depend only on NumPy and SciPy.
`ligninkmc` (and, for the full cross-comparison, the LigninGraphs and LGS backward
generators) are only needed for the `code/cross_*.py` scripts.

## Repository layout

```
code/
  mini_forward_generator.py   The forward KMC generator (the core contribution)
  calibrate_full.py           Nelder-Mead calibration to the three wood anchors
  validate.py                 Leave-one-out cross-validation
  forward_emergence_sweep.py  Emergent S/G trends without fitting
  param_recovery_experiment.py, seed_averaged_refit.py
                              Identifiability and seed-stability of the fit
  *_experiment.py, *_sensitivity*.py
                              The individual mechanism and sensitivity tests
  cross_*.py                  Matched cross-generator comparison harness
docs/experiments.md           Every run behind the paper, with its findings
data/                         Literature reference tables used by the scripts
results/                      Calibration/validation logs and CSV outputs
requirements.txt
CITATION.cff                  How to cite the software and the paper
.zenodo.json                  Deposition metadata for the archived release
LICENSE                       MIT
```

`docs/experiments.md` is the full record: what each script does, the command that
produces each committed log, and what the run showed. The sections below cover the
main line of the paper only.

## Reproducing the paper

All commands are run from the repository root. Runs are deterministic (fixed seeds),
so your numbers should match the committed logs in `results/`.

### 1. Emergent S/G trends (no fitting)

```bash
python code/mini_forward_generator.py --n 120 --reps 25
```

Sweeps the S/G ratio at the default parameters and writes
`results/mini_forward_sweep.csv`. Shows the emergent β-O-4 gradient and the
softwood C5-condensation ceiling that arise from "syringyl has no free C5" alone.

### 2. Calibration (fit the seven constants)

```bash
python code/calibrate_full.py --objective chi2
```

Fits the coupling constants to the softwood/poplar/high-S anchors under the reported
χ² objective and writes `results/mini_forward_full_calibrated.csv`. Reproduces
`results/calibration_full.log` — fitted `k` values, per-linkage residuals, and the
S/G sweep (mean RMSE ≈ 1.6 points; this fit scores χ² 25.44 / sse 0.00574).

Swap `--objective sse` for the unweighted sensitivity check (reproduces
`results/calibration_full_sse.log`).

### 3. Validation (leave-one-out)

```bash
python code/validate.py --objective chi2
```

Re-fits on two anchors and predicts the held-out third for each (held-out RMSE ≈ 3.4).
The `--objective` **must match** the calibration run.

### 4. Optional — the eighth-parameter (α-O-4) test

```bash
python code/calibrate_full.py --objective chi2 --with-ao4
```

Adds a non-radical α-O-4 chain-joining constant and refits from the seven-parameter
optimum. Reproduces `results/calibration_full_ao4.log`: `k_ao4` collapses to its
floor (0.0149) and the fit is unchanged (χ² 25.44 → 25.41) — the model declines the
eighth parameter.

### 5. Identifiability of the fitted constants

```bash
python code/param_recovery_experiment.py
python code/seed_averaged_refit.py
```

Refits synthetic data from known constants to see which of the seven the inventory
actually pins. Four come back determined; `k_bb` and `k_bo4_nuc` do not, and are
reported as undetermined rather than read as chemistry. The second script averages the
objective over seeds and reports that fit alongside the adopted one.

### 6. Optional — cross-generator comparison

The `code/cross_*.py` scripts run the forward generator against the Lignin-KMC,
LigninGraphs, and LGS backward generators on matched parameterizations. These require
the respective generators installed; see the header of each script.

### Other switches on `calibrate_full.py`

`--with-knuc` frees continuous nucleation, `--with-oligobb` allows β-β between
oligomers, `--oo-scale` sets the oligomer–oligomer coupling scale, and `--ref-4o5`
overrides the softwood 4-O-5 reference. Each one and its result is covered in
`docs/experiments.md`.

## Data

`data/` holds the literature reference tables the scripts read — currently the
Balakshin et al. (2016) hardwood inventories used by `code/hardwood_sg_trends.py`,
with sources cited per row. The three calibration anchors are defined in the code and
their provenance is documented in the paper's Supporting Information.

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
  doi     = {10.5281/zenodo.XXXXXXX}
}
```

Each release is archived on Zenodo and gets its own DOI; `10.5281/zenodo.XXXXXXX`
above is the concept DOI, which always resolves to the newest version. Repository
metadata for the deposition lives in `.zenodo.json`.

## Contact

Aldo F. Combariza — <aldo.combariza@unisucre.edu.co>
