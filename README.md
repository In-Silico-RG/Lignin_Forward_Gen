# Lignin_Forward_Gen

A minimal, transparent **forward** kinetic-Monte-Carlo generator of lignin
structures, with calibration and cross-validation against native lignin.

Lignin's value for valorization is governed by its inter-unit linkage
distribution — above all the β-O-4 aryl-ether content, which sets the monomer
yield attainable by depolymerization, and the condensed C–C linkages (5-5, 4-O-5)
that resist it. This tool **predicts** the full linkage inventory (β-O-4, β-5,
β-β, β-1, 5-5, 4-O-5) from monomer composition (S/G ratio) alone, rather than
imposing a target distribution.

## Forward vs. backward

- **Backward** generators (LGS, LigninGraphs) are *given* the target bond
  fractions and sample/reject structures until they match — they reproduce a
  known lignin but cannot predict an unknown one.
- **Forward** generators (this tool; Lignin-KMC) evolve a structure from local
  radical-coupling propensities, so the linkage distribution *emerges* and is
  genuinely predicted.

Every coupling propensity here is an explicit, tunable parameter, so each
deviation from experiment is a nameable piece of chemistry.

## What's here

```
code/
  mini_forward_generator.py   the forward KMC generator (numpy only)
  calibrate.py                fit the beta-side-chain distribution
  calibrate_full.py           fit the full linkage inventory (fraction of total)
  validate.py                 native-only leave-one-out cross-validation
results/                      sweep outputs, calibration + validation logs/CSVs
  fitted_params.json          calibrated coupling parameters
docs/experiments.md           full experimental record (emergence, calibration, validation)
```

## Model in brief

Units carry explicit reactive positions (phenolic O4, C5 [G/H only], C1, and a
one-shot β radical donor). Bonds form by radical cross-coupling: β-O-4 (chain
extension), β-5, β-1, β-β (monomer–monomer), and the condensed 5-5 / 4-O-5
(oligomer–oligomer). Two mechanistic choices are load-bearing:

1. **End-restricted β-coupling** — β-5/β-1 couple only to phenolic end units, so
   their frequency is set by rate constants, not by the large internal-C5 pool.
2. **Decoupled nucleation** — a short nucleation phase forms dimers (the origin of
   β-β), then endwise growth with a standing pool of one monomer keeps chains
   linear (a single large monomer pool forces condensation, because big oligomers
   can only merge via C–C bonds).

The syringyl-blocks-C5 rule (S has no free C5) is built in, so the softwood
condensation ceiling and the S/G → β-O-4 gradient emerge rather than being imposed.

## Results (calibrated)

Fit to three native wood references (softwood, poplar, high-S hardwood); full
linkage inventory reproduced within ~3 percentage points, near-linear topology,
and **mean held-out RMSE 3.0 points** under leave-one-out cross-validation.
Industrial pulping lignins are deliberately excluded from validation (pulping
alters native structure). See `docs/experiments.md`.

## Usage

Requires Python 3 with `numpy` (generator) and `scipy` (calibration/validation).

```bash
python code/mini_forward_generator.py     # S/G sweep vs NMR anchors
python code/calibrate_full.py             # full-inventory calibration
python code/validate.py                   # leave-one-out cross-validation
```

## Related work

Builds on and compares against **Lignin-KMC** (Orella et al., *ACS Sustainable
Chem. Eng.* 2019) with barriers from Gani et al. (*ACS Sustainable Chem. Eng.*
2019), and **LigninGraphs** (Wang et al. 2022). β-O-4 → monomer-yield relationship
from Bourmaud, Luterbacher et al. (*Angew. Chem. Int. Ed.* 2024).

## License

MIT — see [LICENSE](LICENSE).
