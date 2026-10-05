# Changelog

## 1.2.0 (2026-10-04)

Release accompanying the submitted manuscript.

- Analysis scripts for every result, table and figure; shared code in `scripts/lib_tumour.py` and `scripts/lib_tad.py`;
  configurable input, results and figure folders (`poslayers.config`); `run_pipeline.sh` runs everything in order.
- Intervals that account for dependence between gene pairs where implemented: genomic-block bootstraps for eQTL
  sharing, colocalisation, same-TAD effects, the eQTL relation (within tissues), drug slopes and the liver concordance;
  donor x genomic-block bootstraps for the Hi-C results; Freedman-Lane permutation for tumour mutation effects;
  perturbation-clustered errors for CRISPRi. Benjamini-Hochberg q values for the 12 drug comparisons.
- Re-runs replace rows of results tables by key instead of appending (`poslayers.config.upsert_csv`).
- `crispri_test.py` builds the pair table only; the models are fitted in `crispri_analyse.py`, which also reports the
  coupling distribution of the random trans control.
- Figures can be rebuilt from the source data tables without any raw input; the tables are to be deposited on Zenodo.
- Interactive simulator of the decomposition (`docs/`).
- Citation metadata (`CITATION.cff`, `.zenodo.json`); tests of the identity, the estimators, the edge cases of the API and
  the agreement of release versions; CI fails on any undefined name.

## 1.0.0

First public version.
