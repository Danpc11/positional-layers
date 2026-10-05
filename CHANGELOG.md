# Changelog

## 1.2.0 (2026-10-04)

First archived release, accompanying the submitted manuscript.

- Analysis scripts for every result, table and figure; shared code in `scripts/lib_tumour.py` and `scripts/lib_tad.py`;
  configurable input, results and figure folders (`poslayers.config`); `run_pipeline.sh` runs everything in order.
- Inference that accounts for dependence between gene pairs: donor x genomic-block bootstraps, genomic-block bootstraps
  and Freedman-Lane permutation.
- Figures can be rebuilt from the archived source data without any raw input.
- Interactive simulator of the decomposition (`docs/`).
- Citation metadata (`CITATION.cff`, `.zenodo.json`); tests of the identity, the estimators, the edge cases of the API and
  the agreement of release versions; CI fails on any undefined name.

## 1.0.0

First public version.
