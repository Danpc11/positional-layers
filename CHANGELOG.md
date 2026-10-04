# Changelog

## 1.2.0 (2026-10-04)

Release archived on Zenodo with the submitted manuscript.

- Figures: panel letters at the outer left edge of each panel, aligned by row; legends moved automatically when they
  cover data (`figures/code/style.py`). Fig. 2b is a circular overview of all tissues and tumour cohorts, Fig. 2d an
  UpSet plot of pair attributes, Fig. 3f an own-minus-other-tissue heatmap. Schematics (Figs. 1a, 2a, 4a, 6) now state
  hypotheses as hypotheses and use the revised domain scale (7-24 genes).
- New scripts: `pair_attributes.py` (Fig. 2d), `gc_autocorrelation.py` (Extended Data Fig. 1d) and `decay_recovery.py`
  (Fig. 1h). `theory_sim.py` no longer writes `sim_decay_recovery.csv`: it used the fixed-window estimator replaced in
  1.1.0, so running the pipeline reproduced the superseded panel 1h. No figure reads raw
  data any more: every figure can be rebuilt from the deposited source data.
- Removed dead code in `fig1.py` that read the BioMart export and computed an unused autocorrelation.
- Citation metadata: `CITATION.cff`, `.zenodo.json`; licence names the authors; package, pyproject and metadata
  versions are checked against each other by `tests/test_metadata.py`.
- README and simulator text aligned with the revised manuscript (title; GC-associated rather than technical layer; the
  isochore test is out of sample; the saturation account is a hypothesis).

## 1.1.0

Answer to the external code and methods review; see `REVIEW_RESPONSE.md`. Configurable paths, shared libraries instead
of `exec`, producers for every table, re-analyses (out-of-sample isochore law, log-scale spectral fits, Freedman-Lane
permutation, donor x genomic-block bootstraps, eQTL law on a common scale), and a static check for undefined names in CI.

## 1.0.0

First public version.
