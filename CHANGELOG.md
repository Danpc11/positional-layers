# Changelog

## 1.3.3

- Re-running a unit (tissue or cohort) now replaces ALL its earlier rows, also when the new output has fewer rows or none:
  `poslayers.config.replace_groups_csv`, used by `eqtl_law.py`, `eqtl_scale.py`, `orient.py`, `cont_tests.py`,
  `fric2.py` and `dosage_prediction.py`. Before, `upsert_csv` with pair-level keys kept rows that a new filter had
  removed, and `cont_tests.py` skipped the write when a cohort had no rows. `upsert_csv` is kept for incremental writes.
- `boot_hic.py` starts a fresh replicate file unless `POSLAYERS_RESUME=1`, so a run with fewer replicates does not keep
  earlier ones.
- Tests for a filter that removes a pair, a tissue left without pairs, and the documented behaviour of `upsert_csv`.
- Removed a dead first definition of `ar1` in `theory_sim.py` and duplicate imports.

## 1.3.2

- `dosage_prediction.py`: the dosage model (intercept, linear and quadratic terms) is fitted in the training half only
  and applied without refitting to the validation half; the prediction is the covariance of the fitted dosage components
  (validation copy number, training coefficients), so nothing is fitted in the validation half. Ratios change from
  0.99-1.07 to 0.85-1.07 (split shown) and 0.71-1.07 across ten splits. A linear-only prediction is kept as a
  sensitivity analysis. Rows are replaced by cohort, split and lag, so running one cohort keeps the other.
- Fig. 3f: stratum means now have intervals from the tissue-stratified block bootstrap (`eqtl_law_strata.csv`) instead of
  pair-level standard errors.
- `eqtl_law_summary.py`: sensitivity to the block size (5, 10, 20 Mb) and to blocks coordinated across tissues
  (`eqtl_law_block_sensitivity.csv`; Supplementary Table 8D).
- Re-running a script now recomputes everything; `POSLAYERS_RESUME=1` reuses tissues or replicates already written
  (atlas, orient, isochore_law, eqtl_law, eqtl_scale, boot_hic).
- `lib_boot.py`: genes without a GRCh38 position get blocks of their own (they were placed in block 0); block size can be
  set with `POSLAYERS_BLOCK_MB`; `p_block` documented as a normal approximation.
- `run_pipeline.sh acquire` checks the contact and expected files of both cell lines before skipping `hic_extract.py`.
- `*.zip` added to `.gitignore`; release archives are distributed outside the repository.

## 1.3.1

- Fig. 4e,f: intervals of the ratios to the reference group are now computed within each bootstrap replicate, so they
  include the uncertainty of the reference (`ctcf_barrier.py`, `active_gene_barrier.py` write `relative`,
  `relative_ci_low` and `relative_ci_high`).
- `dosage_prediction.py` repeats the out-of-sample prediction for ten random splits of tumours
  (`dosage_prediction_splits.csv`; Supplementary Table 9F); the split shown in Fig. 3b is unchanged.
- Fig. 5 schematic descriptions updated to the revised text.
- Fig. 2b is a ranked dot plot (tissues on one common scale, tumour cis excess on its own scale) instead of the circular
  overview; Fig. 2d names each attribute combination; the Fig. 2f legend sits above the plot.

## 1.3.0

Code for the restructured manuscript (5 figures, 5 Extended Data figures, Supplementary Tables 1-9).

- Shared-source model: `poslayers.source_covariance`, and `dosage_prediction.py`, which predicts the copy-number layer of
  tumours in held-out samples (Fig. 3b).
- Architectural barriers: `ctcf_barrier.py` (coupling across convergent CTCF loop anchors; Fig. 4e) and
  `active_gene_barrier.py` (coupling across an active intervening gene; Fig. 4f), each writing the adjusted means for
  the figure itself; Extended Data Fig. 5.
- Contact exponents for Supplementary Note 1, section 6: `contact_exponent_within.py` (25-kb maps) and the optional
  `schmitt` stage (`reduce_fithic.py`, `contact_exponent_schmitt.py`; 40-kb tissue maps).
- Figures: Fig. 5 is the model; schematic panels are reserved for drawings (Figs 1a, 2a, 3a, 4a, 5). Extended Data
  Fig. 4 no longer has an intervention panel; Extended Data Fig. 5 is new.
- Perturbation analyses moved to `scripts/review/` and the optional `review` stage; they are not in the paper.
  `derive.py` (fixed-hub saturation model) removed.
- Supplementary Tables renumbered: perturbations removed; former 8 and 9 are now 7 and 8; new 9 for the shared-source
  model and the barriers.

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
- Analyses of the scale of the eQTL relation (`eqtl_scale.py`), detectable effects and equivalence for the perturbation
  analyses (`perturbation_detectability.py`), gene-editing predictions under an explicit rule with false-discovery
  control (`edit_predictions.py`), and a direct test of the change in tumour-normal TAD clustering
  (`tad_clustering_contrast.py`); per-quintile intervals for the liver concordance; fixed random seed per tumour cohort.
- Figures can be rebuilt from the source data tables without any raw input; the tables are to be deposited on Zenodo.
- Interactive simulator of the decomposition (`docs/`).
- Citation metadata (`CITATION.cff`, `.zenodo.json`); tests of the identity, the estimators, the edge cases of the API and
  the agreement of release versions; CI fails on any undefined name.

## 1.0.0

First public version.
