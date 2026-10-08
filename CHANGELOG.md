# Changelog

## Unreleased

- Figures reorganised for readability (manuscript v3): Fig. 1 opens with the schematic of the sources removed before
  co-expression is measured and keeps the panels that carry the argument; its validation panels (identity, consensus
  versus tissue-average spectrum, decay-length recovery) form the new Extended Data Fig. 1 (`ed_validation.py`).
  Fig. 2e,f show coupling against physical distance and the two scales in kb/Mb; the Lorentzian spectrum moved to
  Extended Data Fig. 2e. Fig. 3a and Fig. 5 use the new schematics; Fig. 5 is a single-column synthesis
  (`style.save` accepts a width). Extended Data figures are renumbered by first citation (1-7); the script that draws
  each one is unchanged. Schematic images go in `figures/schematics/` (fig1a, fig2a, fig3a, fig4a, fig5).
- `coupling_robustness.py` and `spectral_rigour.py` apply the technical (batch) covariates only when samples outnumber them
  by at least ten, as `atlas.py` does; previously small GTEx tissues could lose all residual variation. The tissues
  analysed in the paper all satisfy the rule, so the reported results are unchanged.
- `dosage_controls.py`: the copy-number prediction with copy number permuted across tumours, with segments shifted along
  chromosomes, and with purity proxies (immune and stromal scores) in the dosage model (Supplementary Table 9G).
- `eqtl_scale.py` also removes 60 expression components, the largest number of hidden factors in the GTEx eQTL model
  (Supplementary Table 4E).
- `spectral_rigour.py`: peaks of the mean profile tested against nulls that preserve local autocorrelation (block
  permutation, AR(1) red noise; null maximum over frequencies), Whittle-likelihood fits of the covariance spectrum, and
  decay scales of coupling against physical distance (Supplementary Note 1, sections 1 and 5; Supplementary Table 8F).
- `tumour_baseline_sensitivity.py`: the primary tumour analysis repeated with a cis-excess baseline beyond the domain
  scale (Supplementary Table 6H); the lag score is shared through `cohesin_freedman_lane_lib.py`.
- Documentation aligned with the article: the "isochore law" is the GC-associated layer and the "eQTL law" is the
  predicted genetic covariance, in the README, the simulator and the package. `gc_layer_correlation` and
  `predicted_genetic_covariance` are the new names; `isochore_law` and `eqtl_law` remain as aliases.
  `saturation_exponent` and the hub model of `fit_contact_law` are marked as not used in the article.
- Duplicated entries removed from `poslayers.__all__`; source-data table count (46) and Extended Data figure count (6)
  updated in the README, the data manifest and `.zenodo.json`.
- `coupling_robustness.py` and Extended Data Fig. 6 (`ed6.py`): coupling of adjacent pairs is compared between Pearson and
  Spearman, between random halves of donors, and with distant pairs matched for expression and GC content, without and
  with removal of 15 expression principal components (Supplementary Note 1, section 7; Supplementary Table 8E).
- Figures are saved at exactly 180 mm width (Nature double column), scaling width and height together, so the font sizes
  in the code are the printed sizes; all text is at least 5.5 pt and raster elements are written at 600 ppi.
- Schematic panels (Figs 1a, 2a, 3a, 4a and 5) are drawn from `figures/schematics/*.png`; the reserved frame is used only
  when a drawing is missing.
- Layout: more height for Fig. 2b and Fig. 3g so that tissue names do not overlap; Fig. 2c legend moved; Fig. 1b label
  states that the landscape share refers to the simulation.
- Extended Data Figs 1a and 2a: legends placed clear of the axes and data (legends can be fixed with `keep_position`);
  Extended Data Fig. 1b: labels of nearby points separated.
- Fig. 5 is drawn in code as vector schematics (numbered layers with their effect on expression, the sign of a shared
  variant's contribution, and coupling arcs for domains, boundaries and active intervening genes).

## 0.1.0 (2026-10-05)

First public release, accompanying the submitted manuscript "Common cis-regulatory inputs shape local gene co-expression
across human tissues" (5 figures, 5 Extended Data figures, Supplementary Tables 1-9).

- `poslayers` package: the spectral identity (periodogram = mean landscape + covariance spectrum), GC correction and the
  out-of-sample GC relation, the eQTL relation, the shared-source covariance and genome simulation, with tests.
- Analysis scripts for every result, table and figure of the paper, run in dependency order by `run_pipeline.sh`
  (stages: acquire, atlas, eqtl, architecture, tumours, simulations, liver, figures); `scripts/SCRIPTS.md` maps each
  script to its inputs, outputs and figure panels.
- Shared-source model tested on copy number: a dosage model fitted in half of the tumours predicts the copy-number
  covariance of the other half without refitting (Fig. 3b).
- Architectural barriers: coupling across convergent CTCF loop anchors and across active intervening genes (Fig. 4e,f).
- Inference that accounts for dependence between gene pairs: genomic-block bootstraps (within tissues and coordinated
  across tissues), joint donor and block bootstraps for Hi-C, Freedman-Lane permutation for tumour mutation effects;
  block-size sensitivity.
- Reproducible re-runs: recomputed units replace all their earlier rows; `POSLAYERS_RESUME=1` reuses finished units.
- Figures rebuild from the source-data tables without raw inputs; schematic panels are reserved for drawings.
- Optional stages: `schmitt` (40-kb tissue Hi-C, Supplementary Note 1) and `review` (perturbation analyses kept for peer
  review, not in the paper).
- Interactive simulator of the decomposition in `docs/`.
