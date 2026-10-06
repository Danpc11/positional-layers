# Changelog

## Unreleased

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
