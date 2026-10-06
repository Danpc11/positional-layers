# Scripts

Run order, inputs and outputs. Input paths are relative to `$POSLAYERS_DATA` (see `data/MANIFEST.md`); everything else is
written to and read from `$POSLAYERS_RESULTS`. `run_pipeline.sh` runs the scripts in this order. Shared code lives in
`lib_tad.py`, `lib_tumour.py`, `lib_boot.py` (genomic-block bootstrap) and `cohesin_freedman_lane_lib.py` (per-tumour positional score) and is imported; no script reads or executes
another script's source. Figure and table references are to the paper (Fig., Extended Data Fig. = ED, Supplementary
Table = ST).

| Stage | Script | Arguments | Reads | Writes | Used in |
| --- | --- | --- | --- | --- | --- |
| acquire | `make_tss_table.py` | — | Ensembl GRCh37 (pyannotables) | `DATA/hic/genes_hg19_tss.tsv` | Hi-C |
| acquire | `hic_extract.py` | — (network) | Rao 2014 .hic files | `DATA/hic/hic_contacts_*`, `hic_expected_*` | Hi-C |
| acquire | `extract.py` | TCGA cohorts | `tcga/GDC-PANCAN.htseq_counts.tsv.zip` | `expr_<cohort>.npz` | tumours |
| acquire | `cn_continuous.py` | — | `tcga/GDC-PANCAN_cnv.tsv`, `expr_*.npz` | `cn_cont_<cohort>.npz` | tumours |
| atlas | `atlas.py` | GTEx count files | `gtex/`, sample attributes, BioMart GC | `atlas_results.csv` | Fig. 1c,f; Fig. 2b; ED 4a; ST1 |
| atlas | `orient.py` | GTEx count files | as above | `pairs_<tissue>.csv.gz`, `orientation_by_tissue.csv` | Fig. 2c; input to pair analyses |
| atlas | `isochore_law.py` | — | `gtex/` | `isochore_v2.csv` | Fig. 1g; ED 4b; ST8A |
| atlas | `spectral_form.py` | tissues | `gtex/` | `lorentz_v2.csv`, `spec_v2_<tissue>.npy` | Fig. 2e,f; ED 1b; ST8B |
| atlas | `robust_families.py` | — | `pairs_*.csv.gz` | `robust_families.csv` | text (family exclusion) |
| atlas | `gc_autocorrelation.py` | — | `biomart_GRCh38_gene_gc.txt` | `gc_autocorrelation.csv` | ED 1d |
| atlas | `coupling_robustness.py` | GTEx tissues | `gtex/`, sample attributes, BioMart GC | `coupling_robustness.csv`, `coupling_robustness_pairs_<tissue>.csv.gz` | ED 6; ST8E; Supplementary Note 1, section 7 |
| atlas | `spectral_rigour.py` | GTEx tissue | `gtex/`, sample attributes, BioMart GC | `spectral_rigour.csv`, `coupling_by_bp_<tissue>.csv` | ST8F; Supplementary Note 1, sections 1 and 5 |
| eqtl | `eqtl_share.py` | signif-pairs parquet files | `gtex_eqtl/` | `shares_<tissue>.csv.gz` | input |
| eqtl | `eqtl_test.py` | — | `pairs_*`, `shares_*` | `eqtl_cis_test.csv`, `eqtl_tissue_specificity.csv` | Fig. 3c,d,g; ED 1c; ST3 |
| eqtl | `coloc.py` | SuSiE parquet files | `gtex_eqtl/` | `coloc_<tissue>.csv.gz` | input |
| eqtl | `coloc_test.py` | — | `pairs_*`, `shares_*`, `coloc_*` | `coloc_cis_test.csv`, `coloc_dose_response.csv` | Fig. 3e; ST3 |
| eqtl | `predict.py` | GTEx count files | `gtex/`, SuSiE | `pred_pairs_<tissue>.csv.gz` | input |
| eqtl | `eqtl_law.py` | — | `pred_pairs_*`, `gtex_eqtl/` | `eqtl_law_v2_pairs.csv` | Fig. 3f; ST4B |
| eqtl | `eqtl_law_summary.py` | — | `eqtl_law_v2_pairs.csv` | `eqtl_law_v2_summary.csv`, `eqtl_law_strata.csv`, `eqtl_law_block_sensitivity.csv` | Fig. 3f; ST4A,D; ST8D |
| eqtl | `eqtl_scale.py` | — | `eqtl_law_v2_pairs.csv`, `gtex_eqtl/` | `eqtl_scale_pairs.csv`, `eqtl_scale_sensitivity.csv` | ST4C,E |
| architecture | `tad_test.py` | — | `pairs_*`, `shares_*`, `TAD-full/` | `tad_cis_test.csv` | Fig. 4b; ST5 |
| architecture | `tad_bootstrap.py` | — | `pairs_*`, `TAD-full/` | `tad_block_bootstrap.csv` | Fig. 2b; Fig. 4b; ST5 |
| architecture | `pair_attributes.py` | — | `pairs_*`, `shares_*`, `TAD-full/` | `upset_pair_attributes.csv` | Fig. 2d |
| architecture | `hic_test.py`, `hic_test_imr90.py` | — | `hic/`, `gtex/` | `hic_coupling_<cell>.csv.gz` | Fig. 4c,d |
| architecture | `boot_hic.py` | cell, tissue, replicates | `hic_coupling_*`, `gtex/` | `boot_hic_<cell>.csv` | input |
| architecture | `boot_hic_summary.py` | — | `boot_hic_*.csv` | `boot_hic_summary.csv` | Fig. 4d; ST7C |
| architecture | `fric1.py` | — | `hic_coupling_*`, `pairs_*`, `TAD-full/` | `friction1_*.csv` | ED 3; ST7A,B |
| architecture | `contact_exponent_within.py` | — | `hic_coupling_*` | `contact_exponent_within.csv` | Supplementary Note 1, section 6 |
| architecture | `ctcf_barrier.py` | `GM12878` or `IMR90` | `hic_coupling_*`, `hic/` loop lists with motifs | `ctcf_barrier_<cell>.csv`, `ctcf_barrier_means_<cell>.csv` | Fig. 4e; ED 5c,d; ST9B,C |
| architecture | `active_gene_barrier.py` | GTEx tissue | `gtex/`, sample attributes, BioMart GC | `active_gene_barrier_<tissue>.csv`, `active_gene_barrier_means_<tissue>.csv` | Fig. 4f; ED 5a,b; ST9D,E |
| tumours | `cohesin_v2.py` | — | `tcga/`, `expr_*` | `cohesin_v2_results.csv` (GISTIC, sensitivity) | ST6B |
| tumours | `cont_tests.py` | comma-separated cohorts | `tcga/`, `expr_*`, `cn_cont_*` | `aneuploidy_continuousCN.csv`, `cohesin_continuousCN.csv` | ED 2a,b; ST6 |
| tumours | `aneuploidy.py`, `replicate.py`, `tad_tumour.py`, `fric2.py` | — | `tcga/`, `expr_*`, `TAD-full/` | `aneuploidy_scaling.csv`, `replication_LAML_GBM.csv`, `tad_tumour_results.csv`, `friction2_*.csv` | Fig. 2b; ED 2c; ST6 |
| tumours | `tad_clustering_contrast.py` | — | `friction2_per_sample_z.csv` | `tad_clustering_contrast.csv` | ED 2c; ST6G |
| tumours | `cohesin_freedman_lane.py` | cohorts, permutations | as `cont_tests.py` | `cohesin_freedman_lane.csv` | Fig. 4g; ST6A |
| tumours | `beataml_freedman_lane.py` | — | `beataml/` | `beataml_freedman_lane.csv` | Fig. 4g; ST6C |
| tumours | `meta_stag2.py` | — | the two tables above | `stag2_meta.csv` | Fig. 4g; ST6C |
| tumours | `tumour_baseline_sensitivity.py` | cohorts, permutations | as `cohesin_freedman_lane.py` | `tumour_baseline_sensitivity.csv` | ST6H |
| tumours | `dosage_prediction.py` | comma-separated cohorts | `expr_*`, `cn_cont_*` | `dosage_prediction.csv`, `dosage_prediction_splits.csv` | Fig. 3b; ST9A,F |
| tumours | `dosage_controls.py` | comma-separated cohorts | `expr_*`, `cn_cont_*` | `dosage_controls.csv` | ST9G |
| simulations | `theory_sim.py` | — | — | `sim_summary.csv`, `sim_lag_profiles.csv` | Fig. 1b |
| simulations | `decay_recovery.py` | — | — | `sim_decay_recovery.csv` | Fig. 1h |
| simulations | `gc_correlated_sim.py` | — | — | `gc_correlated_cis_sim.csv` | ST8C |
| liver | `pilot_spectra.py`, `pilot_domain_tests.py`, `pilot_liver_gc.py` | — | Liver_Spectra outputs (`LIVER_SPECTRA`) | `p1_spectra.csv`, `p4_domain_scale_tests.csv`, `p5a_liver_gc.csv` | Fig. 1d,e; ED 1a |
| schmitt (optional) | `reduce_fithic.py` | Fit-Hi-C folder, output folder, codes | GSE87112 Fit-Hi-C files | `schmitt_gene_pairs/<code>_gene_pairs.csv.gz` | input |
| schmitt (optional) | `contact_exponent_schmitt.py` | codes | `schmitt/schmitt_gene_pairs/`, `gtex/` | `contact_exponent_schmitt.csv` | Supplementary Note 1, section 6 |
| — | `simulate_demo.py` | `--decay`, `--gc-bias` | — | prints | web simulator |

Figures: `figures/code/fig1.py` to `fig5.py`, `ed1_3.py` (Extended Data Figs 1–3), `ed4.py`, `ed5.py` and `ed6.py`. Schematic panels of Figs 1a, 2a, 3a and 4a are drawings in `figures/schematics/`, placed by the figure scripts;
Fig. 5 is drawn in code by `fig5.py`.

## Analyses kept for peer review (`scripts/review/`, stage `review`)

Perturbation analyses that are not part of the paper: `bystander_liver.py` (liver fibrosis), `predict_first.py`,
`test_edits.py`, `edit_cis.py` and `edit_predictions.py` (gene editing, GSE264491), `drug_test.py` (SLAM-seq, GSE100708),
`crispri_test.py` and `crispri_analyse.py` (genome-wide CRISPRi), and `perturbation_detectability.py`.
