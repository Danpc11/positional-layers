# Scripts

Run order, inputs and outputs. All paths are relative to `$POSLAYERS_DATA` (inputs, see `data/MANIFEST.md`) or
`$POSLAYERS_RESULTS` (everything else). `run_pipeline.sh` runs them in this order. Shared code lives in `lib_tad.py` and
`lib_tumour.py` and is imported; no script reads or executes another script's source.

| Stage | Script | Arguments | Reads | Writes | Used in |
| --- | --- | --- | --- | --- | --- |
| acquire | `make_tss_table.py` | — | Ensembl GRCh37 (pyannotables) | `DATA/hic/genes_hg19_tss.tsv` | Hi-C |
| acquire | `hic_extract.py` | — (network) | Rao 2014 .hic files | `DATA/hic/hic_contacts_*`, `hic_expected_*` | Hi-C |
| acquire | `extract.py` | TCGA cohorts | `tcga/GDC-PANCAN.htseq_counts.tsv.zip` | `expr_<cohort>.npz` | tumours |
| acquire | `cn_continuous.py` | — | `tcga/GDC-PANCAN_cnv.tsv`, `expr_*.npz` | `cn_cont_<cohort>.npz` | tumours |
| atlas | `atlas.py` | GTEx count files | `gtex/`, sample attributes, BioMart GC | `atlas_results.csv` | Fig. 1c,e,f; 2b; ED 4a |
| atlas | `orient.py` | GTEx count files | as above | `pairs_<tissue>.csv.gz`, `orientation_by_tissue.csv` | Fig. 2c; input to most pair analyses |
| atlas | `isochore_law.py` | — | `gtex/` | `isochore_v2.csv` | Fig. 1g; ED 4b; ST9A |
| atlas | `spectral_form.py` | tissues | `gtex/` | `lorentz_v2.csv`, `spec_v2_<tissue>.npy` | Fig. 2e,f; ED 1b; ST9B |
| atlas | `robust_families.py` | — | `pairs_*.csv.gz` | `robust_families.csv` | text (family exclusion) |
| eqtl | `eqtl_share.py` | signif-pairs parquet files | `gtex_eqtl/` | `shares_<tissue>.csv.gz` | |
| eqtl | `eqtl_test.py` | — | `pairs_*`, `shares_*` | `eqtl_cis_test.csv`, `eqtl_tissue_specificity.csv` | Fig. 3b,c,f; ED 1c |
| eqtl | `coloc.py` | SuSiE parquet files | `gtex_eqtl/` | `coloc_<tissue>.csv.gz` | |
| eqtl | `coloc_test.py` | — | `pairs_*`, `shares_*`, `coloc_*` | `coloc_cis_test.csv`, `coloc_dose_response.csv` | Fig. 3d |
| eqtl | `predict.py` | GTEx count files | `gtex/`, SuSiE | `pred_pairs_<tissue>.csv.gz` (pair list for the eQTL law); also `isochore_law.csv`, the superseded metric kept for comparison | |
| eqtl | `eqtl_law.py` | — | `pred_pairs_*`, `gtex_eqtl/` | `eqtl_law_v2_pairs.csv` | Fig. 3e; ST4B |
| eqtl | `eqtl_law_summary.py` | — | `eqtl_law_v2_pairs.csv` | `eqtl_law_v2_summary.csv` | ST4A |
| architecture | `tad_test.py` | — | `pairs_*`, `shares_*`, `TAD-full/` | `tad_cis_test.csv` | Fig. 4b; ST5 |
| architecture | `tad_bootstrap.py` | — | `pairs_*`, `TAD-full/` | `tad_block_bootstrap.csv` | Fig. 2b; 4b; ST5 |
| architecture | `pair_attributes.py` | — | `pairs_*`, `shares_*`, `TAD-full/` | `upset_pair_attributes.csv` | Fig. 2d |
| architecture | `hic_test.py`, `hic_test_imr90.py` | — | `hic/`, `gtex/` | `hic_coupling_<cell>.csv.gz` | Fig. 4c,d |
| architecture | `boot_hic.py` | cell, tissue, replicates | `hic_coupling_*`, `gtex/` | `boot_hic_<cell>.csv` | |
| architecture | `boot_hic_summary.py` | — | `boot_hic_*.csv` | `boot_hic_summary.csv` | Fig. 4d,f; ST8D |
| architecture | `derive.py` | — | `hic_coupling_*` | `derivation_hub_model.csv` | Supplementary Note 1.5; ST8C |
| architecture | `fric1.py` | — | `hic_coupling_*`, `pairs_*`, `TAD-full/` | `friction1_*.csv` | ED 3; ST8A,B |
| tumours | `cohesin_v2.py` | — | `tcga/`, `expr_*` | `cohesin_v2_results.csv` (GISTIC, sensitivity) | ST6B |
| tumours | `cont_tests.py` | comma-separated cohorts | `tcga/`, `expr_*`, `cn_cont_*` | `aneuploidy_continuousCN.csv`, `cohesin_continuousCN.csv` | ED 2a,b |
| tumours | `aneuploidy.py`, `replicate.py`, `tad_tumour.py`, `fric2.py` | — | `tcga/`, `expr_*`, `TAD-full/` | `aneuploidy_scaling.csv`, `replication_LAML_GBM.csv`, `tad_tumour_results.csv`, `friction2_tumour_tad_clustering.csv` | Fig. 4e; ED 2c; ST6 |
| tumours | `cohesin_freedman_lane.py` | cohorts, permutations | as `cont_tests.py` | `cohesin_freedman_lane.csv` | Fig. 4g; ST6A |
| tumours | `beataml_freedman_lane.py` | — | `beataml/` | `beataml_freedman_lane.csv` | Fig. 4g; ST6C |
| tumours | `meta_stag2.py` | — | the two tables above | `stag2_meta.csv` | Fig. 4g; ST6C2 |
| perturbations | `bystander_liver.py` | — | `pairs_liver`, `liver_tables/` | `bystander_liver_*.csv`, `bystander_block_bootstrap.csv` | Fig. 5b |
| perturbations | `predict_first.py` → `test_edits.py` → `edit_cis.py` | — | `pairs_whole_blood`, `edit/` | `predictions_*.csv`, `edit_*.csv` | Fig. 5c,d |
| perturbations | `drug_test.py` | read threshold | `slam/`, `beataml/` | `drug_cis_results_thr<k>.csv` | Fig. 5f; ST7 |
| perturbations | `crispri_test.py` → `crispri_analyse.py` | — | `perturb/`, `beataml/` | `crispri_pairs.csv.gz` | Fig. 5e |
| simulations | `theory_sim.py` | — | — | `sim_*.csv` | Fig. 1b,h |
| simulations | `gc_correlated_sim.py` | — | — | `gc_correlated_cis_sim.csv` | ST9C |
| liver | `pilot_spectra.py`, `pilot_domain_tests.py`, `pilot_liver_gc.py` | — | Liver_Spectra outputs (`LIVER_SPECTRA`) | `p1_spectra.csv`, `p4_domain_scale_tests.csv`, `p5a_liver_gc.csv` | Fig. 1c-e; ED 1a |
| — | `simulate_demo.py` | `--decay`, `--gc-bias` | — | prints | web simulator |

Removed in 1.1.0: `predict_cal.py` (replaced by `eqtl_law.py`), `derive2.py` (printed results that were then typed into a
table; replaced by `boot_hic.py`), `cohesin_test.py` (superseded), `lorentz.py` (replaced by `spectral_form.py`, which fits
and evaluates on the same scale) and `beat_test.py` (replaced by `beataml_freedman_lane.py`). See `REVIEW_RESPONSE.md`.
