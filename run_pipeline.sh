#!/usr/bin/env bash
# Reproduces every result, table and figure of the paper in dependency order.
# Inputs: $POSLAYERS_DATA (see data/MANIFEST.md). Outputs: $POSLAYERS_RESULTS (tables) and $POSLAYERS_FIGS (figures).
# Stages can be run alone, e.g. bash run_pipeline.sh tumours. 'all' runs every stage except 'schmitt' and 'review'.
# Full run time is several hours; boot_hic.py and isochore_law.py dominate.
set -euo pipefail
D=${POSLAYERS_DATA:-data}; S=scripts
ONLY=${1:-all}
stage() { [ "$ONLY" = "$1" ] || { [ "$ONLY" = all ] && [ "$1" != schmitt ] && [ "$1" != review ]; }; }
run() { echo "+ python $*"; python "$@"; }

if stage acquire; then                                   # network access; run once
  run $S/make_tss_table.py
  ok=1; for c in GM12878 IMR90; do for k in contacts expected; do [ -f "$D/hic/hic_${k}_${c}.csv.gz" ] || ok=0; done; done
  [ $ok = 1 ] || run $S/hic_extract.py                   # extracts both cell lines
  run $S/extract.py BLCA UCEC STAD COAD GBM LAML
  run $S/cn_continuous.py
fi
if stage atlas; then                                     # Fig. 1c,f,g; Fig. 2; Extended Data Figs 1, 4 and 6; Supplementary Tables 1, 2, 8
  run $S/atlas.py   "$D"/gtex/*_gct.gz
  run $S/orient.py  "$D"/gtex/*_gct.gz
  run $S/isochore_law.py
  run $S/spectral_form.py thyroid nerve_tibial skin_sun_exposed_lower_leg cells_cultured_fibroblasts cells_ebv-transformed_lymphocytes muscle_skeletal
  run $S/robust_families.py
  run $S/gc_autocorrelation.py
  run $S/coupling_robustness.py thyroid cells_ebv-transformed_lymphocytes muscle_skeletal lung
  for t in thyroid nerve_tibial skin_sun_exposed_lower_leg cells_cultured_fibroblasts cells_ebv-transformed_lymphocytes muscle_skeletal; do run $S/spectral_rigour.py $t; done
fi
if stage eqtl; then                                      # Fig. 3c-g; Supplementary Tables 3, 4
  run $S/eqtl_share.py "$D"/gtex_eqtl/*_v11_eQTLs_signif_pairs.parquet
  run $S/eqtl_test.py
  run $S/coloc.py      "$D"/gtex_eqtl/*_v11_eQTLs_SuSiE_summary.parquet
  run $S/coloc_test.py
  run $S/predict.py    "$D"/gtex/*_gct.gz
  run $S/eqtl_law.py
  run $S/eqtl_law_summary.py
  run $S/eqtl_scale.py
fi
if stage architecture; then                              # Fig. 4b-f; Extended Data Figs 3 and 5; Supplementary Tables 5, 7, 9B-E
  run $S/tad_test.py
  run $S/tad_bootstrap.py
  run $S/pair_attributes.py
  run $S/hic_test.py
  run $S/hic_test_imr90.py
  run $S/boot_hic.py GM12878 cells_ebv-transformed_lymphocytes 200
  run $S/boot_hic.py IMR90 cells_cultured_fibroblasts 200
  run $S/boot_hic_summary.py
  run $S/fric1.py
  run $S/contact_exponent_within.py
  for c in GM12878 IMR90; do run $S/ctcf_barrier.py $c; done
  for t in cells_ebv-transformed_lymphocytes thyroid cells_cultured_fibroblasts; do run $S/active_gene_barrier.py $t; done
fi
if stage tumours; then                                   # Fig. 3b; Fig. 4g; Extended Data Fig. 2; Supplementary Tables 6, 9A
  run $S/cohesin_v2.py
  run $S/cont_tests.py BLCA,UCEC,STAD,COAD,GBM
  run $S/aneuploidy.py
  run $S/replicate.py
  run $S/tad_tumour.py
  run $S/fric2.py
  run $S/tad_clustering_contrast.py
  run $S/cohesin_freedman_lane.py BLCA,UCEC 2000
  run $S/beataml_freedman_lane.py
  run $S/meta_stag2.py
  run $S/tumour_baseline_sensitivity.py BLCA,UCEC 2000
  run $S/dosage_prediction.py BLCA,UCEC
  run $S/dosage_controls.py BLCA,UCEC
fi
if stage simulations; then                               # Fig. 1b,h; Supplementary Table 8C
  run $S/theory_sim.py
  run $S/decay_recovery.py
  run $S/gc_correlated_sim.py
fi
if stage liver; then                                     # Fig. 1d,e; Extended Data Fig. 1a; needs LIVER_SPECTRA (data/MANIFEST.md)
  run $S/pilot_spectra.py; run $S/pilot_domain_tests.py; run $S/pilot_liver_gc.py
fi
if stage figures; then
  for f in fig1 fig2 fig3 fig4 ed_validation ed1_3 ed4 ed5 ed6; do run figures/code/$f.py; done
fi
if stage schmitt; then                                   # optional: Supplementary Note 1, section 6 (40-kb tissue Hi-C)
  run $S/reduce_fithic.py "$D"/schmitt/FitHiC_primary_cohort "$D"/schmitt/schmitt_gene_pairs GM12878,imr90
  run $S/contact_exponent_schmitt.py GM12878,imr90
fi
if stage review; then                                    # optional: perturbation analyses kept for peer review, not in the paper
  run $S/review/bystander_liver.py; run $S/review/predict_first.py; run $S/review/test_edits.py; run $S/review/edit_cis.py
  run $S/review/edit_predictions.py; run $S/review/drug_test.py 3; run $S/review/crispri_test.py; run $S/review/crispri_analyse.py
  run $S/review/perturbation_detectability.py
fi
