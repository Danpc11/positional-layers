#!/usr/bin/env bash
# Runs every analysis in dependency order, then every figure. Inputs: $POSLAYERS_DATA (see data/MANIFEST.md).
# Outputs: $POSLAYERS_RESULTS (tables) and $POSLAYERS_FIGS (figures). Stages can be run alone: bash run_pipeline.sh tumours
# Full run time is several hours; boot_hic.py and isochore_law.py dominate.
set -euo pipefail
D=${POSLAYERS_DATA:-data}; S=scripts
stage() { [ $# -eq 0 ] || [ "${ONLY:-all}" = all ] || [ "${ONLY}" = "$1" ]; }
ONLY=${1:-all}
run() { echo "+ python $*"; python "$@"; }

if stage acquire; then                                   # network access; run once
  run $S/make_tss_table.py
  [ -f "$D/hic/hic_contacts_GM12878.csv.gz" ] || run $S/hic_extract.py
  run $S/extract.py BLCA UCEC STAD COAD GBM LAML
  run $S/cn_continuous.py
fi
if stage atlas; then
  run $S/atlas.py   "$D"/gtex/*_gct.gz
  run $S/orient.py  "$D"/gtex/*_gct.gz
  run $S/isochore_law.py
  run $S/spectral_form.py thyroid nerve_tibial skin_sun_exposed_lower_leg cells_cultured_fibroblasts cells_ebv-transformed_lymphocytes muscle_skeletal
  run $S/robust_families.py
  run $S/gc_autocorrelation.py
fi
if stage eqtl; then
  run $S/eqtl_share.py "$D"/gtex_eqtl/*_v11_eQTLs_signif_pairs.parquet
  run $S/eqtl_test.py
  run $S/coloc.py      "$D"/gtex_eqtl/*_v11_eQTLs_SuSiE_summary.parquet
  run $S/coloc_test.py
  run $S/predict.py    "$D"/gtex/*_gct.gz
  run $S/eqtl_law.py
  run $S/eqtl_law_summary.py
fi
if stage architecture; then
  run $S/tad_test.py
  run $S/tad_bootstrap.py
  run $S/pair_attributes.py
  run $S/hic_test.py
  run $S/hic_test_imr90.py
  run $S/boot_hic.py GM12878 cells_ebv-transformed_lymphocytes 200
  run $S/boot_hic.py IMR90 cells_cultured_fibroblasts 200
  run $S/boot_hic_summary.py
  run $S/derive.py
  run $S/fric1.py
fi
if stage tumours; then
  run $S/cohesin_v2.py
  run $S/cont_tests.py BLCA,UCEC,STAD,COAD,GBM
  run $S/aneuploidy.py
  run $S/replicate.py
  run $S/tad_tumour.py
  run $S/fric2.py
  run $S/cohesin_freedman_lane.py BLCA,UCEC 2000
  run $S/beataml_freedman_lane.py
  run $S/meta_stag2.py
fi
if stage perturbations; then
  run $S/bystander_liver.py
  run $S/predict_first.py
  run $S/test_edits.py
  run $S/edit_cis.py
  run $S/drug_test.py 3
  run $S/crispri_test.py
  run $S/crispri_analyse.py
fi
if stage simulations; then
  run $S/theory_sim.py
  run $S/decay_recovery.py
  run $S/gc_correlated_sim.py
fi
if stage liver; then                                     # needs LIVER_SPECTRA (see data/MANIFEST.md)
  run $S/pilot_spectra.py; run $S/pilot_domain_tests.py; run $S/pilot_liver_gc.py
fi
if stage figures; then
  for f in fig1 fig2 fig3 fig4 fig5 fig6 ed1_3 ed4; do run figures/code/$f.py; done
fi
