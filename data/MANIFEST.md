# Inputs

Every script reads raw inputs from `$POSLAYERS_DATA` (default `data/`) and writes every table it produces to
`$POSLAYERS_RESULTS` (default `results/`). Figures go to `$POSLAYERS_FIGS` (default `figures/output/`). Nothing is read from
the current directory. The layout below is what the code expects; file names are those of the public releases.

| Path under `$POSLAYERS_DATA` | Content | Source |
| --- | --- | --- |
| `gtex/gene_reads_adult_gtex_v11_<tissue>_gct.gz` | Gene read counts, 36 tissues (whole blood, oesophagus muscularis and pituitary from v10 as `gene_reads_v10_<tissue>_gct.gz`) | GTEx portal, v11/v10 |
| `GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt` | Sample attributes (RIN, ischaemic time, batches) | GTEx portal, v11 |
| `gtex_eqtl/<Tissue>_v11_eQTLs_signif_pairs.parquet` | Significant cis-eQTL pairs | GTEx portal, v11 |
| `gtex_eqtl/<Tissue>_v11_eQTLs_SuSiE_summary.parquet` | SuSiE credible sets | GTEx portal, v11 |
| `gtex_eqtl/<Tissue>_v11_normalized_expression_bed.gz` | Normalised (inverse-normal) expression used for eQTL mapping | GTEx portal, v11 |
| `biomart_GRCh38_gene_gc.txt` | Ensembl 100 BioMart export: `Gene stable ID`, `Gene % GC content`, `Gene name` (tab-separated) | Ensembl BioMart, GRCh38 |
| `TAD-full/` | TAD stability landscapes (20-bin landscapes and boundaries by stability) | github.com/emcarthur/TAD-stability-heritability |
| `hic/genes_hg19_tss.tsv` | Gene TSS table (hg19) | produced by `scripts/make_tss_table.py` |
| `hic/hic_contacts_<cell>.csv.gz`, `hic/hic_expected_<cell>.csv.gz` | KR contacts between promoters within 2 Mb and expected contact by distance, GM12878 and IMR90 | produced by `scripts/hic_extract.py` from Rao et al. 2014 (GSE63525) |
| `tcga/GDC-PANCAN.htseq_counts.tsv.zip`, `tcga/GDC-PANCAN.mutect2_snv.tsv`, `tcga/GDC-PANCAN.gistic.tsv`, `tcga/GDC-PANCAN_cnv.tsv`, `tcga/GDC-PANCAN_basic_phenotype.tsv` | GDC Pan-Cancer release | UCSC Xena |
| `beataml/beataml_waves1to4_counts_dbgap.txt`, `beataml/beataml_waves1to4_sample_mapping.xlsx`, `beataml/beataml_wes_wv1to4_mutations_dbgap.txt`, `beataml/beataml_wv1to4_clinical.xlsx` | BeatAML2 | biodev.github.io/BeatAML2 |
| `hic/GSE63525_GM12878_primary_replicate_HiCCUPS_looplist_with_motifs_txt.gz`, `hic/GSE63525_IMR90_HiCCUPS_looplist_with_motifs_txt.gz` | HiCCUPS loops with CTCF motif orientation (`scripts/ctcf_barrier.py`) | GEO GSE63525 |
| `schmitt/FitHiC_primary_cohort/` | Optional: Fit-Hi-C output of 40-kb tissue Hi-C maps (`schmitt` stage) | GEO GSE87112 |
| `edit/GSE264491_merged_counts.csv.gz` | Review only: edited erythroblasts | GEO GSE264491 |
| `slam/GSM*_tcount.tsv.gz` | Review only: SLAM-seq T>C counts | GEO GSE100708 |
| `perturb/K562_gwps_normalized_bulk_01.h5ad` | Review only: genome-wide Perturb-seq, pseudobulk | gwps.wi.mit.edu |
| `liver_tables/Table_S6c_stage_effect_per_gene_with_without_composition.csv` | Review only: per-gene fibrosis stage effects in the liver biopsy cohort | output of the liver pipeline, github.com/Danpc11/Liver_Spectra |

The three `scripts/pilot_*.py` scripts (Fig. 1d,e, Extended Data Fig. 1a) run inside the liver pipeline: set `LIVER_SPECTRA`
to a checkout of github.com/Danpc11/Liver_Spectra after running it.

Gene annotation (GRCh38 and GRCh37, Ensembl 100) comes from the `pyannotables` package.

## Source data

The tables that the figure scripts read (46 files, 40 MB) will be deposited on Zenodo, together with
Supplementary Tables 1-9. With them, `bash run_pipeline.sh figures` rebuilds every figure without any of the inputs above:
set `POSLAYERS_RESULTS` to the unpacked folder.
