#!/usr/bin/env bash
# update_repo.sh - brings Danpc11/positional-layers from commit fb0c9e1 to version 1.1.0, answering the code and methods
# review point by point (see REVIEW_RESPONSE.md, written by this script).
#
# Run from the root of the repository:   bash update_repo.sh
# It overwrites the files it lists, deletes the five superseded scripts, then checks the result.
# It does not commit or push.
set -euo pipefail
[ -f pyproject.toml ] && [ -d src/poslayers ] && [ -d scripts ] || { echo "Run this from the root of positional-layers."; exit 1; }
grep -q "name = \"poslayers\"" pyproject.toml || { echo "This does not look like positional-layers."; exit 1; }
echo "Updating files..."
rm -f "scripts/beat_test.py"; echo "  removed  scripts/beat_test.py"
rm -f "scripts/cohesin_test.py"; echo "  removed  scripts/cohesin_test.py"
rm -f "scripts/derive2.py"; echo "  removed  scripts/derive2.py"
rm -f "scripts/lorentz.py"; echo "  removed  scripts/lorentz.py"
rm -f "scripts/predict_cal.py"; echo "  removed  scripts/predict_cal.py"
mkdir -p ".github/workflows"
cat > ".github/workflows/ci.yml" << 'POSLAYERS_FILE_000_END'
name: tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.10", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install -e ".[dev]"
      - run: pytest -q
      # Static check that would have caught the fDATA and rng errors of the review: fail on any undefined name.
      # The three liver pilot scripts import an external module with '*', which pyflakes cannot see into.
      - name: undefined names
        run: |
          out=$(python -m pyflakes scripts/*.py figures/code/*.py src/poslayers/*.py | grep -E "undefined name|may be undefined" | grep -v "common" || true)
          if [ -n "$out" ]; then echo "$out"; exit 1; fi
POSLAYERS_FILE_000_END
echo "  updated  .github/workflows/ci.yml"
cat > ".gitignore" << 'POSLAYERS_FILE_001_END'
__pycache__/
*.py[cod]
.ipynb_checkpoints/
.venv/
env/
*.egg-info/
build/
dist/
!data/.gitkeep
figures/output/
*.npz
*.parquet
*.hic
*.csv.gz
.DS_Store
data/*
!data/MANIFEST.md
results/
POSLAYERS_FILE_001_END
echo "  updated  .gitignore"
cat > "README.md" << 'POSLAYERS_FILE_002_END'
# positional-layers

Code for **"Shared cis-regulation drives local gene co-expression across human tissues"**.

Any positional measurement of gene expression separates exactly into a tissue landscape, one or
more technical layers that inherit their positional structure from clustered attributes of the
genome, and a residual covariance between nearby genes. Two of those layers are not genome
organisation. This repository contains the reference implementation of the decomposition and of
the three quantitative laws, the analysis scripts that reproduce every figure, and an interactive
simulator.

**Interactive simulator:** https://Danpc11.github.io/positional-layers/

## Install

```bash
git clone https://github.com/Danpc11/positional-layers
cd positional-layers
pip install -e ".[analysis,dev]"     # the library alone: pip install -e .
```

Python 3.10 or newer. The library needs numpy, scipy, pandas, statsmodels and matplotlib; the `analysis` extra adds what the
scripts import (pyannotables, pyarrow, anndata, patsy, openpyxl, hic-straw).

## The laws in thirty seconds

```python
from poslayers import simulate_genome, periodogram_identity, isochore_law, eqtl_law, saturation_exponent

sim = simulate_genome(decay_len=5, gc_bias_sd=0.25, seed=1)

# 1. The decomposition is an identity, not an approximation
periodogram_identity(sim["X"])["max_relative_error"]      # ~1e-15

# 2. Isochore law: the covariance a per-sample GC bias creates between two genes, from one number per sample.
#    Tested out of sample in 36 tissues (slope from odd chromosomes, covariance on even ones).
isochore_law(var_b=0.06, gc_i=1.2, gc_j=0.9, sd_i=1.0, sd_j=1.0)

# 3. eQTL law: the correlation a shared variant induces, sign included. In GTEx it predicts the sign and rank of
#    coupling (r = 0.58) but not its scale: colocalised pairs share more covariance than their variants carry.
eqtl_law(freq=0.3, beta1=0.4, beta2=-0.5)                 # negative: opposite effects

# 4. Local elasticity of a saturating response, 1 - occupancy. A HYPOTHESIS for the contact exponent: supported by the
#    expression gradient in lymphoblastoid cells, not significant in fibroblasts.
saturation_exponent(contact=1.0, K=1.0)                   # 0.5
```

## Layout

| Path | What it holds |
| --- | --- |
| `src/poslayers/` | Reference implementation: `decompose.py`, `laws.py`, `simulate.py` |
| `scripts/` | Analysis scripts, one per result, and two shared libraries; see `scripts/SCRIPTS.md` |
| `figures/code/` | One script per main and Extended Data figure |
| `docs/` | The interactive simulator served by GitHub Pages |
| `tests/` | Tests of the identity and of each law |
| `data/` | `MANIFEST.md` only. Inputs are public; see below |

## Data and configuration

Every input is public and none is redistributed here. `data/MANIFEST.md` lists each file, where it comes from and where
the code expects it. Three environment variables control every path; nothing is read from the current directory:

```bash
export POSLAYERS_DATA=/path/to/inputs        # default data/
export POSLAYERS_RESULTS=/path/to/results    # default results/   (every table the scripts write)
export POSLAYERS_FIGS=/path/to/figures       # default figures/output/
```

## Reproducing the analyses and figures

```bash
bash run_pipeline.sh            # everything, in dependency order (several hours)
bash run_pipeline.sh tumours    # one stage: acquire, atlas, eqtl, architecture, tumours, perturbations, simulations, liver, figures
```

`scripts/SCRIPTS.md` gives, for every script, its arguments, the files it reads and writes, and the figure panel it
produces. Shared code is in `scripts/lib_tumour.py` and `scripts/lib_tad.py`. The three liver pilot scripts run inside the
liver pipeline (github.com/Danpc11/Liver_Spectra); point `LIVER_SPECTRA` at a checkout after running it.

## Tests

```bash
pytest -q
```

The tests check the decomposition identity to machine precision, the behaviour of each law, the edge cases of the
API (constant GC, chromosomes shorter than the lag), the decay-length estimator, and two limitations kept visible on
purpose: GC correction removes GC-tracking regulation at domain scale, and comparing a power law with a fixed-size hub
cannot identify a saturating response. CI also fails on any undefined name.

## Review

Version 1.1.0 answers an external review of the code and methods point by point; see `REVIEW_RESPONSE.md`.

## Citation

Pérez-Calixto, D. *et al.* Shared cis-regulation drives local gene co-expression across human tissues (2026).

## Licence

MIT for the code. The public datasets keep the licences of their original sources.
POSLAYERS_FILE_002_END
echo "  updated  README.md"
cat > "REVIEW_RESPONSE.md" << 'POSLAYERS_FILE_003_END'
# Response to the code and methods review of commit fb0c9e1

Version 1.1.0 addresses each point. "Fixed" means changed and verified by running the code; "re-analysed" means the
analysis was redone and the manuscript updated with the new numbers; "documented" means the limitation is now stated in the
code, the tests and the paper.

| # | Point | Status | What changed |
| --- | --- | --- | --- |
| 1 | `fDATA` undefined; tissue names not interpolated | Fixed | The error was introduced by our earlier path-rewriting script, which moved the `f` of f-strings outside the literal. All paths now go through `poslayers.config` (`DATA`, `OUTDIR`, `FIGDIR`) with interpolation kept inside the f-string. CI fails on any undefined name (pyflakes). |
| 2 | Paths contradict the instructions; helper code resolved from the current directory; `exec` | Fixed | No absolute or current-directory paths remain. The seven scripts that read another script's source and ran it with `exec` now import `lib_tumour.py`, `lib_tad.py` or `theory_sim.py`. The refactored tumour code reproduces the earlier estimates exactly. Inputs are organised by source (`data/MANIFEST.md`); file names that differed between producer and consumer (Hi-C, gene editing) were unified. |
| 3 | Figure inputs without producers; `rng` undefined | Fixed | Producers added or recovered: `robust_families.py`, `bystander_liver.py`, `hic_extract.py`, `make_tss_table.py`, `boot_hic_summary.py`, `eqtl_law_summary.py`, and the liver pilot scripts (which run inside the liver pipeline, now declared through `LIVER_SPECTRA`). `saturation_k_by_expression.csv`, previously typed from printed output, is replaced by `boot_hic_summary.csv`. Every figure was built from a fresh results folder. |
| 4 | Meta-analysis hard-coded in Fig. 4 | Fixed | `meta_stag2.py` computes it from the exported estimates; the figure reads `stag2_meta.csv`. The SE was 3.81, not the 3.71 typed into the figure. `SCRIPTS.md` no longer attributes it to `replicate.py`. |
| 5 | Constant GC, `gc_slopes` NaN, long lags fail | Fixed | Rank-aware projection; zero slope for constant GC; lags without pairs return NaN; shape checks. Tests added for each case. |
| 6 | Isochore law compares different quantities | Re-analysed | Common normalisation, linear model on both sides, and an out-of-sample test (slope from odd chromosomes, covariance on even chromosomes): r = 0.83-0.88, observed/predicted 0.98-1.04. The earlier metric overstated the GC component by 5-16%. `isochore_law.py`. |
| 7 | Spectra fitted on one scale, AIC on another | Re-analysed | Fitted and compared on log power, weighted by frequencies per bin, AICc. The two-Lorentzian model wins in all six tissues by a wider margin. `spectral_form.py`. |
| 8 | Fitted exponent does not identify occupancy | Documented, re-analysed | k is reported as descriptive. The tertile gradient holds in lymphoblastoid cells (difference 0.61, 0.40-0.80) but not significantly in fibroblasts (0.34, -0.07 to 0.67). A test now shows that power law versus fixed hub cannot identify saturation: simulated heterogeneous saturation is also fitted better by the hub form. |
| 9 | Pair dependence | Re-analysed | Donor x 10-Mb genomic-block bootstrap for the Hi-C results; genomic-block bootstrap for the TAD effects, the eQTL law and the liver bystanders. Pair-level P values are no longer used for inference. |
| 10 | Colocalisation score treated as a probability; lead effect only; scale mismatch | Re-analysed | Called a score; every shared variant's own effect is propagated; coupling measured on GTEx's inverse-normal expression. r = 0.58 (0.55-0.60); slope 1.98 (1.86-2.12), or 1.35 without PC removal, so the relation is monotonic and correctly signed but not calibrated. `eqtl_law.py`. |
| 11 | Tautological "power law beats hub" test | Fixed | Renamed to a parameter-recovery test; added a test that the hub form is recovered when it generates the data and a characterisation test that the comparison does not identify saturation. |
| 12 | Identity does not attribute the residual to cis regulation; GC correction and GC-correlated biology | Documented, re-analysed | Claims limited to covariance compatible with shared regulation. `gc_correlated_sim.py` shows that GC correction removes 48-71% of GC-tracking regulation at 10 genes and 5-15% between neighbours; domain-scale estimates are lower bounds. |
| 13 | Controls not matched | Documented | Stated as a limitation in the paper; matched controls are future work. |
| 14 | Permutation exchangeability | Re-analysed | Freedman-Lane permutation of reduced-model residuals with the HC3 t statistic of the reported coefficient; a primary specification is declared, the rest are sensitivity analyses. STAG2: -13%, P = 0.020. |
| 15 | Dependencies and documentation | Fixed | `[analysis]` extra lists every import; `SCRIPTS.md` gives arguments, inputs and outputs per script; `data/MANIFEST.md` lists every input; `run_pipeline.sh` runs everything in order. |

## How this version was verified

- `pytest`: 20 tests pass, including one per edge case of point 5 and the mechanism tests of point 11.
- `pyflakes`: no undefined names outside the three liver scripts, which import an external module.
- Every figure script was run against a results folder through `POSLAYERS_DATA`, `POSLAYERS_RESULTS` and `POSLAYERS_FIGS`.
- In a fresh results folder, `robust_families.py`, `bystander_liver.py`, `gc_correlated_sim.py`, `cohesin_freedman_lane.py`, `beataml_freedman_lane.py`, `meta_stag2.py`, `boot_hic_summary.py` and `eqtl_law_summary.py` were run on the real data and Figures 2, 4 and 5 were rebuilt from their output. Tables that also existed before were compared and are identical.
- Not re-run end to end here because of run time: the GTEx atlas, the eQTL scans, the Hi-C bootstrap and the CRISPRi scan. Their code paths are covered by the static check and were run individually during the re-analysis.
POSLAYERS_FILE_003_END
echo "  new      REVIEW_RESPONSE.md"
mkdir -p "data"
cat > "data/MANIFEST.md" << 'POSLAYERS_FILE_004_END'
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
| `edit/GSE264491_merged_counts.csv.gz` | Edited erythroblasts | GEO GSE264491 |
| `slam/GSM*_tcount.tsv.gz` | SLAM-seq T>C counts | GEO GSE100708 |
| `perturb/K562_gwps_normalized_bulk_01.h5ad` | Genome-wide Perturb-seq, pseudobulk | gwps.wi.mit.edu |
| `liver_tables/Table_S6c_stage_effect_per_gene_with_without_composition.csv` | Per-gene fibrosis stage effects in the liver biopsy cohort | output of the liver pipeline, github.com/Danpc11/Liver_Spectra |

The three `scripts/pilot_*.py` scripts (Fig. 1c-e, Extended Data Fig. 1a) run inside the liver pipeline: set `LIVER_SPECTRA`
to a checkout of github.com/Danpc11/Liver_Spectra after running it.

Gene annotation (GRCh38 and GRCh37, Ensembl 100) comes from the `pyannotables` package.
POSLAYERS_FILE_004_END
echo "  new      data/MANIFEST.md"
mkdir -p "figures/code"
cat > "figures/code/ed1_3.py" << 'POSLAYERS_FILE_005_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, DATA, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
# ---- Extended Data Fig. 1: robustness of the cis layer
P4 = pd.read_csv(OUTDIR + 'p4_domain_scale_tests.csv'); Lz = pd.read_csv(OUTDIR + 'lorentz_v2.csv'); E = pd.read_csv(A + 'eqtl_cis_test.csv'); F = pd.read_csv(A + 'robust_families.csv')
import pyannotables as pa
BMg = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
Gg = pa.tables()['homo_sapiens-GRCh38-ensembl100']; Gg = Gg[~Gg.index.duplicated()]; Gg = Gg[Gg.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])].join(BMg[['gc', 'type']], how='inner'); Gg = Gg[Gg.type == 'protein_coding'].sort_values(['Chromosome', 'Start'])
lags_ = np.arange(1, 61); ac = [np.mean([np.mean(((g.gc.values - g.gc.mean()) / g.gc.std())[:-L] * ((g.gc.values - g.gc.mean()) / g.gc.std())[L:]) for c, g in Gg.groupby('Chromosome') if len(g) > L + 5]) for L in lags_]
fig, axs = plt.subplots(1, 4, figsize=(W, W * 0.33)); plt.subplots_adjust(wspace=0.75)
ax = axs[0]; x = np.arange(len(P4)); ax.bar(x - 0.2, P4.cis_L1, 0.4, color=OI['blue'], label='adjacent (cis)'); ax.bar(x + 0.2, P4.domain_L10_30, 0.4, color=OI['red'], label='10–30 genes (domain)')
ax.set_xticks(x); ax.set_xticklabels(['raw', '−comp.', '−props', '−5 PCs', '−10 PCs', '−20 PCs'], fontsize=5.4, rotation=45, ha='right'); ax.set_ylabel('Correlation vs permuted order'); ax.set_ylim(0, 0.3); ax.legend(fontsize=5.2, loc='upper right'); lab(ax, 'a')
ax = axs[1]; dA = Lz.AIC_powerlaw - Lz.AIC_lor2; dB = Lz.AIC_lor1 - Lz.AIC_lor2; ax.scatter(dB, dA, s=16, color=OI['purple'])
for t, b, a in zip(Lz.tissue, dB, dA): ax.annotate(tissue_label(t), (b, a), fontsize=5.0, xytext=(3, 2), textcoords='offset points', ha='left' if b < 40 else 'right')
ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5); ax.set_xlabel('ΔAIC, 1 vs 2 Lorentzians'); ax.set_ylabel('ΔAIC, power law vs 2 Lorentzians'); lab(ax, 'b')
ax = axs[2]; ax.scatter(E.both_eGenes, 100 * E.share_of_cis_from_shared_eQTL, s=12, color=OI['orange']); ax.set_xlabel('Pairs with both genes eGenes'); ax.set_ylabel('Cis explained by shared eQTLs (%)')
lab(ax, 'c')
ax = axs[3]; ax.plot(lags_, ac, color=OI['red'], lw=1.1); ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('GC autocorrelation'); ax.text(0.95, 0.9, f'{ac[0]:.2f} at 1 gene\n{ac[9]:.2f} at 10\n{ac[29]:.2f} at 30', transform=ax.transAxes, ha='right', va='top', fontsize=5.6); lab(ax, 'd')
save(fig, 'ExtData_Fig1_robustness')
# ---- Extended Data Fig. 2: tumours and gene dosage
AN = pd.read_csv(OUTDIR + 'aneuploidy_continuousCN.csv'); FR2 = pd.read_csv(OUTDIR + 'friction2_tumour_tad_clustering.csv')
fig, axs = plt.subplots(1, 3, figsize=(W, W * 0.3)); plt.subplots_adjust(wspace=0.5)
p = AN.pivot(index='cohort', columns='expression', values='mean_far'); x = np.arange(len(p)); ax = axs[0]
ax.bar(x - 0.2, p['raw'], 0.4, color=OI['red'], label='uncorrected'); ax.bar(x + 0.2, p['continuous CN corrected'], 0.4, color=OI['blue'], label='copy-number corrected')
ax.set_xticks(x); ax.set_xticklabels(p.index, fontsize=6); ax.set_ylabel('Long-range floor (20–30 genes)'); ax.set_ylim(0, 0.14); ax.legend(fontsize=5.3, loc='upper center', ncol=2); lab(ax, 'a')
q = AN.pivot(index='cohort', columns='expression', values='rho_far_cna'); ax = axs[1]
ax.bar(x - 0.2, q['raw'], 0.4, color=OI['red']); ax.bar(x + 0.2, q['continuous CN corrected'], 0.4, color=OI['blue']); ax.axhline(0, color='k', lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(q.index, fontsize=6); ax.set_ylabel('Spearman ρ, floor vs CNA burden'); lab(ax, 'b')
ax = axs[2]; x = np.arange(len(FR2)); w = 0.2
ax.bar(x - 1.5 * w, FR2.z_normal, w, color='0.7', label='normal'); ax.bar(x - 0.5 * w, FR2.z_tumour, w, color=OI['red'], label='tumour')
ax.bar(x + 0.5 * w, FR2.z_normal_same_exclusion, w, color='0.7', hatch='///', label='normal, CNA genes removed'); ax.bar(x + 1.5 * w, FR2.z_tumour_CNA_genes_removed, w, color=OI['red'], hatch='///', label='tumour, CNA genes removed')
ax.set_xticks(x); ax.set_xticklabels(FR2.cohort, fontsize=6); ax.set_ylabel('TAD clustering of active genes (z)'); ax.legend(fontsize=5, loc='upper right'); lab(ax, 'c')
save(fig, 'ExtData_Fig2_dosage')
# ---- Extended Data Fig. 3: TAD effect across distance, with and without contact
FR1 = pd.read_csv(OUTDIR + 'friction1_tad_vs_contact.csv'); order = ['(25000.0, 50000.0]', '(50000.0, 100000.0]', '(100000.0, 200000.0]', '(200000.0, 500000.0]', '(500000.0, 1000000.0]', '(1000000.0, 2000000.0]']
labels = ['25–50 kb', '50–100 kb', '100–200 kb', '0.2–0.5 Mb', '0.5–1 Mb', '1–2 Mb']
fig, axs = plt.subplots(1, 2, figsize=(W * 0.7, W * 0.3)); plt.subplots_adjust(wspace=0.45)
for cell, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = FR1[FR1.cell == cell].set_index('distance').reindex(order)
    axs[0].plot(range(6), d.r_same / d.r_diff, 'o-', color=col, ms=3.5, lw=1, label=cell)
    axs[1].plot(range(6), d.TAD_effect, 'o-', color=col, ms=3.5, lw=1, label=f'{cell}, distance only'); axs[1].plot(range(6), d.TAD_effect_given_contact, 's--', color=col, ms=3.5, lw=1, mfc='white', label=f'{cell}, + contact')
for ax in axs: ax.set_xticks(range(6)); ax.set_xticklabels(labels, rotation=30, fontsize=5.8)
axs[0].axhline(1.2, color='0.6', lw=0.8, ls=':'); axs[0].text(5.0, 1.0, 'constant ×1.2 (ref. 8)', fontsize=5.3, color='0.4', ha='right', va='top'); axs[0].set_ylim(0.8, None); axs[0].set_ylabel('Coupling, same TAD / different TAD'); axs[0].legend(fontsize=5.6); lab(axs[0], 'a')
axs[1].set_ylabel('Same-TAD effect on coupling'); axs[1].legend(fontsize=5.2); lab(axs[1], 'b')
save(fig, 'ExtData_Fig3_TAD_vs_contact'); print('ED ok')
POSLAYERS_FILE_005_END
echo "  new      figures/code/ed1_3.py"
mkdir -p "figures/code"
cat > "figures/code/ed4.py" << 'POSLAYERS_FILE_006_END'
"""Extended Data Fig. 4: heatmaps summarising the atlas, the laws and the perturbation rules."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
import schem
R = pd.read_csv(A + 'atlas_results.csv'); L = pd.read_csv(OUTDIR + 'isochore_v2.csv'); L = L.rename(columns={**{f'held_pred_L{k}': f'pred_L{k}' for k in (1, 2, 5, 10, 20, 30)}, **{f'held_obs_L{k}': f'obs_gc_component_L{k}' for k in (1, 2, 5, 10, 20, 30)}}); E = pd.read_csv(A + 'eqtl_cis_test.csv')
fig = plt.figure(figsize=(W, W * 0.95)); gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.55, height_ratios=[1, 1.5])

# a: per-tissue layer budget (tissues x layers)
R = R.sort_values('cis_L1_minusGC_tech', ascending=False).reset_index(drop=True)
cols = ['landscape_share', 'eta2_GCslope_batch', 'domain_L10_30_raw', 'domain_minusGC_tech', 'cis_L1_minusGC_tech', 'cis_decay_length_genes']
names = ['landscape\nshare', 'GC batch\nη²', 'domain,\nnaive', 'domain,\n−GC', 'cis,\nadjacent', 'decay\n(genes)']
M = R[cols].copy()
for c in cols: M[c] = (M[c] - M[c].min()) / (M[c].max() - M[c].min())
ax = fig.add_subplot(gs[0, :])
im = ax.imshow(M.T.values, aspect='auto', cmap='magma')
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=5.2)
ax.set_xticks(range(len(R))); ax.set_xticklabels([tissue_label(t) for t in R.tissue], fontsize=4.4, rotation=90)
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.ax.tick_params(labelsize=5); cb.set_label('scaled within row', fontsize=5.2)
for s in ax.spines.values(): s.set_visible(False)
lab(ax, 'a', -0.05)

# b: isochore law, predicted vs observed by distance (tissues x distances), ratio obs/pred
dists = [1, 2, 5, 10, 20, 30]
Lr = L.copy(); Lr['tissue'] = Lr.tissue.map(tissue_label)
Mb = np.column_stack([Lr[f'obs_gc_component_L{d}'] / Lr[f'pred_L{d}'] for d in dists])
order = np.argsort(np.nanmedian(Mb, 1))[::-1]
ax = fig.add_subplot(gs[1, 0])
im = ax.imshow(Mb[order], aspect='auto', cmap='RdBu_r', vmin=0.6, vmax=1.4)
ax.set_xticks(range(len(dists))); ax.set_xticklabels([f'{d}' for d in dists], fontsize=5.4); ax.set_xlabel('Distance (genes)', fontsize=6)
ax.set_yticks(range(len(Lr))); ax.set_yticklabels(Lr.tissue.values[order], fontsize=4.0)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cb.ax.tick_params(labelsize=5); cb.set_label('observed / predicted (held-out)', fontsize=5.2)
ax.set_title('Isochore law holds at every distance', fontsize=6, loc='left', pad=4)
for s in ax.spines.values(): s.set_visible(False)
lab(ax, 'b', -0.42)

# c: intervention rules as a rule table heatmap
rules = ['shared element\n(HBG1/2 edit)', 'shared enhancer\n(BET inhibitor)', 'shared variant\n(eQTL)', 'disease programme\n(fibrosis)', 'promoter silencing\n(CRISPRi)', 'protein, trans\n(BCL11A edit)']
obs = ['propagates with\ncoupling', 'distance only\n(<50 kb)', 'no cis effect']
G = np.array([[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
ax = fig.add_subplot(gs[1, 1])
im = ax.imshow(G, aspect='auto', cmap='Blues', vmin=0, vmax=1.6)
ax.set_xticks(range(3)); ax.set_xticklabels(obs, fontsize=5.0); ax.set_yticks(range(len(rules))); ax.set_yticklabels(rules, fontsize=5.0)
for i in range(G.shape[0]):
    j = int(np.argmax(G[i])); ax.plot(j, i, marker='o', ms=5, mfc=OI['blue'], mec='white', mew=0.8)
ax.set_xlabel('Observed response of neighbours', fontsize=6); ax.set_ylabel('Intervention acts on', fontsize=6, labelpad=26)
for s in ax.spines.values(): s.set_visible(False)
lab(ax, 'c', -0.42)
save(fig, 'ExtData_Fig4_heatmaps'); print('ed4 ok')
POSLAYERS_FILE_006_END
echo "  updated  figures/code/ed4.py"
mkdir -p "figures/code"
cat > "figures/code/fig1.py" << 'POSLAYERS_FILE_007_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, DATA, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys
import pyannotables as pa
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from theory_sim import simulate
R = pd.read_csv(A + 'atlas_results.csv'); P1 = pd.read_csv(OUTDIR + 'p1_spectra.csv')
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])].join(BM[['gc', 'type']], how='inner'); G = G[G.type == 'protein_coding'].sort_values(['Chromosome', 'Start'])
lags = np.arange(1, 61); ac = []
for L_ in lags:
    v = []
    for c, g in G.groupby('Chromosome'):
        x = (g.gc.values - g.gc.mean()) / g.gc.std()
        if len(x) > L_ + 5: v.append(np.mean(x[:-L_] * x[L_:]))
    ac.append(np.mean(v))
R = pd.read_csv(A + 'atlas_results.csv'); L = pd.read_csv(OUTDIR + 'isochore_v2.csv'); L = L.rename(columns={**{f'held_pred_L{k}': f'pred_L{k}' for k in (1, 2, 5, 10, 20, 30)}, **{f'held_obs_L{k}': f'obs_gc_component_L{k}' for k in (1, 2, 5, 10, 20, 30)}}); P4 = pd.read_csv(OUTDIR + 'p4_domain_scale_tests.csv'); P5 = pd.read_csv(OUTDIR + 'p5a_liver_gc.csv'); S = pd.read_csv(OUTDIR + 'sim_decay_recovery.csv'); P1 = pd.read_csv(OUTDIR + 'p1_spectra.csv')
import schem
fig = plt.figure(figsize=(W, W * 0.85)); gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.5)
# B identity in simulation
Xs, chrs, GC, mu_, z = simulate(seed=1, land_sd=4); c0 = chrs == 0; Xc = Xs[:, c0]; nn = Xc.shape[1]; mbar = Xc.mean(0); Dv = Xc - mbar
lhs = np.mean([np.abs(np.fft.fft(x)) ** 2 for x in Xc], 0); Cl = np.array([np.mean(np.sum(Dv * np.roll(Dv, -L, axis=1), 1)) for L in range(nn)])
land = np.abs(np.fft.fft(mbar)) ** 2; cov = np.real(np.fft.fft(Cl)); k = slice(1, nn // 2)
ax = fig.add_subplot(gs[0, 0:2]); schem.decomposition(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 2]); ax.loglog(lhs[k], (land + cov)[k], '.', ms=1.4, color='0.35', rasterized=True); lim = [lhs[k].min(), lhs[k].max()]; ax.plot(lim, lim, color=OI['red'], lw=0.7)
ax.set_xlabel('Mean periodogram'); ax.set_ylabel('|M(f)|² + S(f)'); ax.text(0.95, 0.06, f'landscape {np.sum(land) / np.sum(lhs):.0%}', transform=ax.transAxes, fontsize=6, ha='right'); lab(ax, 'b')
# D universal peaks real vs random gene order
ax = fig.add_subplot(gs[1, 0]); ax.scatter(R.universal_peaks_random_order, R.universal_peaks_real_order, s=10, color=OI['blue'], lw=0)
m = max(R.universal_peaks_real_order.max(), R.universal_peaks_random_order.max()) * 1.05; ax.plot([0, m], [0, m], 'k:', lw=0.7)
ax.set_xlabel('Peaks, random gene order'); ax.set_ylabel('Peaks, real gene order'); ax.text(0.05, 0.9, f'median ratio {np.median(R.universal_peaks_real_order / R.universal_peaks_random_order):.2f}', transform=ax.transAxes, fontsize=6); lab(ax, 'c')
# E liver: consensus spectrum vs landscape spectrum
ax = fig.add_subplot(gs[1, 1]); ax.loglog(P1.w_static, P1.cons_full, '.', ms=1, color='0.45', rasterized=True); u = P1.universal.fillna(False).astype(bool)
ax.loglog(P1.w_static[u], P1.cons_full[u], '.', ms=3, color=OI['red'], label='186 universal peaks')
ax.set_xlabel('Landscape spectrum'); ax.set_ylabel('Consensus spectrum'); ax.legend(loc='upper left', bbox_to_anchor=(0, 0.92), markerscale=2)
ax.text(0.05, 0.9, f"r = {np.corrcoef(np.log(P1.w_static), np.log(P1.cons_full))[0, 1]:.3f}", transform=ax.transAxes, fontsize=6); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 2]); lbl = ['none', '−comp.', '−GC', '−5 PCs']
cis = [P4.cis_L1.iloc[0], P5.cis_L1.iloc[0], P5.cis_L1.iloc[1], P4.cis_L1.iloc[3]]; dom = [P4.domain_L10_30.iloc[0], P5.domain_L10_30.iloc[0], P5.domain_L10_30.iloc[1], P4.domain_L10_30.iloc[3]]
x = np.arange(4); ax.bar(x - 0.2, cis, 0.4, color=OI['blue'], label='adjacent genes (cis)'); ax.bar(x + 0.2, dom, 0.4, color=OI['red'], label='10–30 genes (domain)')
ax.set_xticks(x); ax.set_xticklabels(lbl, fontsize=5.6, rotation=25); ax.set_ylabel('Correlation vs permuted order'); ax.set_ylim(0, 0.34); ax.legend(loc='upper right', fontsize=5.4); lab(ax, 'e')
ax = fig.add_subplot(gs[2, 0])
for _, r in R.iterrows(): ax.plot([0, 1], [r.domain_L10_30_raw, r.domain_minusGC], color='0.75', lw=0.5)
ax.plot([0, 1], [R.domain_L10_30_raw.median(), R.domain_minusGC.median()], 'o-', color=OI['red'], lw=1.4, ms=4, label='median, 36 tissues')
ax.set_xticks([0, 1]); ax.set_xticklabels(['Naive', 'GC-corrected']); ax.set_xlim(-0.3, 1.3); ax.set_ylabel('Domain correlation (10–30 genes)'); ax.legend(loc='lower left', fontsize=5.4)
lab(ax, 'f')
ax = fig.add_subplot(gs[2, 1])
for lg, col in [(1, OI['blue']), (5, OI['green']), (20, OI['red'])]:
    ax.scatter(L[f'pred_L{lg}'], L[f'obs_gc_component_L{lg}'], s=9, color=col, lw=0, label=f'{lg} gene' + ('s' if lg > 1 else '') + f'  (r = {np.corrcoef(L[f"pred_L{lg}"], L[f"obs_gc_component_L{lg}"])[0, 1]:.2f})')
m = max(L[[c for c in L.columns if c.startswith('obs')]].max().max(), L[[c for c in L.columns if c.startswith('pred')]].max().max()) * 1.05
ax.plot([0, m], [0, m], 'k:', lw=0.7); ax.set_xlabel('Predicted, GC slope from odd chromosomes'); ax.set_ylabel('Observed, even chromosomes'); ax.legend(loc='upper left', fontsize=5.5)
lab(ax, 'g')
ax = fig.add_subplot(gs[2, 2]); g = S.groupby('true_lambda')[['naive', 'corrected']].median()
ax.plot(g.index, g.naive, 'o-', color=OI['red'], ms=3, lw=1, label='naive'); ax.plot(g.index, g.corrected, 's-', color=OI['blue'], ms=3, lw=1, label='GC-corrected')
ax.plot([1, 13], [1, 13], 'k:', lw=0.7, label='identity'); ax.set_xlim(0, 13.5); ax.set_ylim(0, 17); ax.set_xlabel('True cis length (genes)'); ax.set_ylabel('Estimated length (genes)'); ax.legend(loc='lower right', fontsize=5.6)
lab(ax, 'h')
save(fig, 'Fig1_nonpositional_layers'); print('Fig1_nonpositional_layers ok')
POSLAYERS_FILE_007_END
echo "  updated  figures/code/fig1.py"
mkdir -p "figures/code"
cat > "figures/code/fig2.py" << 'POSLAYERS_FILE_008_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
R = pd.read_csv(A + 'atlas_results.csv'); O = pd.read_csv(A + 'orientation_by_tissue.csv'); F = pd.read_csv(A + 'robust_families.csv'); Lz = pd.read_csv(OUTDIR + 'lorentz_v2.csv')
import schem
fig = plt.figure(figsize=(W, W * 0.55)); gs = fig.add_gridspec(2, 4, hspace=0.55, wspace=0.65)

ax = fig.add_subplot(gs[0, 0:2]); schem.two_scales(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 2]); s = R.sort_values('cis_L1_minusGC_tech'); y = np.arange(len(s)); cl = s.tissue.str.startswith('cells')
ax.barh(y, s.cis_L1_minusGC_tech, color=np.where(cl, OI['orange'], OI['blue']), height=0.75)
ax.set_yticks(y[::3]); ax.set_yticklabels([tissue_label(t) for t in s.tissue][::3], fontsize=5.3); ax.set_xlabel('Correlation of adjacent genes\n(GC- and batch-corrected)')
ax.text(1.0, 1.01, 'orange: cell lines', transform=ax.transAxes, ha='right', va='bottom', fontsize=5.8, color=OI['orange']); lab(ax, 'b', -0.62)
ax = fig.add_subplot(gs[0, 3]); O['excess'] = O.mean_r - O.random_pair_floor
order = ['0-1kb', '1-5kb', '5-20kb', '20-100kb', '100-500kb', '>500kb']; mids = [0.5, 3, 12, 50, 250, 1000]
for o, col in [('divergent', OI['green']), ('tandem', OI['blue']), ('convergent', OI['red'])]:
    m = O[O.orientation == o].groupby('dist_bin').excess.median().reindex(order); ax.plot(mids, m.values, 'o-', color=col, ms=3, lw=1, label=o)
ax.set_xscale('log'); ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Correlation above random pairs'); ax.legend(loc='upper right'); lab(ax, 'c')
ax = fig.add_subplot(gs[1, 0]); cats = [('r_same_stem', 'same gene family'), ('r_readthrough', 'read-through'), ('r_clusters', 'classic clusters'), ('r_all', 'all pairs'), ('r_clean', 'clean pairs')]
med = [F[c].median() for c, _ in cats]; lo = [F[c].quantile(0.25) for c, _ in cats]; hi = [F[c].quantile(0.75) for c, _ in cats]
cols = ['0.6', '0.6', '0.6', OI['blue'], OI['green']]; x = np.arange(len(cats)); ax.bar(x, med, color=cols, width=0.65); ax.errorbar(x, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)], fmt='none', ecolor='k', lw=0.7, capsize=2)
ax.set_xticks(x); ax.set_xticklabels([n for _, n in cats], rotation=30, ha='right', fontsize=6); ax.set_ylabel('Adjacent-pair correlation'); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 1]); sp = np.load(OUTDIR + 'spec_v2_thyroid.npy'); ax.loglog(sp[0], sp[1], 'o', ms=2.5, color='0.25', label='thyroid, 684 samples'); ax.loglog(sp[0], sp[2], color=OI['blue'], lw=1.2, label='two Lorentzians')
ax.loglog(sp[0], sp[3], color=OI['red'], lw=0.9, ls='--', label='power law'); ax.set_xlabel('Spatial frequency (cycles per gene)'); ax.set_ylabel('Covariance spectrum S(f)'); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize=5.2, columnspacing=0.8)
lab(ax, 'e')
ax = fig.add_subplot(gs[1, 2:4]); xs = np.arange(len(Lz))
ax.scatter(xs - 0.12, Lz.lam_acf_short, s=14, color=OI['blue'], label='short, autocorrelation'); ax.scatter(xs + 0.12, Lz.lam_spec_short, s=14, marker='s', color=OI['sky'], label='short, spectrum')
ax.scatter(xs - 0.12, Lz.lam_acf_long, s=14, color=OI['red'], label='long, autocorrelation'); ax.set_yscale('log')
ax.set_xticks(xs); ax.set_xticklabels([tissue_label(t) for t in Lz.tissue], rotation=35, ha='right', fontsize=5.8); ax.set_ylabel('Decay length (genes)'); ax.set_ylim(0.6, 300); ax.legend(loc='upper right', fontsize=5.4, ncol=3)
lab(ax, 'f')
save(fig, 'Fig2_cis_layer'); print('Fig2_cis_layer ok')
POSLAYERS_FILE_008_END
echo "  updated  figures/code/fig2.py"
mkdir -p "figures/code"
cat > "figures/code/fig3.py" << 'POSLAYERS_FILE_009_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys
E = pd.read_csv(A + 'eqtl_cis_test.csv'); Dz = pd.read_csv(A + 'coloc_dose_response.csv'); SP = pd.read_csv(A + 'eqtl_tissue_specificity.csv'); Cc = pd.read_csv(A + 'coloc_cis_test.csv')
PC = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv').rename(columns={'pred_r': 'pred_genetic_r', 'obs_r_int_pc15': 'obs_r'})
import schem
rng = np.random.default_rng(0)
fig = plt.figure(figsize=(W, W * 0.55)); gs = fig.add_gridspec(2, 4, hspace=0.6, wspace=0.65)

ax = fig.add_subplot(gs[0, 0:2]); schem.shared_variant(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 2]); cats = [('r_not_both_eGenes', 'not both\neGenes', '0.6'), ('r_eGenes_no_share', 'no shared\nvariant', OI['sky']), ('r_share_same', 'shared,\nsame sign', OI['blue'])]
for i, (c, n, col) in enumerate(cats): ax.scatter(i + rng.uniform(-0.15, 0.15, len(E)), E[c], s=6, color=col, lw=0); ax.hlines(E[c].median(), i - 0.25, i + 0.25, color='k', lw=1.2)
ax.set_xticks(range(3)); ax.set_xticklabels([n for _, n, _ in cats], fontsize=5.4); ax.set_ylabel('Adjacent-pair correlation'); lab(ax, 'b')
ax = fig.add_subplot(gs[0, 3]); s = E.sort_values('delta_same_adj'); y = np.arange(len(s))
ax.scatter(s.delta_same_adj, y, s=8, color=OI['blue'], label='same direction'); ax.scatter(s.delta_opposite_only_adj, y, s=8, color=OI['red'], label='opposite direction only')
ax.axvline(0, color='k', lw=0.6); ax.set_yticks([]); ax.set_ylabel('36 tissues'); ax.set_xlabel('Effect on coupling (adjusted)')
lab(ax, 'c', -0.08)
ax = fig.add_subplot(gs[1, 0]); order = ['0', '0-0.1', '0.1-0.5', '0.5-0.8', '>0.8']; g = Dz.groupby('p_coloc_bin').mean_r
med = [g.get_group(b).median() for b in order]; q1 = [g.get_group(b).quantile(0.25) for b in order]; q3 = [g.get_group(b).quantile(0.75) for b in order]
ax.errorbar(range(5), med, yerr=[np.subtract(med, q1), np.subtract(q3, med)], fmt='o-', color=OI['purple'], ms=4, capsize=2, lw=1.2)
ax.set_xticks(range(5)); ax.set_xticklabels(['none', '<0.1', '0.1–0.5', '0.5–0.8', '>0.8'], fontsize=5.4, rotation=30); ax.set_xlabel('Colocalisation score'); ax.set_ylabel('Adjacent-pair correlation')
lab(ax, 'd')
ax = fig.add_subplot(gs[1, 1]); PC['bin'] = pd.qcut(PC.pred_genetic_r, 8, duplicates='drop'); b = PC.groupby('bin', observed=True).agg(p=('pred_genetic_r', 'mean'), o=('obs_r', 'mean'), se=('obs_r', lambda x: x.std() / np.sqrt(len(x))))
ax.errorbar(b.p, b.o, yerr=1.96 * b.se, fmt='o', color=OI['blue'], ms=4, capsize=2); sl = np.polyfit(PC.pred_genetic_r, PC.obs_r, 1); xx = np.linspace(b.p.min(), b.p.max(), 10)
ax.set_xscale('symlog', linthresh=0.01); ax.set_xlabel('Predicted genetic correlation, Σ w·2p(1−p)β₁β₂'); ax.set_ylabel('Observed coupling (normalized expr.)')
ax.text(0.05, 0.9, f'slope {sl[0]:.2f}; r = {np.corrcoef(PC.pred_genetic_r, PC.obs_r)[0, 1]:.2f}', transform=ax.transAxes, fontsize=6); lab(ax, 'e')
ax = fig.add_subplot(gs[1, 2:4]); ax.scatter(SP.b_other, SP.b_own, s=4, color=OI['green'], alpha=0.5, lw=0); m = max(SP.b_own.max(), SP.b_other.max()); mn = min(SP.b_own.min(), SP.b_other.min())
ax.plot([mn, m], [mn, m], 'k:', lw=0.7); ax.set_xlabel('Effect of sharing in another tissue'); ax.set_ylabel('Effect of sharing in the same tissue')
lab(ax, 'f')
save(fig, 'Fig3_genetic_component'); print('Fig3_genetic_component ok')
POSLAYERS_FILE_009_END
echo "  updated  figures/code/fig3.py"
mkdir -p "figures/code"
cat > "figures/code/fig4.py" << 'POSLAYERS_FILE_010_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
T = pd.read_csv(A + 'tad_cis_test.csv'); TT = pd.read_csv(OUTDIR + 'tad_tumour_results.csv'); CO = pd.read_csv(OUTDIR + 'cohesin_freedman_lane.csv'); MT = pd.read_csv(OUTDIR + 'stag2_meta.csv').iloc[0]; TB = pd.read_csv(OUTDIR + 'tad_block_bootstrap.csv').set_index('tissue')
BA = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv'); AN = pd.read_csv(OUTDIR + 'aneuploidy_continuousCN.csv')
H = {c: pd.read_csv(A + f'hic_coupling_{c}.csv.gz') for c in ['GM12878', 'IMR90']}
import schem
fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 4, hspace=0.75, wspace=0.9)

ax = fig.add_subplot(gs[0, 0]); schem.saturation(ax); lab(ax, 'a', -0.25)
ax = fig.add_subplot(gs[0, 1]); s = T.sort_values('b_same_tad'); y = np.arange(len(s))
tb = TB.reindex(s.tissue); ax.errorbar(s.b_same_tad, y + 0.15, xerr=[s.b_same_tad.values - tb.ci_low.values, tb.ci_high.values - s.b_same_tad.values], fmt='o', ms=3, color=OI['blue'], lw=0.6, capsize=0, label='all pairs (95% block CI)'); ax.scatter(s.b_same_tad_no_shared_eQTL, y - 0.15, s=12, color=OI['sky'], marker='s', label='no shared eQTL')
ax.axvline(0, color='k', lw=0.6); ax.set_yticks(y); ax.set_yticklabels([tissue_label(t) for t in s.tissue], fontsize=5.3); ax.set_xlabel('Same-TAD increase in coupling\n(at equal distance)'); ax.set_xlim(-0.01, 0.11); ax.legend(loc='lower left', bbox_to_anchor=(-0.05, 1.0), fontsize=5.0, ncol=1, handletextpad=0.2, borderaxespad=0)
lab(ax, 'b', -0.6)
bins = [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6]; mids = [37.5, 75, 150, 350, 750, 1500]
ax = fig.add_subplot(gs[0, 2]); kk = {}
for c, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = H[c].copy(); d['db'] = pd.cut(d.tss_distance, bins); lo, hi, mr, mc = [], [], [], []
    for b, g in d.groupby('db', observed=True):
        q = pd.qcut(g.oe.rank(method='first'), 5, labels=False); lo.append(g.r[q == 0].mean()); hi.append(g.r[q == 4].mean()); mr.append(g.r.mean()); mc.append(g.contact_KR.mean())
    ax.plot(mids, hi, 'o-', color=col, ms=3, lw=1.1, label=f'{c}, most contact'); ax.plot(mids, lo, 'o--', color=col, ms=3, lw=0.9, mfc='white', label=f'{c}, least contact'); kk[c] = (mc, mr)
ax.set_xscale('log'); ax.set_xlabel('Distance between promoters (kb)'); ax.set_ylabel('Coupling'); ax.set_ylim(-0.005, 0.16); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=1, fontsize=5); lab(ax, 'c')
ax = fig.add_subplot(gs[0, 3])
for c, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    mc, mr = kk[c]; k = np.polyfit(np.log10(mc), np.log10(np.clip(mr, 1e-4, None)), 1)[0]; ax.loglog(mc, mr, 'o-', color=col, ms=3.5, lw=1, label=f'{c}: k = {k:.2f}')
ax.set_xlabel('Mean Hi-C contact (KR)'); ax.set_ylabel('Mean coupling'); ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), fontsize=5.0, borderaxespad=0); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 0]); r = TT.drop_duplicates('cohort')
x = np.arange(len(r)); ax.bar(x - 0.2, r.within_wt, 0.4, color=OI['green'], label='same TAD'); ax.bar(x + 0.2, r.stable_wt, 0.4, color='0.6', label='across stable boundary')
ax.set_xticks(x); ax.set_xticklabels(['Bladder', 'Endometrium'], fontsize=6.3); ax.set_ylabel('Cis excess (adjacent genes)'); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), fontsize=5.2, ncol=1); ax.set_ylim(0, 0.22)
lab(ax, 'e', -0.3)
ax = fig.add_subplot(gs[1, 1]); BS = pd.read_csv(OUTDIR + 'boot_hic_summary.csv').set_index(['cell', 'stat']); xk = np.arange(3)
for cell, col, off in [('GM12878', OI['blue'], -0.15), ('IMR90', OI['red'], 0.15)]:
    d = BS.loc[cell].loc[['k_t1', 'k_t2', 'k_t3']]; ax.errorbar(xk + off, d.estimate, yerr=[d.estimate - d.ci_low, d.ci_high - d.estimate], fmt='o', color=col, ms=4, capsize=2, label=cell)
ax.set_xticks(xk); ax.set_xticklabels(['low', 'mid', 'high']); ax.set_xlabel('Expression of the pair (tertile)'); ax.set_ylabel('Exponent k in coupling ~ contact^k'); ax.set_ylim(0, 1.05); ax.legend(loc='lower left', fontsize=5.4)
lab(ax, 'f')
ax = fig.add_subplot(gs[1, 2:4]); rows = []
for _, q in CO.iterrows(): rows.append((f"{'STAG2' if q.gene == 'STAG2' else 'CTCF'} {'BLCA' if q.cohort == 'BLCA' else 'UCEC'} {'trunc' if q['class'] == 'truncating' else 'any'}{', strict' if q.tmb_q == 0.7 else ''}", q.adj_pct, q.ci_low, q.ci_high, OI['red'] if q.gene == 'STAG2' else '0.5'))
for _, q in BA.iterrows(): rows.append((f"{q.group.replace(', ', ' ').replace('any coding', 'any').replace('truncating', 'trunc')} AML", q.adj_pct, q.ci_low, q.ci_high, OI['orange']))
rows.append(('STAG2 BLCA + AML, meta', MT.meta_pct, MT.ci_low, MT.ci_high, 'k'))
for i, (n, e, l, h, col) in enumerate(rows[::-1]): ax.errorbar(e, i, xerr=[[e - l], [h - e]], fmt='o', color=col, ms=3, capsize=1.5, lw=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([n for n, *_ in rows[::-1]], fontsize=5.2); ax.yaxis.tick_right(); ax.spines['right'].set_visible(True); ax.spines['left'].set_visible(False); ax.axvline(0, color='k', lw=0.6); ax.set_xlabel('Change in cis excess (%), adjusted')
lab(ax, 'g', -0.12)
save(fig, 'Fig4_architecture_cohesin'); print('Fig4_architecture_cohesin ok')
POSLAYERS_FILE_010_END
echo "  updated  figures/code/fig4.py"
mkdir -p "figures/code"
cat > "figures/code/fig5.py" << 'POSLAYERS_FILE_011_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys
import pyannotables as pa
from scipy import stats
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
BY = pd.read_csv(OUTDIR + 'bystander_liver_t_stage.csv'); ED = pd.read_csv(OUTDIR + 'edit_effects.csv', index_col=0)
WB = pd.read_csv(A + 'pairs_whole_blood.csv.gz'); DR = pd.read_csv(OUTDIR + 'drug_cis_results_thr3.csv'); CR = pd.read_csv(OUTDIR + 'crispri_pairs.csv.gz')
import schem
fig = plt.figure(figsize=(W * 0.92, W * 0.9)); gs = fig.add_gridspec(3, 6, hspace=0.6, wspace=1.4)
def quint(d, a, b):
    q = pd.qcut(d.r, 5, labels=False); return [d.r[q == i].mean() for i in range(5)], [stats.pearsonr(d[a][q == i], d[b][q == i])[0] for i in range(5)]
ax = fig.add_subplot(gs[0, 0:3]); schem.perturbation_rules(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 3:6]); x, y = quint(BY, 't1', 't2'); ax.plot(x, y, 'o-', color=OI['purple'], ms=4, lw=1.2); ax.plot([-0.3, 0.6], [-0.3, 0.6], 'k:', lw=0.6)
ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Coupling in healthy liver (GTEx)'); ax.set_ylabel('Concordance of fibrosis effects\n(437 biopsies)'); lab(ax, 'b')
ax = fig.add_subplot(gs[1, 0:2]); pos = ((G.Start + G.End) / 2).reindex(ED.index); ch = G.Chromosome.astype(str).reindex(ED.index)
for (c, p0, col, tcol, nm) in [('11', 5.25e6, OI['red'], 't_HBG', 'HBG1/2 promoter edit'), ('2', 60.45e6, OI['blue'], 't_BCL11A', 'BCL11A enhancer edit (Casgevy)')]:
    w = (ch == c) & ((pos - p0).abs() < 1.5e6); xx = (pos[w] - p0) / 1e6; yy = ED.loc[w, tcol]
    ax.scatter(xx, yy, s=7, color=col, alpha=0.75, lw=0, label=nm)
    for gid in ED.index[w][np.abs(yy.values) > 4]: ax.annotate(ED.at[gid, 'sym'], (xx[gid], yy[gid]), fontsize=4.8, xytext=(2, 1), textcoords='offset points', color=col)
ax.axhline(0, color='0.7', lw=0.5); ax.set_xlabel('Distance from edited site (Mb)'); ax.set_ylabel('Response (t)'); ax.legend(loc='upper right', fontsize=5.2); lab(ax, 'c')
ax = fig.add_subplot(gs[1, 2:4]); P = WB[WB.dist > 0]
for tcol, col, nm in [('t_HBG', OI['red'], 'HBG1/2'), ('t_BCL11A', OI['blue'], 'BCL11A')]:
    d = P.assign(t1=ED[tcol].reindex(P.g1).values, t2=ED[tcol].reindex(P.g2).values).dropna(); x, y = quint(d, 't1', 't2'); ax.plot(x, y, 'o-', color=col, ms=4, lw=1.2, label=nm)
ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Coupling in healthy blood (GTEx)'); ax.set_ylabel('Concordance of edit responses'); ax.legend(loc='upper left'); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 4:6]); CR = CR[np.isfinite(CR.resp) & np.isfinite(CR.r)]; CR['resp'] = CR.resp.clip(-20, 20); C2 = CR[(CR.type == 'cis') & (CR.kd > 2)].copy()
C2['db'] = pd.cut(C2.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['<10 kb', '10–50 kb', '50–200 kb', '0.2–0.5 Mb', '0.5–1 Mb']); qs = C2.r.quantile([1 / 3, 2 / 3]).values
C2['rt'] = pd.cut(C2.r, [-1, qs[0], qs[1], 1], labels=['low', 'mid', 'high']); tab = C2.pivot_table(index='db', columns='rt', values='resp', aggfunc='mean', observed=True)
xx = np.arange(len(tab))
for j, (c, col) in enumerate(zip(['low', 'mid', 'high'], ['0.75', OI['sky'], OI['blue']])): ax.bar(xx + (j - 1) * 0.26, tab[c], 0.26, color=col, label=f'{c} coupling')
ax.axhline(0, color='k', lw=0.6); ax.set_xticks(xx); ax.set_xticklabels(tab.index, fontsize=5.6, rotation=25); ax.set_ylabel('Neighbour response (robust z)'); ax.legend(loc='lower right', fontsize=5.6)
lab(ax, 'e', -0.12)
ax = fig.add_subplot(gs[2, 0:4]); DR['grp'] = np.where(DR.perturbation.str.contains('JQ1'), 'BET inhibitor (JQ1)', np.where(DR.perturbation.str.contains('BRD4'), 'BRD4 degradation', np.where(DR['class'] == 'signalling', 'signalling inhibitor', 'CDK9 inhibitor')))
cmap = {'BET inhibitor (JQ1)': OI['green'], 'BRD4 degradation': OI['sky'], 'CDK9 inhibitor': OI['orange'], 'signalling inhibitor': '0.55'}
d = DR.sort_values(['grp', 'slope_on_baseline_coupling']); y = np.arange(len(d))
ax.barh(y, d.slope_on_baseline_coupling, color=[cmap[g] for g in d.grp], height=0.7); ax.axvline(0, color='k', lw=0.6)
ax.set_yticks(y); ax.set_yticklabels([f'{a} · {b}' for a, b in zip(d.cell, d.perturbation)], fontsize=5); ax.set_xlabel('Propagation to coupled neighbours (slope)')
from matplotlib.patches import Patch; ax.legend(handles=[Patch(color=v, label=k) for k, v in cmap.items()], loc='upper right', fontsize=5.4); lab(ax, 'f', -0.55)
save(fig, 'Fig5_perturbations'); print('Fig5_perturbations ok')
POSLAYERS_FILE_011_END
echo "  updated  figures/code/fig5.py"
mkdir -p "figures/code"
cat > "figures/code/fig6.py" << 'POSLAYERS_FILE_012_END'
import os
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import W, lab, os, plt, save, sys
import schem
fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.18, height_ratios=[1, 1.25])
ax = fig.add_subplot(gs[0, :]); schem.layers_diagram(ax); lab(ax, 'a', -0.02)
ax = fig.add_subplot(gs[1, 0]); schem.mechanism_diagram(ax); lab(ax, 'b', -0.04)
ax = fig.add_subplot(gs[1, 1]); schem.perturbation_rules(ax); lab(ax, 'c', -0.04)
save(fig, 'Fig6_model'); print('fig6 ok')
POSLAYERS_FILE_012_END
echo "  updated  figures/code/fig6.py"
mkdir -p "figures/code"
cat > "figures/code/style.py" << 'POSLAYERS_FILE_013_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np, pandas as pd
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7.5,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.2, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42, 'savefig.dpi': 300})
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'yellow': '#F0E442', 'grey': '#7F7F7F', 'black': '#000000'}
W = 180 / 25.4
OUT = FIGDIR
def lab(ax, s, dx=-0.16, dy=1.24): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left') if False else ax.text(dx, dy, s.lower(), transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left')
def save(fig, name): fig.savefig(OUT + name + '.pdf', bbox_inches='tight'); fig.savefig(OUT + name + '.png', dpi=220, bbox_inches='tight'); plt.close(fig)
A = OUTDIR
NAMES = {'cells_ebv-transformed_lymphocytes': 'EBV lymphocytes', 'cells_ebv_transformed_lymphocytes': 'EBV lymphocytes', 'cells_cultured_fibroblasts': 'Fibroblasts', 'esophagus_gastroesophageal_junction': 'Oesophagus, GEJ',
         'esophagus_mucosa': 'Oesophagus, mucosa', 'esophagus_muscularis': 'Oesophagus, muscularis', 'skin_sun_exposed_lower_leg': 'Skin', 'brain_frontal_cortex_ba9': 'Brain, frontal cortex',
         'brain_caudate_basal_ganglia': 'Brain, caudate', 'brain_cerebellar_hemisphere': 'Brain, cerebellum', 'small_intestine_terminal_ileum': 'Small intestine', 'adipose_visceral_omentum': 'Adipose, visceral',
         'adipose_subcutaneous': 'Adipose, subcutaneous', 'heart_left_ventricle': 'Heart, ventricle', 'heart_atrial_appendage': 'Heart, atrium', 'colon_transverse': 'Colon, transverse', 'colon_sigmoid': 'Colon, sigmoid'}
def tissue_label(t):
    t = t.replace('-', '_') if t not in NAMES else t
    return NAMES.get(t, NAMES.get(t.replace('_', '-'), t.replace('_', ' ').capitalize()))
POSLAYERS_FILE_013_END
echo "  updated  figures/code/style.py"
cat > "pyproject.toml" << 'POSLAYERS_FILE_014_END'
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "poslayers"
version = "1.1.0"
description = "Decomposition of positional signals in transcriptomic data"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
# the library itself (src/poslayers)
dependencies = ["numpy>=1.26", "scipy>=1.13", "pandas>=2.2", "statsmodels>=0.14", "matplotlib>=3.8"]

[project.optional-dependencies]
# everything the analysis scripts in scripts/ and the figure code in figures/code/ import
analysis = ["pyannotables>=0.5", "pyarrow>=14", "anndata>=0.10", "patsy>=0.5", "openpyxl>=3.1", "hic-straw>=1.3"]
dev = ["pytest>=7.4", "pyflakes>=3.1"]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
POSLAYERS_FILE_014_END
echo "  updated  pyproject.toml"
cat > "run_pipeline.sh" << 'POSLAYERS_FILE_015_END'
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
  run $S/gc_correlated_sim.py
fi
if stage liver; then                                     # needs LIVER_SPECTRA (see data/MANIFEST.md)
  run $S/pilot_spectra.py; run $S/pilot_domain_tests.py; run $S/pilot_liver_gc.py
fi
if stage figures; then
  for f in fig1 fig2 fig3 fig4 fig5 fig6 ed1_3 ed4; do run figures/code/$f.py; done
fi
POSLAYERS_FILE_015_END
echo "  new      run_pipeline.sh"
mkdir -p "scripts"
cat > "scripts/SCRIPTS.md" << 'POSLAYERS_FILE_016_END'
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
| atlas | `robust_families.py` | — | `pairs_*.csv.gz` | `robust_families.csv` | Fig. 2d |
| eqtl | `eqtl_share.py` | signif-pairs parquet files | `gtex_eqtl/` | `shares_<tissue>.csv.gz` | |
| eqtl | `eqtl_test.py` | — | `pairs_*`, `shares_*` | `eqtl_cis_test.csv`, `eqtl_tissue_specificity.csv` | Fig. 3b,c,f; ED 1c |
| eqtl | `coloc.py` | SuSiE parquet files | `gtex_eqtl/` | `coloc_<tissue>.csv.gz` | |
| eqtl | `coloc_test.py` | — | `pairs_*`, `shares_*`, `coloc_*` | `coloc_cis_test.csv`, `coloc_dose_response.csv` | Fig. 3d |
| eqtl | `predict.py` | GTEx count files | `gtex/`, SuSiE | `pred_pairs_<tissue>.csv.gz` (pair list for the eQTL law); also `isochore_law.csv`, the superseded metric kept for comparison | |
| eqtl | `eqtl_law.py` | — | `pred_pairs_*`, `gtex_eqtl/` | `eqtl_law_v2_pairs.csv` | Fig. 3e; ST4B |
| eqtl | `eqtl_law_summary.py` | — | `eqtl_law_v2_pairs.csv` | `eqtl_law_v2_summary.csv` | ST4A |
| architecture | `tad_test.py` | — | `pairs_*`, `shares_*`, `TAD-full/` | `tad_cis_test.csv` | Fig. 4b; ST5 |
| architecture | `tad_bootstrap.py` | — | `pairs_*`, `TAD-full/` | `tad_block_bootstrap.csv` | Fig. 4b; ST5 |
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
POSLAYERS_FILE_016_END
echo "  updated  scripts/SCRIPTS.md"
mkdir -p "scripts"
cat > "scripts/aneuploidy.py" << 'POSLAYERS_FILE_017_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm
from lib_tumour import CHR, prepare

def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
rows = []
for c in ['BLCA', 'UCEC', 'STAD', 'COAD', 'GBM']:
    D, Draw, chrs, samp, cov = prepare(c, codes=('01', '03'), return_raw=True); cna = cov.cna_burden.values
    for name, M in [('raw', Draw), ('copy-number corrected', D)]:
        far = lagscore(M, chrs, (20, 30)); near = lagscore(M, chrs, (1, 2, 3))
        rows.append({'cohort': c, 'expression': name, 'tumours': len(far), 'mean_far': far.mean(), 'mean_near': near.mean(),
                     'rho_far_cna': pd.Series(far).corr(pd.Series(cna), method='spearman'), 'rho_cis_excess_cna': pd.Series(near - far).corr(pd.Series(cna), method='spearman'), 'cohort_mean_cna': cna.mean()})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'aneuploidy_scaling.csv', index=False); pd.set_option('display.width', 220); print(R.round(3).to_string(index=False))
raw = R[R.expression == 'raw']; print('\nacross cohorts (raw): Spearman(mean far floor, mean CNA burden) = %.2f' % raw.mean_far.corr(raw.cohort_mean_cna, method='spearman'))
POSLAYERS_FILE_017_END
echo "  updated  scripts/aneuploidy.py"
mkdir -p "scripts"
cat > "scripts/atlas.py" << 'POSLAYERS_FILE_018_END'
"""Positional-layer atlas across GTEx tissues: landscape, technical isochore (GC) and cis covariance."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd
from scipy import stats
OUT = OUTDIR + 'atlas_results.csv'
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene start (bp)': 'start', 'Chromosome/scaffold name': 'chr'})
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = BM[BM.chr.astype(str).isin(CHR)].drop_duplicates('gid').set_index('gid'); BM['chr'] = BM.chr.astype(str)
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
lags = np.array([1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 150, 200, 300])
def prof(D, chrs, perm=False):
    out = {l: [] for l in lags}
    for c in CHR:
        idx = np.where(chrs == c)[0]
        if perm: idx = np.random.default_rng(3).permutation(idx)
        Z = D[idx]; Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
        for l in lags:
            if len(Z) > l + 5: out[l] += list((Z[:-l] * Z[l:]).mean(1))
    return np.array([np.mean(out[l]) for l in lags])
def excess(D, chrs): return prof(D, chrs) - prof(D, chrs, perm=True)
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
def pgram(x): n = len(x); return np.abs(np.fft.rfft(x)[1:n // 2 + 1]) ** 2
def whiten(P, n):
    k = np.arange(1, len(P) + 1); lx = np.log10(k / n); b = np.polyfit(lx, np.log10(P + 1e-12), 1); w = P / 10 ** (b[0] * lx + b[1]); return w / w.mean()
def n_universal(Yc, chrs, perm_seed=None):
    tot = 0
    for c in CHR:
        idx = np.where(chrs == c)[0]
        if perm_seed is not None: idx = np.random.default_rng(perm_seed).permutation(idx)
        Z = Yc[idx]; n = len(idx); W = np.array([whiten(pgram(Z[:, i]), n) for i in range(Z.shape[1])])
        tot += int((np.mean(W > 3, 0) >= 0.9).sum())
    return tot
done = set(pd.read_csv(OUT).tissue) if os.path.exists(OUT) else set()
for f in sys.argv[1:]:
    t = os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', '')
    if t in done: continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
    if len(samp) < 40: print(t, 'skipped: too few samples', len(samp)); continue
    C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
    g = BM.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); mu = Yc.mean(1); D = Yc - mu[:, None]
    share = np.sum(mu ** 2) * D.shape[1] / (np.sum(mu ** 2) * D.shape[1] + np.sum(D ** 2))
    sub = np.random.default_rng(0).choice(D.shape[1], min(120, D.shape[1]), replace=False)
    up_real = n_universal(Yc[:, sub], chrs); up_perm = n_universal(Yc[:, sub], chrs, perm_seed=11)
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); slope = np.array([np.polyfit(gcz, D[:, i], 1)[0] for i in range(D.shape[1])])
    lab = a.SMNABTCH.astype(str).values; grp = [slope[lab == b] for b in np.unique(lab) if (lab == b).sum() >= 3]
    eta2 = sum(len(x) * (x.mean() - slope.mean()) ** 2 for x in grp) / np.sum((np.concatenate(grp) - slope.mean()) ** 2) if len(grp) > 1 else np.nan
    e_raw = excess(D, chrs); Dg = resid_genes(D, np.column_stack([gcz, gcz ** 2])); e_gc = excess(Dg, chrs)
    Dgt = resid_samples(Dg, tech) if tech.shape[1] < D.shape[1] - 10 else Dg; e_gct = excess(Dgt, chrs)
    m = (lags <= 30) & (e_gct > 0.003); b = np.polyfit(lags[m], np.log(e_gct[m]), 1) if m.sum() >= 3 else [np.nan]
    row = {'tissue': t, 'version': 'v10' if 'v10' in f else 'v11', 'samples': D.shape[1], 'genes': D.shape[0], 'landscape_share': share,
           'universal_peaks_real_order': up_real, 'universal_peaks_random_order': up_perm,
           'rho_GCslope_RIN': stats.spearmanr(slope, rin)[0], 'eta2_GCslope_batch': eta2,
           'cis_L1_raw': e_raw[0], 'domain_L10_30_raw': e_raw[6:10].mean(), 'long_L150_300_raw': e_raw[10:].mean(),
           'cis_L1_minusGC': e_gc[0], 'domain_minusGC': e_gc[6:10].mean(),
           'cis_L1_minusGC_tech': e_gct[0], 'cis_L2_minusGC_tech': e_gct[1], 'cis_L5_minusGC_tech': e_gct[4], 'domain_minusGC_tech': e_gct[6:10].mean(), 'cis_decay_length_genes': -1 / b[0]}
    pd.DataFrame([row]).to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False); print(t, D.shape, 'ok', flush=True)
POSLAYERS_FILE_018_END
echo "  updated  scripts/atlas.py"
mkdir -p "scripts"
cat > "scripts/beataml_freedman_lane.py" << 'POSLAYERS_FILE_019_END'
"""Replication of the cohesin effect on the cis layer in BeatAML2 (open data). Per-sample cis-excess score as in TCGA:
neighbour products at 1-3 genes minus 20-30 genes, on GC-corrected expression. Covariates: blasts, monocytic score, log TMB, sex,
specimen type. 10,000 label permutations on covariate-residualised scores."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.api as sm
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
C = pd.read_csv(DATA + 'beataml/beataml_waves1to4_counts_dbgap.txt', sep='\t'); C['gid'] = C.stable_id.str.split('.').str[0]; sym = dict(zip(C.gid, C.display_label.astype(str)))
C = C.drop_duplicates('gid').set_index('gid')[[x for x in C.columns if x.startswith('BA')]]
mp = pd.read_excel(DATA + 'beataml/beataml_waves1to4_sample_mapping.xlsx'); cl = pd.read_excel(DATA + 'beataml/beataml_wv1to4_clinical.xlsx')
mp = mp[(mp.rna_control != 'yes') & mp.dbgap_rnaseq_sample.notna() & mp.dbgap_dnaseq_sample.notna()] if 'rna_control' in mp else mp
pairs = mp[['dbgap_rnaseq_sample', 'dbgap_dnaseq_sample']].dropna().drop_duplicates('dbgap_rnaseq_sample')
m = pd.read_csv(DATA + 'beataml/beataml_wes_wv1to4_mutations_dbgap.txt', sep='\t', low_memory=False); tmb = m.groupby('dbgap_sample_id').size()
pairs = pairs[pairs.dbgap_rnaseq_sample.isin(C.columns) & pairs.dbgap_dnaseq_sample.isin(tmb.index)]
C = C[pairs.dbgap_rnaseq_sample.values]; C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True)
gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
lag = lambda L: np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0)
e = np.mean([lag(L) for L in (1, 2, 3)], 0) - np.mean([lag(L) for L in (20, 30)], 0)
syms = np.array([sym.get(x, '') for x in g.index]); zsc = lambda gs: np.nanmean([(Y[syms == s][0] - Y[syms == s][0].mean()) / (Y[syms == s][0].std() + 1e-9) for s in gs if (syms == s).any()], 0)
dna = pairs.dbgap_dnaseq_sample.values; rna = pairs.dbgap_rnaseq_sample.values
clin = cl.drop_duplicates('dbgap_rnaseq_sample').set_index('dbgap_rnaseq_sample').reindex(rna)
cov = pd.DataFrame({'monocytic': zsc(['CD14', 'LYZ', 'CSF1R', 'FCGR1A', 'CD68']) - zsc(['CD34', 'KIT', 'PROM1']), 'log_tmb': np.log10(tmb.reindex(dna).values + 1),
                    'blasts': pd.to_numeric(clin['%.Blasts.in.BM'], errors='coerce').fillna(pd.to_numeric(clin['%.Blasts.in.PB'], errors='coerce')).values,
                    'male': (clin.consensus_sex.astype(str).str.lower() == 'male').astype(float).values, 'pb_specimen': clin.specimenType.astype(str).str.contains('Peripheral', case=False).astype(float).values}, index=rna)
cov['blasts'] = cov.blasts.fillna(cov.blasts.median())
COH = ['STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']; TR = ['frameshift_variant', 'stop_gained', 'splice_acceptor_variant', 'splice_donor_variant']
rows = []
for name, genes, cls in [('STAG2, any coding', ['STAG2'], None), ('STAG2, truncating', ['STAG2'], TR), ('cohesin, any coding', COH, None), ('cohesin, truncating', COH, TR)]:
    mm = m[m.symbol.isin(genes)]; mm = mm[mm.variant_classification.isin(cls)] if cls else mm
    mut = pd.Index(dna).isin(set(mm.dbgap_sample_id)).astype(float)
    Xs = (cov - cov.mean()) / cov.std(); Xd = sm.add_constant(pd.DataFrame({'mutant': mut}, index=cov.index).join(Xs))
    fit = sm.OLS(e, Xd).fit(cov_type='HC3'); base = e[mut == 0].mean(); r = e - sm.OLS(e, sm.add_constant(Xs.values)).fit().fittedvalues
    red = sm.OLS(e, sm.add_constant(Xs.values)).fit(); Xf = Xd.values; t_obs = fit.tvalues['mutant']; rng = np.random.default_rng(1)
    tnull = np.array([sm.OLS(red.fittedvalues + rng.permutation(red.resid), Xf).fit(cov_type='HC3').tvalues[1] for _ in range(2000)])
    rows.append({'group': name, 'n_samples': len(e), 'n_mutant': int(mut.sum()), 'cis_excess_wt': base, 'raw_pct': 100 * (e[mut == 1].mean() - base) / base,
                 'adj_pct': 100 * fit.params['mutant'] / base, 'ci_low': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high': 100 * fit.conf_int().loc['mutant', 1] / base, 'se_pct': 100 * fit.bse['mutant'] / base, 'p_HC3': fit.pvalues['mutant'], 'p_freedman_lane': (np.sum(np.abs(tnull) >= abs(t_obs)) + 1) / 2001})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'beataml_freedman_lane.csv', index=False); pd.set_option('display.width', 200); print(R.round(4).to_string(index=False))
POSLAYERS_FILE_019_END
echo "  new      scripts/beataml_freedman_lane.py"
mkdir -p "scripts"
cat > "scripts/boot_hic.py" << 'POSLAYERS_FILE_020_END'
"""Review point 9: donor x genomic-block bootstrap for the Hi-C results (Fig. 4c,d,f).
Each replicate resamples GTEx donors with replacement (re-standardising expression, then recomputing every pair's coupling) AND
resamples 10-Mb genomic blocks with replacement (blocks defined on the first gene's position). Statistics per replicate:
  contact_coef : coefficient of log2 O/E contact in r ~ B-spline(log10 distance, 5 df) + log_oe
  k_distance   : log-log slope of mean coupling on mean contact across the six distance bins (as in Fig. 4d)
  k_contact    : power-law exponent fitted to mean coupling over 20 contact quantiles (as in Supplementary Note 1)
  k_t1..k_t3   : the same exponent within expression tertiles (Fig. 4f)
Resumable: appends one row per replicate to boot_hic_<cell>.csv."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, pyannotables as pa
from scipy.optimize import curve_fit
from patsy import dmatrix
cell, tissue, B = sys.argv[1], sys.argv[2], int(sys.argv[3])
OUT = OUTDIR + f'boot_hic_{cell}.csv'; done = len(pd.read_csv(OUT)) if os.path.exists(OUT) else 0
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
C = pd.read_csv(DATA + f'gtex/gene_reads_adult_gtex_v11_{tissue}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description')
C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]; samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); expr = Y.mean(1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[C.index, 'gc'].values; gz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gz)), gz, gz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T                  # technical residuals, computed once on all donors
pos = {g: i for i, g in enumerate(C.index)}
H = pd.read_csv(OUTDIR + f'hic_coupling_{cell}.csv.gz', dtype={'chr': str}); H = H[H.g1.isin(pos) & H.g2.isin(pos) & (H.contact_KR > 0)].reset_index(drop=True)
i1 = H.g1.map(pos).values; i2 = H.g2.map(pos).values
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
start = G.Start.reindex(H.g1).values; block = (H.chr.astype(str) + ':' + (np.nan_to_num(start) // 1e7).astype(int).astype(str)).values
ub, binv = np.unique(block, return_inverse=True); members = [np.where(binv == k)[0] for k in range(len(ub))]
X0 = np.asarray(dmatrix('bs(x, df=5)', {'x': H.log_d.values}, return_type='dataframe')); Xf = np.column_stack([X0, H.log_oe.values])
dbin = pd.cut(H.tss_distance, [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6], labels=False).values
emin = np.minimum(expr[i1], expr[i2]); tert = pd.qcut(emin, 3, labels=False)
pw = lambda c, r0, A, k: r0 + A * c ** k
def k_fit(c, r, nq):
    q = pd.qcut(np.log10(c), nq, labels=False, duplicates='drop'); g = pd.DataFrame({'c': c, 'r': r, 'q': q}).groupby('q').agg(c=('c', 'median'), r=('r', 'mean'), s=('r', 'sem'))
    try: return curve_fit(pw, g.c, g.r, p0=[0.0, 0.002, 0.5], sigma=g.s, bounds=([-0.05, 0, 0.01], [0.05, 5, 3]), maxfev=20000)[0][2]
    except Exception: return np.nan
def stats_for(r, idx):
    out = {}
    beta = np.linalg.lstsq(Xf[idx], r[idx], rcond=None)[0]; out['contact_coef'] = beta[-1]
    dd = pd.DataFrame({'b': dbin[idx], 'r': r[idx], 'c': H.contact_KR.values[idx], 'd': H.tss_distance.values[idx]}).groupby('b').agg(r=('r', 'mean'), c=('c', 'mean'))
    out['k_distance'] = np.polyfit(np.log10(dd.c), np.log10(dd.r.clip(lower=1e-4)), 1)[0]
    out['k_contact'] = k_fit(H.contact_KR.values[idx], r[idx], 20)
    for t in range(3):
        j = idx[tert[idx] == t]; out[f'k_t{t + 1}'] = k_fit(H.contact_KR.values[j], r[j], 15)
    return out
def coupling(cols):
    M = D[:, cols]; Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9); r = np.empty(len(H))
    for s in range(0, len(H), 150000): r[s:s + 150000] = (Z[i1[s:s + 150000]] * Z[i2[s:s + 150000]]).mean(1)
    return r
rng = np.random.default_rng(2024 + done)
if done == 0:
    est = stats_for(coupling(np.arange(D.shape[1])), np.arange(len(H))); est['rep'] = -1
    pd.DataFrame([est]).to_csv(OUT, index=False); print('point estimate', {k: round(v, 4) for k, v in est.items()}, flush=True); done = 1
for b in range(done, B + 1):
    cols = rng.integers(0, D.shape[1], D.shape[1]); idx = np.concatenate([members[k] for k in rng.integers(0, len(members), len(members))])
    st = stats_for(coupling(cols), idx); st['rep'] = b
    pd.DataFrame([st]).to_csv(OUT, mode='a', header=False, index=False)
    if b % 10 == 0: print(b, {k: round(v, 3) for k, v in st.items()}, flush=True)
POSLAYERS_FILE_020_END
echo "  new      scripts/boot_hic.py"
mkdir -p "scripts"
cat > "scripts/boot_hic_summary.py" << 'POSLAYERS_FILE_021_END'
"""Summarise the donor x genomic-block bootstrap of boot_hic.py into estimates and 95% percentile intervals (Fig. 4d,f; ST8D)."""
import numpy as np, pandas as pd
from poslayers.config import OUTDIR
rows = []
for cell in ('GM12878', 'IMR90'):
    d = pd.read_csv(OUTDIR + f'boot_hic_{cell}.csv'); pt = d[d.rep == -1].iloc[0]; b = d[d.rep >= 0]
    for k in ['contact_coef', 'k_distance', 'k_contact', 'k_t1', 'k_t2', 'k_t3']:
        lo, hi = np.nanpercentile(b[k], [2.5, 97.5])
        rows.append({'cell': cell, 'stat': k, 'estimate': pt[k], 'ci_low': lo, 'ci_high': hi, 'boot_se': b[k].std(), 'B': b[k].notna().sum()})
    g = b.k_t1 - b.k_t3
    rows.append({'cell': cell, 'stat': 'k_t1_minus_k_t3', 'estimate': pt.k_t1 - pt.k_t3, 'ci_low': np.nanpercentile(g, 2.5), 'ci_high': np.nanpercentile(g, 97.5), 'boot_se': g.std(), 'B': g.notna().sum()})
    rows.append({'cell': cell, 'stat': 'P(k_t1>k_t2>k_t3) across replicates', 'estimate': np.mean((b.k_t1 > b.k_t2) & (b.k_t2 > b.k_t3)), 'ci_low': np.nan, 'ci_high': np.nan, 'boot_se': np.nan, 'B': len(b)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'boot_hic_summary.csv', index=False); print(R.round(3).to_string(index=False))
POSLAYERS_FILE_021_END
echo "  new      scripts/boot_hic_summary.py"
mkdir -p "scripts"
cat > "scripts/bystander_liver.py" << 'POSLAYERS_FILE_022_END'
"""Disease bystanders (Fig. 5b): concordance of fibrosis stage effects between adjacent genes versus their baseline coupling in GTEx liver.
Inputs: pairs_liver.csv.gz (atlas.py) and the per-gene stage effects of the liver biopsy cohort (DATA/liver_tables/Table_S6c_...csv).
The 95% interval resamples 10-Mb genomic blocks (2,000 replicates); donor-level data are not available, so it does not reflect donor sampling."""
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.formula.api as smf
from scipy import stats
from poslayers.config import DATA, OUTDIR
P = pd.read_csv(OUTDIR + 'pairs_liver.csv.gz'); P = P[P.dist > 0]
T = pd.read_csv(DATA + 'liver_tables/Table_S6c_stage_effect_per_gene_with_without_composition.csv').set_index('gene_id')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
rows = []
for col in ('t_stage', 't_stage_comp_adjusted'):
    d = P.assign(t1=T[col].reindex(P.g1).values, t2=T[col].reindex(P.g2).values).dropna(subset=['t1', 't2']).reset_index(drop=True)
    d['prod'] = (d.t1 - d.t1.mean()) / d.t1.std() * (d.t2 - d.t2.mean()) / d.t2.std()
    d.to_csv(OUTDIR + f'bystander_liver_{col}.csv', index=False)
    m = smf.ols('prod ~ r + C(orientation) + np.log10(dist)', data=d).fit()
    blk = (G.Chromosome.reindex(d.g1).astype(str).values + ':' + (G.Start.reindex(d.g1).fillna(0).values // 1e7).astype(int).astype(str))
    grp = [np.where(blk == b)[0] for b in np.unique(blk)]; rng = np.random.default_rng(11)
    sl = lambda ix: np.polyfit(d.r.values[ix], d['prod'].values[ix], 1)[0]
    est = sl(np.arange(len(d))); bs = np.array([sl(np.concatenate([grp[k] for k in rng.integers(0, len(grp), len(grp))])) for _ in range(2000)])
    rows.append({'stage_effect': col, 'pairs': len(d), 'blocks': len(grp), 'slope': est, 'slope_adjusted_orientation_distance': m.params['r'],
                 'ci_low': np.percentile(bs, 2.5), 'ci_high': np.percentile(bs, 97.5), 'boot_se': bs.std(), 'p_normal_approx': 2 * stats.norm.sf(abs(est / bs.std()))})
pd.DataFrame(rows).to_csv(OUTDIR + 'bystander_block_bootstrap.csv', index=False); print(pd.DataFrame(rows).round(4).to_string(index=False))
POSLAYERS_FILE_022_END
echo "  new      scripts/bystander_liver.py"
mkdir -p "scripts"
cat > "scripts/cn_continuous.py" << 'POSLAYERS_FILE_023_END'
"""Continuous gene-level copy number from GDC segment log2 ratios, for the cohorts already extracted."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
mid = ((G.Start + G.End) / 2).values; gch = G.Chromosome.astype(str).values; gid = G.index.values
for c in ['BLCA', 'UCEC', 'STAD', 'COAD', 'GBM']:
    z = np.load(OUTDIR + f'expr_{c}.npz', allow_pickle=True); samp = set(z['samples']); genes = set(z['genes'])
    sel = np.isin(gid, list(genes)); gi, gm, gc = gid[sel], mid[sel], gch[sel]
    parts = [ch[ch['sample'].isin(samp)] for ch in pd.read_csv(DATA + 'tcga/GDC-PANCAN_cnv.tsv', sep='\t', chunksize=1_000_000, dtype={'Chrom': str})]
    S = pd.concat(parts); S['Chrom'] = S.Chrom.str.replace('chr', ''); smp = sorted(S['sample'].unique())
    M = np.zeros((len(gi), len(smp)), dtype=np.float32)
    for k, (s, d) in enumerate(S.groupby('sample')):
        for ch, dd in d.groupby('Chrom'):
            idx = np.where(gc == ch)[0]
            if not len(idx): continue
            dd = dd.sort_values('Start'); j = np.searchsorted(dd.Start.values, gm[idx], side='right') - 1
            ok = (j >= 0) & (gm[idx] <= dd.End.values[np.clip(j, 0, None)]); M[idx[ok], smp.index(s)] = dd.value.values[j[ok]]
    np.savez_compressed(OUTDIR + f'cn_cont_{c}.npz', CN=M, genes=gi, samples=np.array(smp)); print(c, M.shape, 'segments', len(S), flush=True)
POSLAYERS_FILE_023_END
echo "  updated  scripts/cn_continuous.py"
mkdir -p "scripts"
cat > "scripts/cohesin_freedman_lane.py" << 'POSLAYERS_FILE_024_END'
"""Review point 14: Freedman-Lane permutation for the cohesin/CTCF effect on cis excess.
Reduced model y ~ covariates gives fitted values and residuals; each permutation shuffles the RESIDUALS, adds them back to the
reduced fit, refits the FULL model y* ~ mutant + covariates and records the HC3 t statistic of 'mutant'. The p-value compares the
observed HC3 t with that null, so the permutation tests the same estimand that is reported, and the covariate design is kept.
Specification pre-declared as primary: any coding mutation, tumours above the 90th percentile of mutation burden excluded.
The other three specifications are sensitivity analyses."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os
from lib_tumour import CHR, coding, prepare, trunc
import numpy as np, pandas as pd, statsmodels.api as sm
def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
B = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
rows = []
for c, gm in [x for x in [('BLCA', ['STAG2']), ('UCEC', ['CTCF'])] if x[0] in sys.argv[1].split(',')]:
    D, chrs, samp, cov = prepare(c, cn='continuous'); e = lagscore(D, chrs, (1, 2, 3)) - lagscore(D, chrs, (20, 30))
    np.save(OUTDIR + f'cisexcess_{c}.npy', e); cov.assign(cis_excess=e).to_csv(OUTDIR + f'tumour_scores_{c}.csv')
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, hq); y = e[keep]
            m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(gm)].Sample_ID)).astype(float)[keep]
            Xs = cov[keep]; Xs = (Xs - Xs.mean()) / Xs.std()
            Xfull = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xs.index).join(Xs)).values; Xred = sm.add_constant(Xs.values)
            full = sm.OLS(y, Xfull).fit(cov_type='HC3'); t_obs = full.tvalues[1]; base = y[m == 0].mean()
            red = sm.OLS(y, Xred).fit(); fit_r, res_r = red.fittedvalues, red.resid
            rng = np.random.default_rng(1); tnull = np.empty(B)
            for b in range(B): tnull[b] = sm.OLS(fit_r + rng.permutation(res_r), Xfull).fit(cov_type='HC3').tvalues[1]
            rows.append({'cohort': c, 'gene': gm[0], 'class': mclass, 'tmb_q': hq, 'primary': mclass == 'all coding' and hq == 0.9,
                         'n_mut': int(m.sum()), 'n': len(y), 'adj_pct': 100 * full.params[1] / base,
                         'ci_low': 100 * full.conf_int()[1, 0] / base, 'ci_high': 100 * full.conf_int()[1, 1] / base,
                         'p_HC3': full.pvalues[1], 'p_freedman_lane': (np.sum(np.abs(tnull) >= abs(t_obs)) + 1) / (B + 1),
                         'corr_mutant_with_covariates_max': float(np.max(np.abs([np.corrcoef(m, Xs[k])[0, 1] for k in Xs.columns])))})
            print(rows[-1], flush=True)
out = OUTDIR + 'cohesin_freedman_lane.csv'
pd.DataFrame(rows).to_csv(out, mode='a', header=not os.path.exists(out), index=False)
POSLAYERS_FILE_024_END
echo "  new      scripts/cohesin_freedman_lane.py"
mkdir -p "scripts"
cat > "scripts/cohesin_v2.py" << 'POSLAYERS_FILE_025_END'
"""Per-tumour cis-excess score and covariate-adjusted test of cohesin/CTCF loss.
Score: for each tumour, mean product of pooled-standardised residual expression of genes 1-3 positions apart minus that of
genes 20-30 apart (GC- and own-copy-number-corrected expression). Covariates: expression-based immune and stromal scores
(purity proxy), copy-number burden, log mutation burden, cohort-specific subtype score."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm
from lib_tumour import coding, prepare, scores, trunc
rows = []
for c, genes_mut in [('BLCA', ['STAG2']), ('UCEC', ['CTCF']), ('UCEC', ['STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']), ('STAD', ['CTCF', 'STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']), ('COAD', ['CTCF', 'STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL'])]:
    D, chrs, samp, cov = prepare(c); e = scores(D, chrs)
    for mclass, src in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            b = cov.log_tmb.values; keep = b <= np.quantile(b, hq)
            mut = pd.Index(samp).isin(set(src[src.gene.isin(genes_mut)].Sample_ID)).astype(float)
            y = e[keep]; m = mut[keep]; Xc = cov[keep]
            if m.sum() < 10: continue
            Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join((Xc - Xc.mean()) / Xc.std()))
            fit = sm.OLS(y, Xd).fit(cov_type='HC3'); base = y[m == 0].mean()
            # permutation of the mutant label on covariate-residualised scores
            r = y - sm.OLS(y, sm.add_constant(((Xc - Xc.mean()) / Xc.std()).values)).fit().fittedvalues
            obs = r[m == 1].mean() - r[m == 0].mean(); rng = np.random.default_rng(0)
            null = np.array([(lambda p: r[p == 1].mean() - r[p == 0].mean())(rng.permutation(m)) for _ in range(10000)])
            rows.append({'cohort': c, 'genes': '/'.join(genes_mut) if len(genes_mut) < 3 else 'cohesin+CTCF' if 'CTCF' in genes_mut else 'cohesin', 'mutation_class': mclass, 'tmb_quantile_kept': hq,
                         'n_mutant': int(m.sum()), 'n_wildtype': int((m == 0).sum()), 'mean_cis_excess_wt': base,
                         'raw_diff_pct': 100 * (y[m == 1].mean() - base) / base, 'adj_effect': fit.params['mutant'], 'adj_effect_pct': 100 * fit.params['mutant'] / base,
                         'ci_low_pct': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high_pct': 100 * fit.conf_int().loc['mutant', 1] / base,
                         'p_adj_HC3': fit.pvalues['mutant'], 'p_perm_10000': (np.sum(np.abs(null) >= abs(obs)) + 1) / 10001})
    print(c, '/'.join(genes_mut)[:20], 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'cohesin_v2_results.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
POSLAYERS_FILE_025_END
echo "  updated  scripts/cohesin_v2.py"
mkdir -p "scripts"
cat > "scripts/coloc.py" << 'POSLAYERS_FILE_026_END'
"""Approximate colocalisation of adjacent genes with SuSiE credible sets: P(same causal variant) = max over credible-set
pairs of sum_v pip1(v)*pip2(v); direction from the shared variant with the highest pip product (afc sign)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
norm = lambda s: s.lower().replace('-', '_')
PAIRS = pd.concat([pd.read_csv(f, usecols=['g1', 'g2']) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')]).drop_duplicates(); genes = set(PAIRS.g1) | set(PAIRS.g2)
atlas = {norm(os.path.basename(f)[6:-7]) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
for f in sys.argv[1:]:
    t = norm(os.path.basename(f).replace('_v11_eQTLs_SuSiE_summary.parquet', '')); out = OUTDIR + f'coloc_{t}.csv.gz'
    if t not in atlas or os.path.exists(out): continue
    E = pq.read_table(f, columns=['phenotype_id', 'variant_id', 'pip', 'cs_id', 'afc']).to_pandas(); E['g'] = E.phenotype_id.str.split('.').str[0]; E = E[E.g.isin(genes)]
    CS = {g: [(dict(zip(c.variant_id, c.pip)), dict(zip(c.variant_id, np.sign(c.afc)))) for _, c in d.groupby('cs_id')] for g, d in E.groupby('g')}
    rows = []
    for g1, g2 in zip(PAIRS.g1.values, PAIRS.g2.values):
        a, b = CS.get(g1), CS.get(g2)
        if a is None or b is None: rows.append((g1, g2, a is not None, b is not None, 0.0, np.nan)); continue
        best, sgn = 0.0, np.nan
        for p1, s1 in a:
            for p2, s2 in b:
                sh = p1.keys() & p2.keys()
                if not sh: continue
                pr = sum(p1[v] * p2[v] for v in sh)
                if pr > best:
                    best = pr; vs = sorted(sh, key=lambda v: -p1[v] * p2[v]); sgn = np.nan
                    for v in vs:
                        if np.isfinite(s1.get(v, np.nan)) and np.isfinite(s2.get(v, np.nan)): sgn = s1[v] * s2[v]; break
        rows.append((g1, g2, True, True, best, sgn))
    S = pd.DataFrame(rows, columns=['g1', 'g2', 'fm1', 'fm2', 'p_coloc', 'direction']); S.to_csv(out, index=False)
    print(t, f'fine-mapped pairs {int((S.fm1 & S.fm2).sum())} | coloc>=0.5 {int((S.p_coloc >= 0.5).sum())}', flush=True)
POSLAYERS_FILE_026_END
echo "  updated  scripts/coloc.py"
mkdir -p "scripts"
cat > "scripts/coloc_test.py" << 'POSLAYERS_FILE_027_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import glob, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
norm = lambda s: s.lower().replace('-', '_'); BINS = [-np.inf, 0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
Cc = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'coloc_*.csv.gz')}
rows, dose = [], []
for t in sorted(set(P) & set(S) & set(Cc)):
    d = P[t].merge(S[t], on=['g1', 'g2']).merge(Cc[t], on=['g1', 'g2']); d = d[(d.dist > 0) & d.fm1 & d.fm2].copy()
    d['bin'] = pd.cut(d.dist, BINS).astype(str)
    d['coloc_same'] = ((d.p_coloc >= 0.5) & (d.direction == 1)).astype(int); d['coloc_opp'] = ((d.p_coloc >= 0.5) & (d.direction == -1)).astype(int)
    d['sig_same'] = d.share_same.astype(int); d['sig_opp'] = (d.share_any & ~d.share_same).astype(int)
    m1 = smf.ols('r ~ sig_same + sig_opp + C(bin) + C(orientation)', data=d).fit()
    m2 = smf.ols('r ~ coloc_same + coloc_opp + sig_same + sig_opp + C(bin) + C(orientation)', data=d).fit()
    rows.append({'tissue': t, 'fine_mapped_pairs': len(d), 'n_coloc_same': int(d.coloc_same.sum()), 'n_coloc_opp': int(d.coloc_opp.sum()),
                 'r_coloc_same': d.r[d.coloc_same == 1].mean(), 'r_sig_shared_not_coloc': d.r[(d.sig_same == 1) & (d.p_coloc < 0.1)].mean(), 'r_no_sharing': d.r[(d.sig_same == 0) & (d.sig_opp == 0)].mean(),
                 'b_sig_same_alone': m1.params['sig_same'], 'b_coloc_same': m2.params['coloc_same'], 'p_coloc_same': m2.pvalues['coloc_same'],
                 'b_coloc_opp': m2.params['coloc_opp'], 'b_sig_same_given_coloc': m2.params['sig_same']})
    sm = d[(d.direction == 1) | (d.p_coloc == 0)]
    for lo, hi, lab in [(-1, 1e-9, '0'), (1e-9, 0.1, '0-0.1'), (0.1, 0.5, '0.1-0.5'), (0.5, 0.8, '0.5-0.8'), (0.8, 1.01, '>0.8')]:
        x = sm[(sm.p_coloc > lo) & (sm.p_coloc <= hi)]; dose.append({'tissue': t, 'p_coloc_bin': lab, 'n': len(x), 'mean_r': x.r.mean()})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'coloc_cis_test.csv', index=False)
pd.set_option('display.width', 250)
print(R[['tissue', 'fine_mapped_pairs', 'n_coloc_same', 'n_coloc_opp', 'r_coloc_same', 'r_sig_shared_not_coloc', 'r_no_sharing', 'b_sig_same_alone', 'b_coloc_same', 'b_coloc_opp', 'b_sig_same_given_coloc']].round(3).to_string(index=False))
print('\nmedians:'); print(R.drop(columns='tissue').median().round(4).to_string())
print(f"\ncoloc_same > 0 in {(R.b_coloc_same > 0).sum()}/{len(R)} tissues (P<0.05 in {(R.p_coloc_same < 0.05).sum()}); coloc_opp < 0 in {(R.b_coloc_opp < 0).sum()}/{len(R)}")
print(f"significant-variant sharing coefficient: alone {R.b_sig_same_alone.median():.3f} -> with colocalisation in the model {R.b_sig_same_given_coloc.median():.3f}")
D = pd.DataFrame(dose); D.to_csv(OUTDIR + 'coloc_dose_response.csv', index=False)
print('\ndose-response (same-direction or none), median over tissues of mean r:'); print(D.groupby('p_coloc_bin', sort=False).agg(n=('n', 'sum'), median_r=('mean_r', 'median')).round(3).to_string())
POSLAYERS_FILE_027_END
echo "  updated  scripts/coloc_test.py"
mkdir -p "scripts"
cat > "scripts/cont_tests.py" << 'POSLAYERS_FILE_028_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm
from lib_tumour import CHR, coding, prepare, trunc
import sys, os
WANT = sys.argv[1].split(',')
def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
an, te = [], []
for c, gm in [x for x in [('BLCA', ['STAG2']), ('UCEC', ['CTCF']), ('STAD', None), ('COAD', None), ('GBM', None)] if x[0] in WANT]:
    D, D_RAW, chrs, samp, cov = prepare(c, cn='continuous', return_raw=True); cna = cov.cna_burden.values
    for nm, M in [('raw', D_RAW), ('continuous CN corrected', D)]:
        far = lagscore(M, chrs, (20, 30)); near = lagscore(M, chrs, (1, 2, 3))
        an.append({'cohort': c, 'expression': nm, 'mean_far': far.mean(), 'rho_far_cna': pd.Series(far).corr(pd.Series(cna), method='spearman'), 'rho_cisexcess_cna': pd.Series(near - far).corr(pd.Series(cna), method='spearman')})
    if gm is None: continue
    e = lagscore(D, chrs, (1, 2, 3)) - lagscore(D, chrs, (20, 30))
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, hq); y = e[keep]; m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(gm)].Sample_ID)).astype(float)[keep]
            Xc = cov[keep]; Xs = (Xc - Xc.mean()) / Xc.std(); Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join(Xs))
            fit = sm.OLS(y, Xd).fit(cov_type='HC3'); base = y[m == 0].mean(); r = y - sm.OLS(y, sm.add_constant(Xs.values)).fit().fittedvalues
            obs = r[m == 1].mean() - r[m == 0].mean(); rng = np.random.default_rng(0); null = np.array([(lambda p: r[p == 1].mean() - r[p == 0].mean())(rng.permutation(m)) for _ in range(10000)])
            te.append({'cohort': c, 'gene': gm[0], 'class': mclass, 'tmb_q': hq, 'n_mut': int(m.sum()), 'adj_pct': 100 * fit.params['mutant'] / base,
                       'ci_low': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high': 100 * fit.conf_int().loc['mutant', 1] / base, 'p_perm': (np.sum(np.abs(null) >= abs(obs)) + 1) / 10001})
    print(c, 'done', flush=True)
A = pd.DataFrame(an); T = pd.DataFrame(te)
A.to_csv(OUTDIR + 'aneuploidy_continuousCN.csv', mode='a', header=not os.path.exists(OUTDIR + 'aneuploidy_continuousCN.csv'), index=False)
if len(T): T.to_csv(OUTDIR + 'cohesin_continuousCN.csv', mode='a', header=not os.path.exists(OUTDIR + 'cohesin_continuousCN.csv'), index=False)
pd.set_option('display.width', 200); print(A.round(3).to_string(index=False)); print(T.round(3).to_string(index=False))
POSLAYERS_FILE_028_END
echo "  updated  scripts/cont_tests.py"
mkdir -p "scripts"
cat > "scripts/crispri_analyse.py" << 'POSLAYERS_FILE_029_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, statsmodels.formula.api as smf
D = pd.read_csv(OUTDIR + 'crispri_pairs.csv.gz'); D = D[np.isfinite(D.resp) & np.isfinite(D.r)]; D['resp'] = D.resp.clip(-20, 20); D['rk'] = D.r * D.kd
CI = D[D.type == 'cis'].copy(); TR = D[D.type == 'trans'].copy()
CI['dbin'] = pd.cut(CI.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['a<10kb', 'b10-50kb', 'c50-200kb', 'd200-500kb', 'e0.5-1Mb']).astype(str)
print(f'cis pairs {len(CI):,} | trans pairs {len(TR):,} | perturbations {D.pert.nunique():,}')
gc = lambda d: {'groups': d.pert.astype('category').cat.codes}
mc = smf.ols('resp ~ rk + r + kd + C(dbin) + C(dbin):kd', data=CI).fit(cov_type='cluster', cov_kwds=gc(CI))
mt = smf.ols('resp ~ rk + r + kd', data=TR).fit(cov_type='cluster', cov_kwds=gc(TR))
for nm, m in (('cis', mc), ('trans', mt)):
    ci = m.conf_int().loc['rk']; print(f"{nm:5s}: coupling x knockdown {m.params['rk']:.3f} (95% CI {ci[0]:.3f} to {ci[1]:.3f}), P = {m.pvalues['rk']:.1e}")
# by distance: does coupling matter beyond the KRAB spreading distance?
for b, d in CI.groupby('dbin'):
    m = smf.ols('resp ~ rk + r + kd', data=d).fit(cov_type='cluster', cov_kwds=gc(d)); print(f"  {b[1:]:>10}: n={len(d):,} slope {m.params['rk']:.3f} P={m.pvalues['rk']:.1e} | mean response {d.resp.mean():.3f}")
S = CI[CI.kd > 2].copy(); qs = S.r.quantile([1 / 3, 2 / 3]).values
S['rt'] = pd.cut(S.r, [-1, qs[0], qs[1], 1], labels=['low r', 'mid r', 'high r']); TS = TR[TR.kd > 2].copy(); TS['rt'] = pd.cut(TS.r, [-1, qs[0], qs[1], 1], labels=['low r', 'mid r', 'high r'])
print('\nstrong knockdowns (>4-fold): mean response by coupling tertile')
tab = S.pivot_table(index='dbin', columns='rt', values='resp', aggfunc='mean', observed=True); tab.loc['trans (other chromosomes)'] = TS.groupby('rt', observed=True).resp.mean(); print(tab.round(3).to_string())
POSLAYERS_FILE_029_END
echo "  updated  scripts/crispri_analyse.py"
mkdir -p "scripts"
cat > "scripts/crispri_test.py" << 'POSLAYERS_FILE_030_END'
"""Does CRISPRi silencing of a gene propagate to its cis-coupled neighbours? (Replogle 2022 K562 genome-wide Perturb-seq)
Response: z-normalised pseudobulk (relative to non-targeting controls), robustly rescaled per perturbation.
Baseline coupling: correlation across 615 BeatAML2 leukaemias (GC-corrected), for target-gene pairs.
Cis: measured genes within 1 Mb of the target. Trans control: genes on other chromosomes with the same coupling range."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, anndata as ad, pyannotables as pa, statsmodels.formula.api as smf
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
G['chr'] = G.Chromosome.astype(str); G['tss'] = np.where(G.Strand.astype(str) == '+', G.Start, G.End)
A = ad.read_h5ad(DATA + 'perturb/K562_gwps_normalized_bulk_01.h5ad'); X = np.asarray(A.X, dtype=np.float32)
obs = A.obs.copy(); obs['target'] = obs.index.str.split('_').str[-1]; genes = np.array(A.var.index)
keep = (obs.fold_expr < 0.5) & obs.target.isin(G.index) & ~obs.index.str.contains('non-targeting') & (obs.num_cells_filtered >= 25)
obs = obs[keep]; X = X[keep.values]; print('effective perturbations (>=50% knockdown):', len(obs))
med = np.median(X, 1, keepdims=True); mad = np.median(np.abs(X - med), 1, keepdims=True) * 1.4826 + 1e-6; R = (X - med) / mad    # robust per-perturbation scaling
# baseline coupling from BeatAML2
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
B = pd.read_csv(DATA + 'beataml/beataml_waves1to4_counts_dbgap.txt', sep='\t', low_memory=False); B['gid'] = B.stable_id.str.split('.').str[0]
B = B.drop_duplicates('gid').set_index('gid')[[x for x in B.columns if x.startswith('BA')]]; B = B.loc[B.index.intersection(BM.index)]; B = B[B.median(axis=1) >= 5]
Yb = np.log2(B.values / B.values.sum(0) * 1e6 + 1); Db = Yb - Yb.mean(1, keepdims=True); gcz = BM.loc[B.index, 'gc'].values; gcz = (gcz - gcz.mean()) / gcz.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); Db = Db - G1 @ np.linalg.lstsq(G1, Db, rcond=None)[0]; Zb = (Db - Db.mean(1, keepdims=True)) / (Db.std(1, keepdims=True) + 1e-9)
bpos = {g: i for i, g in enumerate(B.index)}; vpos = {g: i for i, g in enumerate(genes)}
meas = [g for g in genes if g in G.index and g in bpos]; mg = G.loc[meas]
rng = np.random.default_rng(0); rows = []
for k, (pid, o) in enumerate(obs.iterrows()):
    t = o.target
    if t not in bpos: continue
    tc, tt = G.at[t, 'chr'], G.at[t, 'tss']; kd = -np.log2(max(o.fold_expr, 0.01))
    same = mg[(mg.chr == tc) & ((mg.tss - tt).abs() <= 1e6) & (mg.index != t)]
    other = mg[mg.chr != tc]; other = other.iloc[rng.choice(len(other), min(60, len(other)), replace=False)]
    for typ, sub in (('cis', same), ('trans', other)):
        if not len(sub): continue
        r = (Zb[[bpos[g] for g in sub.index]] * Zb[bpos[t]]).mean(1); resp = R[k, [vpos[g] for g in sub.index]]
        d = (sub.tss - tt).abs().values if typ == 'cis' else np.full(len(sub), np.nan)
        rows.append(pd.DataFrame({'pert': pid, 'target': t, 'type': typ, 'kd': kd, 'r': r, 'resp': resp, 'dist': d}))
D = pd.concat(rows, ignore_index=True); D.to_csv(OUTDIR + 'crispri_pairs.csv.gz', index=False)
D['rk'] = D.r * D.kd; C = D[D.type == 'cis'].copy(); T = D[D.type == 'trans']
C['dbin'] = pd.cut(C.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['<10kb', '10-50kb', '50-200kb', '200-500kb', '0.5-1Mb'])
print(f'cis pairs {len(C):,} | trans pairs {len(T):,}')
mc = smf.ols('resp ~ rk + r + kd + C(dbin) + C(dbin):kd', data=C).fit(cov_type='cluster', cov_kwds={'groups': C.pert.astype('category').cat.codes})
mt = smf.ols('resp ~ rk + r + kd', data=T).fit(cov_type='cluster', cov_kwds={'groups': T.pert.astype('category').cat.codes})
print(f"cis:   coupling x knockdown {mc.params['rk']:.3f} (95% CI {mc.conf_int().loc['rk', 0]:.3f} to {mc.conf_int().loc['rk', 1]:.3f}), P = {mc.pvalues['rk']:.1e}")
print(f"trans: coupling x knockdown {mt.params['rk']:.3f} (95% CI {mt.conf_int().loc['rk', 0]:.3f} to {mt.conf_int().loc['rk', 1]:.3f}), P = {mt.pvalues['rk']:.1e}")
print('\nmean response of cis neighbours by distance (all) and by coupling tertile, strong knockdowns (kd>2):')
S = C[C.kd > 2].copy(); S['rt'] = pd.qcut(S.r, 3, labels=['low r', 'mid r', 'high r'])
print(S.pivot_table(index='dbin', columns='rt', values='resp', aggfunc='mean', observed=True).round(3).to_string())
St = T[T.kd > 2].copy(); St['rt'] = pd.cut(St.r, [-1, S.r.quantile(1 / 3), S.r.quantile(2 / 3), 1], labels=['low r', 'mid r', 'high r'])
print('trans genes, same coupling cut-offs:', St.groupby('rt', observed=True).resp.mean().round(3).to_dict())
POSLAYERS_FILE_030_END
echo "  updated  scripts/crispri_test.py"
mkdir -p "scripts"
cat > "scripts/derive.py" << 'POSLAYERS_FILE_031_END'
"""Test of the hub model: r = r0 + A*(1 + (c*/c)^(2/3))^(-3/2) versus a pure power law r = r0 + A*c^k.
Pairs binned by Hi-C contact (all distances 25 kb-2 Mb); local slope of the coupling-contact relation across contact."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
hub = lambda c, r0, A, cs: r0 + A * (1 + (cs / c) ** (2 / 3)) ** (-1.5)
pw = lambda c, r0, A, k: r0 + A * c ** k
out = []
for cell in ['GM12878', 'IMR90']:
    H = pd.read_csv(OUTDIR + f'hic_coupling_{cell}.csv.gz'); H = H[H.contact_KR > 0]
    H['cb'] = pd.qcut(np.log10(H.contact_KR), 20, labels=False)
    b = H.groupby('cb').agg(c=('contact_KR', 'median'), r=('r', 'mean'), se=('r', lambda x: x.std() / np.sqrt(len(x))), d=('tss_distance', 'median'), n=('r', 'size'))
    w = 1 / b.se ** 2
    p1, _ = curve_fit(hub, b.c, b.r, p0=[0.005, 0.2, 500], sigma=b.se, maxfev=50000, bounds=([-0.05, 0, 1], [0.05, 5, 1e6]))
    p2, _ = curve_fit(pw, b.c, b.r, p0=[0.0, 0.002, 0.5], sigma=b.se, maxfev=50000, bounds=([-0.05, 0, 0.01], [0.05, 5, 3]))
    chi1 = np.sum(w * (b.r - hub(b.c, *p1)) ** 2); chi2 = np.sum(w * (b.r - pw(b.c, *p2)) ** 2)
    x = (p1[2] / b.c) ** (2 / 3); b['k_model'] = x / (1 + x)
    # empirical local slope of (r - r0) vs contact, from adjacent contact bins
    y = np.log(np.clip(b.r - p1[0], 1e-5, None)); lc = np.log(b.c); b['k_empirical'] = np.gradient(y, lc)
    print(f'== {cell}: hub model r0={p1[0]:.4f} A={p1[1]:.3f} c*={p1[2]:.0f}  chi2={chi1:.1f} | power law k={p2[2]:.2f} chi2={chi2:.1f}  (20 bins)')
    dstar = np.interp(np.log(p1[2]), np.log(b.c.values[::-1]), b.d.values[::-1])
    print(f'   distance at which contact = c* (separation = hub size): ~{dstar / 1e3:.0f} kb')
    print(b[['n', 'd', 'c', 'r', 'k_model', 'k_empirical']].round(4).to_string())
    b['cell'] = cell; out.append(b)
pd.concat(out).to_csv(OUTDIR + 'derivation_hub_model.csv')
POSLAYERS_FILE_031_END
echo "  updated  scripts/derive.py"
mkdir -p "scripts"
cat > "scripts/drug_test.py" << 'POSLAYERS_FILE_032_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys
THR = float(sys.argv[1]) if len(sys.argv) > 1 else 10
import glob, os, re, numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
G['chr'] = G.Chromosome.astype(str); G['len'] = G.End - G.Start
# ---- SLAM-seq: newly synthesised reads per gene (Entrez), mapped to Ensembl by 3' UTR overlap on the same strand
def load(f):
    d = pd.read_csv(f, sep='\t'); return d.groupby('Name').agg(chr=('Chromosome', 'first'), s=('Start', 'min'), e=('End', 'max'), strand=('Strand', 'first'), tc=('TcReadCount', 'sum'), tot=('ReadCount', 'sum'))
files = {re.search(r'Sample(\d+)', f).group(1): f for f in glob.glob(DATA + 'slam/GSM*_tcount.tsv.gz')}
ref = load(files['26']); ref['chr'] = ref.chr.str.replace('chr', '')
emap = {}
for c, d in ref.groupby('chr'):
    g = G[G.chr == c]
    for ez, r in d.iterrows():
        hit = g[(g.Start <= r.e) & (g.End >= r.s) & (g.Strand == r.strand)]
        if len(hit) == 1: emap[ez] = hit.index[0]
TC = pd.DataFrame({k: load(f).tc for k, f in files.items()}).fillna(0); TC = TC[TC.index.isin(emap)]; TC.index = [emap[i] for i in TC.index]; TC = TC.groupby(level=0).sum()
print('genes mapped:', len(TC))
# ---- baseline coupling in myeloid leukaemia (BeatAML2), adjacent genes in GRCh38 order, GC-corrected
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
B = pd.read_csv(DATA + 'beataml/beataml_waves1to4_counts_dbgap.txt', sep='\t'); B['gid'] = B.stable_id.str.split('.').str[0]
B = B.drop_duplicates('gid').set_index('gid')[[x for x in B.columns if x.startswith('BA')]]; B = B.loc[B.index.intersection(G.index).intersection(BM.index)]; B = B[B.median(axis=1) >= 10]
g = G.loc[B.index].sort_values(['chr', 'Start']); B = B.loc[g.index]
Y = np.log2(B.values / B.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[g.index, 'gc'].values; gcz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]; Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
same = g.chr.values[:-1] == g.chr.values[1:]; i = np.where(same)[0]
P = pd.DataFrame({'g1': g.index[i], 'g2': g.index[i + 1], 'r': (Z[i] * Z[i + 1]).mean(1), 'dist': g.Start.values[i + 1] - g.End.values[i]}); P = P[P.dist > 0]
print('baseline pairs:', len(P), '| median r %.3f' % P.r.median())
# ---- drug responses
COMP = [('MOLM-13', 'JQ1 (exp 1)', 'chromatin/transcription', ['26', '27', '28'], ['29', '30', '31']), ('MOLM-13', 'JQ1 (exp 2)', 'chromatin/transcription', ['35', '36', '37'], ['38', '39', '40']),
        ('MV4-11', 'JQ1', 'chromatin/transcription', ['20', '21', '22'], ['23', '24', '25']), ('K562', 'JQ1', 'chromatin/transcription', ['14', '15', '16'], ['17', '18', '19']),
        ('K562', 'BRD4 degradation', 'chromatin/transcription', ['8', '9', '10'], ['11', '12', '13']), ('MOLM-13', 'NVP-2 6 nM (CDK9i)', 'chromatin/transcription', ['35', '36', '37'], ['41', '42', '43']),
        ('MOLM-13', 'NVP-2 60 nM (CDK9i)', 'chromatin/transcription', ['35', '36', '37'], ['44', '45']), ('K562', 'flavopiridol (CDK9i)', 'chromatin/transcription', ['1'], ['2']),
        ('K562', 'nilotinib (BCR-ABL)', 'signalling', ['3'], ['7']), ('K562', 'trametinib (MEK)', 'signalling', ['3'], ['5']), ('K562', 'MK-2206 (AKT)', 'signalling', ['3'], ['4']), ('K562', 'MK-2206 + trametinib', 'signalling', ['3'], ['6'])]
lc = np.log10(G.reindex(TC.index)['len'].values.astype(float)); rows = []
for cell, drug, cls, ctrl, trt in COMP:
    X = np.log2(TC[ctrl + trt] / TC[ctrl + trt].sum() * 1e6 + 1); keep = TC[ctrl].mean(1) >= THR
    a, b = X[trt].values, X[ctrl].values; lfc = a.mean(1) - b.mean(1)
    stat = lfc / np.sqrt(a.var(1, ddof=1) / a.shape[1] + b.var(1, ddof=1) / b.shape[1] + 0.01) if len(trt) > 1 and len(ctrl) > 1 else lfc
    s = pd.Series(stat, index=TC.index)[keep]; base = b.mean(1)[keep.values]
    Xc = np.column_stack([np.ones(keep.sum()), lc[keep.values], base, base ** 2]); ok = np.isfinite(Xc).all(1)
    res = pd.Series(np.nan, index=s.index); res[ok] = s[ok].values - Xc[ok] @ np.linalg.lstsq(Xc[ok], s[ok].values, rcond=None)[0]   # remove length/expression-level trends
    d = P.assign(t1=res.reindex(P.g1).values, t2=res.reindex(P.g2).values).dropna()
    z1 = (d.t1 - d.t1.mean()) / d.t1.std(); z2 = (d.t2 - d.t2.mean()) / d.t2.std(); prod = z1 * z2
    slope, _, _, p, _ = stats.linregress(d.r, prod); q = pd.qcut(d.r, 5, labels=False)
    rows.append({'cell': cell, 'perturbation': drug, 'class': cls, 'replicates': f'{len(ctrl)}v{len(trt)}', 'pairs': len(d), 'slope_on_baseline_coupling': slope, 'p': p,
                 'concordance_low_coupling': np.corrcoef(d.t1[q == 0], d.t2[q == 0])[0, 1], 'concordance_high_coupling': np.corrcoef(d.t1[q == 4], d.t2[q == 4])[0, 1]})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + f'drug_cis_results_thr{int(THR)}.csv', index=False); pd.set_option('display.width', 220); print(R.round(3).to_string(index=False))
print('\nmedian slope by class:'); print(R.groupby('class').slope_on_baseline_coupling.median().round(3).to_string())
print('Mann-Whitney chromatin vs signalling P = %.3f' % stats.mannwhitneyu(R[R['class'] == 'chromatin/transcription'].slope_on_baseline_coupling, R[R['class'] == 'signalling'].slope_on_baseline_coupling).pvalue)
POSLAYERS_FILE_032_END
echo "  updated  scripts/drug_test.py"
mkdir -p "scripts"
cat > "scripts/edit_cis.py" << 'POSLAYERS_FILE_033_END'
"""Do therapeutic CRISPR edits perturb the cis neighbours of the edited locus, in proportion to their baseline coupling?
Baseline coupling: GTEx whole blood (GC- and technically corrected), computed earlier. Edits: BCL11A erythroid enhancer
(chr2) and HBG1/2 promoters (chr11, beta-globin locus). Predictions are made from genomic position alone, before looking
at the edited-cell data."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
C = pd.read_csv(DATA + 'edit/GSE264491_merged_counts.csv.gz', index_col=0)
grp = ['unedited'] * 3 + ['BCL11A_enhancer'] * 3 + ['HBG1_2_promoter'] * 3
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
ann = G.loc[G.index.intersection(C.index), ['Chromosome', 'Start', 'End', 'gene_name']]; ann.columns = ['chr', 'start', 'end', 'sym']; ann['chr'] = ann.chr.astype(str)
C = C.loc[ann.index]; C = C[C.median(axis=1) >= 10]; ann = ann.loc[C.index]
cpm = np.log2(C / C.sum(0) * 1e6 + 1)
def lfc(g): return cpm[[c for c, x in zip(C.columns, grp) if x == g]].mean(1) - cpm[[c for c, x in zip(C.columns, grp) if x == 'unedited']].mean(1)
L = pd.DataFrame({'BCL11A_enhancer': lfc('BCL11A_enhancer'), 'HBG1_2_promoter': lfc('HBG1_2_promoter')})
sd = cpm[[c for c, x in zip(C.columns, grp) if x == 'unedited']].std(1)
# baseline coupling from GTEx whole blood
P = pd.read_csv(OUTDIR + 'pairs_whole_blood.csv.gz'); nb = {}
for a, b, r in zip(P.g1, P.g2, P.r): nb.setdefault(a, []).append((b, r)); nb.setdefault(b, []).append((a, r))
def neighbours(target_sym, k=12):
    tid = ann.index[ann.sym == target_sym]
    if not len(tid): return pd.DataFrame()
    t = tid[0]; c = ann.loc[t, 'chr']; sub = ann[ann.chr == c].sort_values('start'); pos = list(sub.index).index(t)
    out = []
    for j in range(max(0, pos - k), min(len(sub), pos + k + 1)):
        if j == pos: continue
        gid = sub.index[j]; r = dict(nb.get(t, [])).get(gid, np.nan)
        out.append({'gene': sub.loc[gid, 'sym'], 'gid': gid, 'rank_from_target': j - pos, 'dist_kb': (sub.loc[gid, 'start'] - ann.loc[t, 'end']) / 1e3,
                    'baseline_coupling_r': r, 'lfc_BCL11A_edit': L.loc[gid, 'BCL11A_enhancer'], 'lfc_HBG_edit': L.loc[gid, 'HBG1_2_promoter']})
    return pd.DataFrame(out)
print('=== Neighbourhood of BCL11A (chr2), edited target of the approved therapy ===')
NB = neighbours('BCL11A'); print(NB.round(3).to_string(index=False))
print('\n=== Neighbourhood of HBG1 (chr11, beta-globin locus), target of reni-cel ===')
NH = neighbours('HBG1'); print(NH.round(3).to_string(index=False))
# local-window enrichment vs matched background: |lfc| of genes within +/-k positions of the edited gene
def window_test(sym, col, k=10, nperm=20000):
    tid = ann.index[ann.sym == sym][0]; c = ann.loc[tid, 'chr']; sub = ann[ann.chr == c].sort_values('start'); pos = list(sub.index).index(tid)
    idx = [sub.index[j] for j in range(max(0, pos - k), min(len(sub), pos + k + 1)) if j != pos]
    obs = L.loc[idx, col].abs().mean(); rng = np.random.default_rng(0); allg = L.index
    null = np.array([L.loc[rng.choice(allg, len(idx), replace=False), col].abs().mean() for _ in range(nperm)])
    return obs, null.mean(), (np.sum(null >= obs) + 1) / (nperm + 1), len(idx)
rows = []
for sym in ['BCL11A', 'HBG1']:
    for col in ['BCL11A_enhancer', 'HBG1_2_promoter']:
        o, n, p, k = window_test(sym, col); rows.append({'edited_locus_window': sym, 'edit': col, 'n_neighbours': k, 'mean_|lfc|_window': o, 'mean_|lfc|_random': n, 'P': p})
W = pd.DataFrame(rows); print('\n=== Local perturbation around each edited locus ===\n' + W.round(4).to_string(index=False))
# genome-wide: does baseline coupling to a strongly changed gene predict a neighbour's change?
res = []
for col in ['BCL11A_enhancer', 'HBG1_2_promoter']:
    d = P.assign(l1=L[col].reindex(P.g1).values, l2=L[col].reindex(P.g2).values, s1=sd.reindex(P.g1).values, s2=sd.reindex(P.g2).values).dropna()
    d = d[d.dist > 0]; big = d[(d.l1.abs() > 0.5) | (d.l2.abs() > 0.5)].copy()
    drv = np.where(big.l1.abs() >= big.l2.abs(), big.l1, big.l2); nbv = np.where(big.l1.abs() >= big.l2.abs(), big.l2, big.l1)
    res.append({'edit': col, 'n_pairs_with_driver': len(big), 'pearson_pred_vs_obs': stats.pearsonr(big.r * drv, nbv)[0], 'slope': np.polyfit(big.r * drv, nbv, 1)[0],
                'corr_lfc_neighbours_all_pairs': stats.pearsonr(d.l1, d.l2)[0], 'corr_shuffled': stats.pearsonr(d.l1, np.random.default_rng(1).permutation(d.l2.values))[0]})
R = pd.DataFrame(res); print('\n=== Bystander model genome-wide ===\n' + R.round(3).to_string(index=False))
NB.to_csv(OUTDIR + 'edit_BCL11A_neighbours.csv', index=False); NH.to_csv(OUTDIR + 'edit_HBG_neighbours.csv', index=False)
W.to_csv(OUTDIR + 'edit_window_tests.csv', index=False); R.to_csv(OUTDIR + 'edit_bystander_model.csv', index=False)
POSLAYERS_FILE_033_END
echo "  updated  scripts/edit_cis.py"
mkdir -p "scripts"
cat > "scripts/eqtl_law.py" << 'POSLAYERS_FILE_034_END'
"""Review point 10: eQTL law on a common scale, with every shared variant's own effect.
Observed coupling: correlation of the two genes in GTEx's normalized (inverse-normal) expression, the scale on which GTEx slopes are
estimated, after removing 15 expression PCs as a stand-in for the hidden factors in the eQTL model (raw INT also kept).
Predicted genetic correlation: sum over pairs of credible sets, and over each shared variant v, of
    PIP_1(v) * PIP_2(v) * 2 p_v (1 - p_v) * beta_1(v) * beta_2(v)
using each variant's own GTEx slope for each gene. The overlap max over credible-set pairs of sum_v PIP_1 PIP_2 is reported as a
colocalisation SCORE, not a posterior probability of a shared causal variant."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
U = DATA + 'gtex_eqtl/'; OUT = OUTDIR + 'eqtl_law_v2_pairs.csv'
done = set(pd.read_csv(OUT).tissue) if os.path.exists(OUT) else set()
NAME = {'nerve_tibial': 'Nerve_Tibial', 'thyroid': 'Thyroid', 'cells_cultured_fibroblasts': 'Cells_Cultured_fibroblasts', 'artery_tibial': 'Artery_Tibial',
        'whole_blood': 'Whole_Blood', 'testis': 'Testis', 'skin_sun_exposed_lower_leg': 'Skin_Sun_Exposed_Lower_leg', 'esophagus_mucosa': 'Esophagus_Mucosa',
        'adipose_subcutaneous': 'Adipose_Subcutaneous', 'lung': 'Lung'}
for t, N in NAME.items():
    if t in done: continue
    P = pd.read_csv(OUTDIR + f'pred_pairs_{t}.csv.gz'); P = P[(P.dist > 0) & (P.p_coloc >= 0.1)]; genes = set(P.g1) | set(P.g2)
    # observed coupling on GTEx's normalized expression
    bed = pd.read_csv(U + f'{N}_v11_normalized_expression_bed.gz', sep='\t'); bed['gid'] = bed.iloc[:, 3].str.split('.').str[0]
    X = bed.drop_duplicates('gid').set_index('gid').iloc[:, 4:].astype(float)
    Xc = X.values - X.values.mean(1, keepdims=True)
    U_, s_, Vt = np.linalg.svd(Xc, full_matrices=False); Xp = Xc - (U_[:, :15] * s_[:15]) @ Vt[:15]
    idx = {g: i for i, g in enumerate(X.index)}
    def corr(M, a, b):
        x, y = M[a] - M[a].mean(), M[b] - M[b].mean(); return float(x @ y / np.sqrt((x @ x) * (y @ y)))
    # credible sets and slopes
    E = pq.read_table(U + f'{N}_v11_eQTLs_SuSiE_summary.parquet', columns=['phenotype_id', 'variant_id', 'pip', 'cs_id']).to_pandas()
    E['gid'] = E.phenotype_id.str.split('.').str[0]; E = E[E.gid.isin(genes)]
    CS = {g: [dict(zip(c.variant_id, c.pip)) for _, c in d.groupby('cs_id')] for g, d in E.groupby('gid')}
    need = set(E.variant_id); pf = pq.ParquetFile(U + f'{N}_v11_eQTLs_signif_pairs.parquet'); slope, af = {}, {}
    for rg in range(pf.num_row_groups):
        S = pf.read_row_group(rg, columns=['phenotype_id', 'variant_id', 'slope', 'af']).to_pandas(); S = S[S.variant_id.isin(need)]
        S['gid'] = S.phenotype_id.str.split('.').str[0]; S = S[S.gid.isin(genes)]
        slope.update({(a, b): s for a, b, s in zip(S.gid, S.variant_id, S.slope)}); af.update(dict(zip(S.variant_id, S.af)))
    rows = []
    for g1, g2 in zip(P.g1, P.g2):
        if g1 not in idx or g2 not in idx: continue
        score, pred, lead_pred, cov_w, tot_w = 0.0, 0.0, 0.0, 0.0, 0.0
        for p1 in CS.get(g1, []):
            for p2 in CS.get(g2, []):
                sh = p1.keys() & p2.keys()
                if not sh: continue
                score = max(score, sum(p1[v] * p2[v] for v in sh))
                lv = max(sh, key=lambda v: p1[v] * p2[v]); ws = sum(p1[v] * p2[v] for v in sh)
                for v in sh:
                    w = p1[v] * p2[v]; tot_w += w
                    if (g1, v) in slope and (g2, v) in slope:
                        p = af[v]; pred += w * 2 * p * (1 - p) * slope[(g1, v)] * slope[(g2, v)]; cov_w += w
                if (g1, lv) in slope and (g2, lv) in slope:
                    p = af[lv]; lead_pred += ws * 2 * p * (1 - p) * slope[(g1, lv)] * slope[(g2, lv)]
        if tot_w == 0: continue
        rows.append({'tissue': t, 'g1': g1, 'g2': g2, 'coloc_score': score, 'pred_r': pred, 'pred_r_lead_only': lead_pred,
                     'weight_with_slopes': cov_w / tot_w, 'obs_r_int_pc15': corr(Xp, idx[g1], idx[g2]), 'obs_r_int_raw': corr(Xc, idx[g1], idx[g2])})
    R = pd.DataFrame(rows); R.to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
    print(t, len(R), 'pairs | slope coverage %.2f | r(pred, obs_pc15) %.2f' % (R.weight_with_slopes.mean(), R[['pred_r', 'obs_r_int_pc15']].corr().iloc[0, 1]), flush=True)
    del bed, X, Xc, Xp, U_, Vt
POSLAYERS_FILE_034_END
echo "  new      scripts/eqtl_law.py"
mkdir -p "scripts"
cat > "scripts/eqtl_law_summary.py" << 'POSLAYERS_FILE_035_END'
"""eQTL law across tissues (ST4A): Pearson r and the slope of observed on predicted coupling with tissue intercepts, on both observed
scales and for all shared variants versus the lead variant only, with 95% intervals from 1,000 resamples of 10-Mb genomic blocks."""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR
d = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
d['block'] = d.tissue + ':' + G.Chromosome.reindex(d.g1).astype(str).values + ':' + (G.Start.reindex(d.g1).fillna(0).values // 1e7).astype(int).astype(str)
def slope(df, y, x): X = np.column_stack([df[x].values, pd.get_dummies(df.tissue).values.astype(float)]); return np.linalg.lstsq(X, df[y].values, rcond=None)[0][0]
blocks = d.block.unique(); grp = {b: np.where(d.block.values == b)[0] for b in blocks}; rng = np.random.default_rng(7); res = []
for y in ('obs_r_int_pc15', 'obs_r_int_raw'):
    for x in ('pred_r', 'pred_r_lead_only'):
        bs, rs = [], []
        for _ in range(1000):
            s = d.iloc[np.concatenate([grp[b] for b in rng.choice(blocks, len(blocks))])]; bs.append(slope(s, y, x)); rs.append(s[[x, y]].corr().iloc[0, 1])
        res.append({'observed_scale': y, 'prediction': x, 'pairs': len(d), 'pearson_r': d[[x, y]].corr().iloc[0, 1], 'r_ci': '%.2f-%.2f' % tuple(np.percentile(rs, [2.5, 97.5])),
                    'slope': slope(d, y, x), 'slope_ci': '%.2f-%.2f' % tuple(np.percentile(bs, [2.5, 97.5]))})
R = pd.DataFrame(res); R.to_csv(OUTDIR + 'eqtl_law_v2_summary.csv', index=False); print(R.round(3).to_string(index=False))
POSLAYERS_FILE_035_END
echo "  new      scripts/eqtl_law_summary.py"
mkdir -p "scripts"
cat > "scripts/eqtl_share.py" << 'POSLAYERS_FILE_036_END'
"""Per tissue: for every adjacent gene pair of the atlas, are both genes eGenes, and do they share a significant eQTL variant
(any direction; same direction of effect)? Reads only the needed parquet columns."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
PAIRS = pd.concat([pd.read_csv(f, usecols=['g1', 'g2']) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')]).drop_duplicates()
genes = set(PAIRS.g1) | set(PAIRS.g2)
norm = lambda s: s.lower().replace('-', '_')
for f in sys.argv[1:]:
    t = norm(os.path.basename(f).replace('_v11_eQTLs_signif_pairs.parquet', '')); out = OUTDIR + f'shares_{t}.csv.gz'
    if os.path.exists(out): continue
    pf = pq.ParquetFile(f); sign = {}; nE = 0
    for rg in range(pf.num_row_groups):          # stream row groups to bound memory
        E = pf.read_row_group(rg, columns=['phenotype_id', 'variant_id', 'slope']).to_pandas()
        E['g'] = E.phenotype_id.str.split('.').str[0]; E = E[E.g.isin(genes)]; nE += len(E)
        for g, d in E.groupby('g'):
            sign.setdefault(g, {}).update(zip(d.variant_id.values, np.sign(d.slope.values).astype(np.int8)))
        del E
    E = range(nE)
    rows = []
    for g1, g2 in zip(PAIRS.g1.values, PAIRS.g2.values):
        a, b = sign.get(g1), sign.get(g2)
        if a is None or b is None: rows.append((g1, g2, a is not None, b is not None, False, False, 0)); continue
        sh = a.keys() & b.keys(); same = sum(1 for v in sh if a[v] == b[v])
        rows.append((g1, g2, True, True, len(sh) > 0, same > 0, len(sh)))
    S = pd.DataFrame(rows, columns=['g1', 'g2', 'egene1', 'egene2', 'share_any', 'share_same', 'n_shared'])
    S.to_csv(out, index=False); print(t, f'{len(E):,} eQTL pairs | both eGenes {int((S.egene1 & S.egene2).sum())} | share {int(S.share_any.sum())}', flush=True)
POSLAYERS_FILE_036_END
echo "  updated  scripts/eqtl_share.py"
mkdir -p "scripts"
cat > "scripts/eqtl_test.py" << 'POSLAYERS_FILE_037_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import glob, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
norm = lambda s: s.lower().replace('-', '_')
BINS = [-np.inf, 0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
tis = sorted(set(P) & set(S)); print('tissues with coupling and eQTL:', len(tis))
rows = []
for t in tis:
    d = P[t].merge(S[t], on=['g1', 'g2']); d = d[d.dist > 0]                       # non-overlapping pairs only (avoids read sharing)
    d['bin'] = pd.cut(d.dist, BINS).astype(str); be = d[d.egene1 & d.egene2].copy()
    be['same'] = be.share_same.astype(int); be['opp_only'] = (be.share_any & ~be.share_same).astype(int)
    m = smf.ols('r ~ same + opp_only + C(bin) + C(orientation)', data=be).fit()
    grp = lambda q: d.loc[q, 'r'].mean()
    r_all = d.r.mean(); delta = m.params['same']; f_same = be.same.sum() / len(d)
    rows.append({'tissue': t, 'pairs': len(d), 'both_eGenes': len(be), 'frac_share_same_of_all_pairs': f_same,
                 'r_share_same': grp(d.egene1 & d.egene2 & d.share_same), 'r_eGenes_no_share': grp(d.egene1 & d.egene2 & ~d.share_any), 'r_not_both_eGenes': grp(~(d.egene1 & d.egene2)),
                 'delta_same_adj': delta, 'p_same': m.pvalues['same'], 'delta_opposite_only_adj': m.params['opp_only'],
                 'share_of_cis_from_shared_eQTL': f_same * delta / r_all})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'eqtl_cis_test.csv', index=False)
pd.set_option('display.width', 250); print(R.round(3).to_string(index=False))
print('\nmedian over tissues:'); print(R.drop(columns='tissue').median().round(4).to_string())
# tissue specificity: own-tissue sharing vs other-tissue sharing, on pairs that are eGene pairs in both tissues
spec = []
for t in tis:
    d0 = P[t][P[t].dist > 0][['g1', 'g2', 'r', 'dist', 'orientation']]; d0['bin'] = pd.cut(d0.dist, BINS).astype(str)
    own = S[t][['g1', 'g2', 'egene1', 'egene2', 'share_same']].rename(columns={'share_same': 'own', 'egene1': 'oe1', 'egene2': 'oe2'})
    for u in tis:
        if u == t: continue
        oth = S[u][['g1', 'g2', 'egene1', 'egene2', 'share_same']].rename(columns={'share_same': 'other'})
        x = d0.merge(own, on=['g1', 'g2']).merge(oth, on=['g1', 'g2']); x = x[x.oe1 & x.oe2 & x.egene1 & x.egene2]
        if len(x) < 500: continue
        x['own'] = x.own.astype(int); x['other'] = x.other.astype(int)
        m = smf.ols('r ~ own + other + C(bin) + C(orientation)', data=x).fit()
        spec.append({'coupling_tissue': t, 'eqtl_tissue': u, 'n': len(x), 'b_own': m.params['own'], 'b_other': m.params['other']})
SP = pd.DataFrame(spec); SP.to_csv(OUTDIR + 'eqtl_tissue_specificity.csv', index=False)
print(f'\ntissue specificity (joint models, {len(SP)} tissue pairs): own-tissue sharing coefficient median {SP.b_own.median():.3f}; other-tissue sharing coefficient median {SP.b_other.median():.3f}; own > other in {100*(SP.b_own > SP.b_other).mean():.0f}% of tissue pairs')
POSLAYERS_FILE_037_END
echo "  updated  scripts/eqtl_test.py"
mkdir -p "scripts"
cat > "scripts/extract.py" << 'POSLAYERS_FILE_038_END'
"""Stream the GDC-PANCAN HTSeq matrix (log2(count+1)) from the zip and keep protein-coding genes for selected cohorts."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, subprocess, numpy as np, pandas as pd
cohorts = sys.argv[1:]
BP = pd.read_csv(DATA + 'tcga/GDC-PANCAN_basic_phenotype.tsv', sep='\t')
pcol = [c for c in BP.columns if 'project' in c.lower()][0]; BP['proj'] = BP[pcol].astype(str).str.replace('TCGA-', '')
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False)
pc = set(BM.loc[BM['Gene type'] == 'protein_coding', 'Gene stable ID'])
p = subprocess.Popen(['unzip', '-p', DATA + 'tcga/GDC-PANCAN.htseq_counts.tsv.zip', 'GDC-PANCAN.htseq_counts.tsv'], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
hdr = p.stdout.readline().rstrip('\n').split('\t')
samp2proj = dict(zip(BP['sample'], BP.proj))
cols = {c: [i for i, s in enumerate(hdr) if i > 0 and samp2proj.get(s) == c and s[13:15] in ('01', '03', '11')] for c in cohorts}
for c in cohorts: print(c, len(cols[c]), 'samples', flush=True)
allidx = np.array(sorted(set(i for v in cols.values() for i in v))); pos = {i: k for k, i in enumerate(allidx)}
genes, rows = [], []
for n, line in enumerate(p.stdout):
    g = line[:line.index('\t')].split('.')[0]
    if g not in pc: continue
    v = np.array(line.rstrip('\n').split('\t'), dtype=object)[allidx].astype(np.float32); genes.append(g); rows.append(v)
X = np.vstack(rows); print('genes kept', len(genes), flush=True)
for c in cohorts:
    sel = [pos[i] for i in cols[c]]
    np.savez_compressed(OUTDIR + f'expr_{c}.npz', X=X[:, sel], genes=np.array(genes), samples=np.array([hdr[i] for i in cols[c]]))
print('done')
POSLAYERS_FILE_038_END
echo "  updated  scripts/extract.py"
mkdir -p "scripts"
cat > "scripts/fric1.py" << 'POSLAYERS_FILE_039_END'
"""Friction 1: is the TAD effect explained by 3D contact? All gene pairs within 2 Mb (GM12878 / IMR-90 Hi-C, GTEx EBV / fibroblasts)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from lib_tad import tads, assign
rows = []
for cell, tadname in [('GM12878', 'GM12878_lymphoblastoid_Lieberman'), ('IMR90', 'IMR90_fetalLungFibroblast_Lieberman')]:
    H = pd.read_csv(OUTDIR + f'hic_coupling_{cell}.csv.gz'); a = assign(tads(tadname))
    H['t1'] = a.reindex(H.g1).values; H['t2'] = a.reindex(H.g2).values; H = H[(H.t1 >= 0) & (H.t2 >= 0)].copy(); H['same_tad'] = (H.t1 == H.t2).astype(int)
    H['db'] = pd.cut(H.tss_distance, [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6]).astype(str)
    for b, d in H.groupby('db'):
        if d.same_tad.sum() < 50 or (1 - d.same_tad).sum() < 50: continue
        m0 = smf.ols('r ~ same_tad + log_d', data=d).fit(); m1 = smf.ols('r ~ same_tad + log_d + log_oe', data=d).fit()
        rows.append({'cell': cell, 'distance': b, 'pairs': len(d), 'frac_same_tad': d.same_tad.mean(), 'r_same': d.r[d.same_tad == 1].mean(), 'r_diff': d.r[d.same_tad == 0].mean(),
                     'TAD_effect': m0.params['same_tad'], 'p_TAD': m0.pvalues['same_tad'], 'TAD_effect_given_contact': m1.params['same_tad'], 'p_TAD_given_contact': m1.pvalues['same_tad'],
                     'contact_effect_given_TAD': m1.params['log_oe'], 'p_contact_given_TAD': m1.pvalues['log_oe']})
    m0 = smf.ols('r ~ same_tad + bs(log_d, df=5)', data=H).fit(); m1 = smf.ols('r ~ same_tad + bs(log_d, df=5) + log_oe', data=H).fit()
    print(f"{cell} all distances: TAD effect {m0.params['same_tad']:.4f} (t {m0.tvalues['same_tad']:.1f}) -> given contact {m1.params['same_tad']:.4f} (t {m1.tvalues['same_tad']:.1f}); contact given TAD t {m1.tvalues['log_oe']:.1f}")
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'friction1_tad_vs_contact.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
# short distances (< 20 kb, unresolved by 25-kb Hi-C): adjacent pairs in GTEx, same vs different TAD
P = []
for tissue, tadname in [('cells_ebv-transformed_lymphocytes', 'GM12878_lymphoblastoid_Lieberman'), ('cells_cultured_fibroblasts', 'IMR90_fetalLungFibroblast_Lieberman')]:
    d = pd.read_csv(OUTDIR + f'pairs_{tissue}.csv.gz'); a = assign(tads(tadname)); d['t1'] = a.reindex(d.g1).values; d['t2'] = a.reindex(d.g2).values
    d = d[(d.t1 >= 0) & (d.t2 >= 0) & (d.dist > 0)]; d['same'] = d.t1 == d.t2; d['db'] = pd.cut(d.dist, [0, 5e3, 2e4, 5e4, 2e5]).astype(str)
    for b, x in d.groupby('db'): P.append({'tissue': tissue, 'distance': b, 'n_same': int(x.same.sum()), 'n_diff': int((~x.same).sum()), 'r_same': x.r[x.same].mean(), 'r_diff': x.r[~x.same].mean()})
print(pd.DataFrame(P).round(3).to_string(index=False)); pd.DataFrame(P).to_csv(OUTDIR + 'friction1_short_adjacent.csv', index=False)
POSLAYERS_FILE_039_END
echo "  updated  scripts/fric1.py"
mkdir -p "scripts"
cat > "scripts/fric2.py" << 'POSLAYERS_FILE_040_END'
"""Friction 2: is the elevated clustering of active genes within TADs in tumours (Liang et al. 2026) a copy-number effect?
Clustering z per sample: among the top-3,000 expressed protein-coding genes, number of TADs with >= 2 active genes versus
200 random gene sets of the same size (as in their bootstrap test). TADs: segments between conserved boundaries
(stability percentile >= 0.5, McArthur & Capra), genes placed with GRCh37 coordinates."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import glob, numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
B = pd.read_csv(glob.glob(DATA + 'TAD-full/*/data/boundariesByStability/100kbBookendBoundaries_mainText/100kbBookendBoundaries_byStability.bed')[0], sep='\t')
B = B[B.stability_percentile >= 0.5]; B['chr'] = B.chr.str.replace('chr', ''); B['mid'] = (B['loc'] + B['loc2']) / 2
G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()]; G37 = G37[G37.Chromosome.astype(str).isin([str(i) for i in range(1, 23)] + ['X'])]
mid = (G37.Start + G37.End) / 2; tad = pd.Series(-1, index=G37.index); off = 0
for c, bb in B.groupby('chr'):
    idx = G37.index[G37.Chromosome.astype(str) == c]; k = np.searchsorted(np.sort(bb['mid'].values), mid[idx].values); tad[idx] = off + k; off += len(bb) + 1
rng = np.random.default_rng(0); K = 3000
def cluster_z(Y, ok_mask=None):
    """Y genes x samples (log CPM), tad index per gene in TI; returns z per sample"""
    z = []
    for s in range(Y.shape[1]):
        cand = np.where(ok_mask[:, s])[0] if ok_mask is not None else np.arange(Y.shape[0])
        if len(cand) < 1.5 * K: z.append(np.nan); continue
        top = cand[np.argsort(-Y[cand, s])[:K]]; obs = np.sum(np.bincount(TI[top]) >= 2)
        null = [np.sum(np.bincount(TI[rng.choice(cand, K, replace=False)]) >= 2) for _ in range(200)]
        z.append((obs - np.mean(null)) / (np.std(null) + 1e-9))
    return np.array(z)
rows = []
for c in ['BLCA', 'UCEC', 'STAD', 'COAD']:
    z_ = np.load(OUTDIR + f'expr_{c}.npz', allow_pickle=True); X = z_['X']; genes = z_['genes']; samp = z_['samples']
    keep = pd.Index(genes).isin(tad.index[tad >= 0]); X = X[keep]; genes = genes[keep]; TI = tad.reindex(genes).values.astype(int)
    Y = np.log2(np.clip(2 ** X - 1, 0, None) / np.clip(2 ** X - 1, 0, None).sum(0) * 1e6 + 1)
    typ = np.array([s[13:15] for s in samp]); T = typ == '01'; N = typ == '11'
    CNz = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CN = pd.DataFrame(CNz['CN'], index=CNz['genes'], columns=CNz['samples'])
    tum = np.where(T & pd.Index(samp).isin(CN.columns))[0]; nor = np.where(N)[0]
    CNt = CN.reindex(index=genes, columns=samp[tum]).fillna(0).values; burden = np.mean(np.abs(CNt) > 0.2, 0)
    zN = cluster_z(Y[:, nor]); zT = cluster_z(Y[:, tum]); M = np.abs(CNt) <= 0.2; zT_clean = cluster_z(Y[:, tum], ok_mask=M)
    zN_masked = cluster_z(Y[:, nor], ok_mask=M[:, rng.choice(M.shape[1], len(nor), replace=True)])   # same gene exclusion as random tumours
    rows.append({'cohort': c, 'normals': len(nor), 'tumours': len(tum), 'z_normal': np.median(zN), 'z_tumour': np.median(zT), 'P_tumour_vs_normal': stats.mannwhitneyu(zT, zN).pvalue,
                 'rho_z_vs_CNA_burden': stats.spearmanr(zT, burden)[0], 'P_rho': stats.spearmanr(zT, burden)[1],
                 'z_tumour_CNA_genes_removed': np.nanmedian(zT_clean), 'z_normal_same_exclusion': np.nanmedian(zN_masked), 'P_clean_vs_normal_matched': stats.mannwhitneyu(zT_clean[np.isfinite(zT_clean)], zN_masked[np.isfinite(zN_masked)]).pvalue})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'friction2_tumour_tad_clustering.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
POSLAYERS_FILE_040_END
echo "  updated  scripts/fric2.py"
mkdir -p "scripts"
cat > "scripts/gc_correlated_sim.py" << 'POSLAYERS_FILE_041_END'
"""Review point 12: GC correction when part of the real cis signal tracks GC (Supplementary Table 9C).
A per-sample regulatory factor acts on each gene in proportion to its locally smoothed GC content, making up 0-75% of the cis variance;
we report the fraction of the TRUE cis lag correlation that survives GC correction."""
import numpy as np, pandas as pd
from poslayers import simulate_genome
from poslayers.decompose import gc_correct, lag_profile_linear
from poslayers.config import OUTDIR
rows = []
for share in (0.0, 0.25, 0.5, 0.75):
    out = []
    for rep in range(3):
        s = simulate_genome(n_genes=1000, n_chrom=4, n_samples=200, decay_len=5, gc_bias_sd=0.25, cis_sd=0.6, seed=50 + rep)
        rng = np.random.default_rng(90 + rep); gc = s['gc']; cis = s['cis']
        bio = np.outer(rng.normal(0, 1, cis.shape[0]), np.convolve(gc, np.ones(9) / 9, 'same')); bio *= cis.std() / bio.std()
        new = np.sqrt(1 - share) * cis + np.sqrt(share) * bio
        truth = lag_profile_linear(new - new.mean(0), 1000, 4, 10); kept = lag_profile_linear(gc_correct(new - new.mean(0), gc), 1000, 4, 10)
        out.append((truth[1], truth[10], kept[1] / truth[1], kept[10] / truth[10]))
    rows.append((share, *np.mean(out, 0)))
R = pd.DataFrame(rows, columns=['gc_tracking_share', 'true_lag1', 'true_lag10', 'surviving_lag1', 'surviving_lag10'])
R.to_csv(OUTDIR + 'gc_correlated_cis_sim.csv', index=False); print(R.round(3).to_string(index=False))
POSLAYERS_FILE_041_END
echo "  new      scripts/gc_correlated_sim.py"
mkdir -p "scripts"
cat > "scripts/hic_extract.py" << 'POSLAYERS_FILE_042_END'
"""Extracts Hi-C contacts between gene promoters (TSS) up to 2 Mb apart, streaming only the needed regions
from the public Rao et al. 2014 in situ Hi-C maps (hg19). Nothing large is downloaded.

Usage:
    python scripts/make_tss_table.py     # writes DATA/hic/genes_hg19_tss.tsv
    python scripts/hic_extract.py
Output: DATA/hic/hic_contacts_<cell>.csv.gz and DATA/hic/hic_expected_<cell>.csv.gz (a few tens of MB each).
"""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import gzip, sys, numpy as np, pandas as pd, hicstraw
os.makedirs(DATA + 'hic', exist_ok=True)
MAPS = {'GM12878': 'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',
        'IMR90': 'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic'}
RES, MAXD, WIN = 25000, 2_000_000, 6_000_000          # bin size, maximum distance, window streamed at a time
G = pd.read_csv(DATA + 'hic/genes_hg19_tss.tsv', sep='\t', dtype={'chr': str})
for name, url in MAPS.items():
    print(f'== {name}', flush=True); hic = hicstraw.HiCFile(url)
    lengths = {c.name: c.length for c in hic.getChromosomes()}
    rows, expected = [], []
    for c, gc in G.groupby('chr'):
        if c not in lengths: continue
        try: mzd = hic.getMatrixZoomData(c, c, 'observed', 'KR', 'BP', RES)
        except Exception as e: print('  skip', c, e); continue
        band = {}
        for start in range(0, lengths[c], WIN - MAXD):
            end = min(start + WIN, lengths[c])
            for r in mzd.getRecords(start, end, start, end):
                d = r.binY - r.binX
                if 0 <= d <= MAXD and np.isfinite(r.counts): band[(r.binX, r.binY)] = r.counts
            if end >= lengths[c]: break
        # expected contact per distance (zeros included): sum of counts / number of bin pairs at that distance
        nb = lengths[c] // RES + 1; s = {}
        for (x, y), v in band.items(): s[y - x] = s.get(y - x, 0.0) + v
        for d in range(0, MAXD // RES + 1): expected.append((c, d * RES, s.get(d * RES, 0.0) / max(nb - d, 1)))
        gc = gc.sort_values('tss'); b = (gc.tss.values // RES) * RES; ids = gc.gene_id.values
        for i in range(len(ids)):                                  # each gene with every gene up to 2 Mb downstream
            j = i + 1
            while j < len(ids) and gc.tss.values[j] - gc.tss.values[i] <= MAXD:
                x, y = min(b[i], b[j]), max(b[i], b[j])
                rows.append((ids[i], ids[j], c, int(gc.tss.values[j] - gc.tss.values[i]), band.get((x, y), 0.0), int(y - x)))
                j += 1
        print(f'  chr{c}: {len(band):,} contacts in band, {len(rows):,} gene pairs so far', flush=True)
    pd.DataFrame(rows, columns=['g1', 'g2', 'chr', 'tss_distance', 'contact_KR', 'bin_distance']).to_csv(DATA + f'hic/hic_contacts_{name}.csv.gz', index=False)
    pd.DataFrame(expected, columns=['chr', 'bin_distance', 'expected_KR']).to_csv(DATA + f'hic/hic_expected_{name}.csv.gz', index=False)
    print(f'  written hic_contacts_{name}.csv.gz', flush=True)
POSLAYERS_FILE_042_END
echo "  new      scripts/hic_extract.py"
mkdir -p "scripts"
cat > "scripts/hic_test.py" << 'POSLAYERS_FILE_043_END'
"""Distance vs 3D contact (GM12878 in situ Hi-C, KR, 25 kb) as predictors of cis coupling in GTEx EBV lymphocytes."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.formula.api as smf
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
C = pd.read_csv(DATA + 'gtex/gene_reads_adult_gtex_v11_cells_ebv-transformed_lymphocytes_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description')
C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]; samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[C.index, 'gc'].values; gcz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T
Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9); pos = {g: i for i, g in enumerate(C.index)}
H = pd.read_csv(DATA + 'hic/hic_contacts_GM12878.csv.gz', dtype={'chr': str}); E = pd.read_csv(DATA + 'hic/hic_expected_GM12878.csv.gz', dtype={'chr': str})
H = H[H.g1.isin(pos) & H.g2.isin(pos)].copy()
i1 = H.g1.map(pos).values; i2 = H.g2.map(pos).values; r = np.empty(len(H))
for s in range(0, len(H), 200000): r[s:s + 200000] = (Z[i1[s:s + 200000]] * Z[i2[s:s + 200000]]).mean(1)
H['r'] = r; H = H.merge(E, on=['chr', 'bin_distance'], how='left'); H = H[(H.bin_distance > 0) & (H.expected_KR > 0) & (H.tss_distance > 0)]
H['oe'] = H.contact_KR / H.expected_KR; H['log_oe'] = np.log2(H.oe + 0.05); H['log_d'] = np.log10(H.tss_distance); H['log_c'] = np.log10(H.contact_KR + 0.1)
print(f'gene pairs within 2 Mb with Hi-C and coupling: {len(H):,}')
# decay laws
H['dbin'] = pd.cut(H.tss_distance, [25e3, 50e3, 100e3, 200e3, 500e3, 1e6, 2e6])
dec = H.groupby('dbin', observed=True).agg(n=('r', 'size'), mean_r=('r', 'mean'), mean_contact=('contact_KR', 'mean'), mid=('tss_distance', 'median'))
print(dec.round(4).to_string())
print('log-log slope: coupling %.2f | contact %.2f' % (np.polyfit(np.log10(dec.mid), np.log10(dec.mean_r.clip(lower=1e-4)), 1)[0], np.polyfit(np.log10(dec.mid), np.log10(dec.mean_contact), 1)[0]))
# dissociation: at equal distance, does observed/expected contact predict coupling?
rows = []
for b, d in H.groupby('dbin', observed=True):
    q = pd.qcut(d.oe.rank(method='first'), 5, labels=False)
    rows.append({'distance': str(b), 'pairs': len(d), 'spearman_r_vs_OE': stats.spearmanr(d.oe, d.r)[0], 'r_lowest_OE_quintile': d.r[q == 0].mean(), 'r_highest_OE_quintile': d.r[q == 4].mean()})
print(pd.DataFrame(rows).round(4).to_string(index=False))
m1 = smf.ols('r ~ bs(log_d, df=5)', data=H).fit(); m2 = smf.ols('r ~ bs(log_d, df=5) + log_oe', data=H).fit(); m3 = smf.ols('r ~ bs(log_c, df=5)', data=H).fit(); m4 = smf.ols('r ~ bs(log_c, df=5) + bs(log_d, df=5)', data=H).fit()
print(f"\nR2 distance only {m1.rsquared:.4f} | + O/E contact {m2.rsquared:.4f} (O/E coef {m2.params['log_oe']:.4f}, t {m2.tvalues['log_oe']:.1f})")
print(f"R2 contact only {m3.rsquared:.4f} | contact + distance {m4.rsquared:.4f} | AIC distance {m1.aic:.0f} vs contact {m3.aic:.0f}")
H.drop(columns='dbin').to_csv(OUTDIR + 'hic_coupling_GM12878.csv.gz', index=False)
POSLAYERS_FILE_043_END
echo "  updated  scripts/hic_test.py"
mkdir -p "scripts"
cat > "scripts/hic_test_imr90.py" << 'POSLAYERS_FILE_044_END'
"""Distance vs 3D contact (GM12878 in situ Hi-C, KR, 25 kb) as predictors of cis coupling in GTEx EBV lymphocytes."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.formula.api as smf
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
C = pd.read_csv(DATA + 'gtex/gene_reads_adult_gtex_v11_cells_cultured_fibroblasts_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description')
C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]; samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[C.index, 'gc'].values; gcz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T
Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9); pos = {g: i for i, g in enumerate(C.index)}
H = pd.read_csv(DATA + 'hic/hic_contacts_IMR90.csv.gz', dtype={'chr': str}); E = pd.read_csv(DATA + 'hic/hic_expected_IMR90.csv.gz', dtype={'chr': str})
H = H[H.g1.isin(pos) & H.g2.isin(pos)].copy()
i1 = H.g1.map(pos).values; i2 = H.g2.map(pos).values; r = np.empty(len(H))
for s in range(0, len(H), 200000): r[s:s + 200000] = (Z[i1[s:s + 200000]] * Z[i2[s:s + 200000]]).mean(1)
H['r'] = r; H = H.merge(E, on=['chr', 'bin_distance'], how='left'); H = H[(H.bin_distance > 0) & (H.expected_KR > 0) & (H.tss_distance > 0)]
H['oe'] = H.contact_KR / H.expected_KR; H['log_oe'] = np.log2(H.oe + 0.05); H['log_d'] = np.log10(H.tss_distance); H['log_c'] = np.log10(H.contact_KR + 0.1)
print(f'gene pairs within 2 Mb with Hi-C and coupling: {len(H):,}')
# decay laws
H['dbin'] = pd.cut(H.tss_distance, [25e3, 50e3, 100e3, 200e3, 500e3, 1e6, 2e6])
dec = H.groupby('dbin', observed=True).agg(n=('r', 'size'), mean_r=('r', 'mean'), mean_contact=('contact_KR', 'mean'), mid=('tss_distance', 'median'))
print(dec.round(4).to_string())
print('log-log slope: coupling %.2f | contact %.2f' % (np.polyfit(np.log10(dec.mid), np.log10(dec.mean_r.clip(lower=1e-4)), 1)[0], np.polyfit(np.log10(dec.mid), np.log10(dec.mean_contact), 1)[0]))
# dissociation: at equal distance, does observed/expected contact predict coupling?
rows = []
for b, d in H.groupby('dbin', observed=True):
    q = pd.qcut(d.oe.rank(method='first'), 5, labels=False)
    rows.append({'distance': str(b), 'pairs': len(d), 'spearman_r_vs_OE': stats.spearmanr(d.oe, d.r)[0], 'r_lowest_OE_quintile': d.r[q == 0].mean(), 'r_highest_OE_quintile': d.r[q == 4].mean()})
print(pd.DataFrame(rows).round(4).to_string(index=False))
m1 = smf.ols('r ~ bs(log_d, df=5)', data=H).fit(); m2 = smf.ols('r ~ bs(log_d, df=5) + log_oe', data=H).fit(); m3 = smf.ols('r ~ bs(log_c, df=5)', data=H).fit(); m4 = smf.ols('r ~ bs(log_c, df=5) + bs(log_d, df=5)', data=H).fit()
print(f"\nR2 distance only {m1.rsquared:.4f} | + O/E contact {m2.rsquared:.4f} (O/E coef {m2.params['log_oe']:.4f}, t {m2.tvalues['log_oe']:.1f})")
print(f"R2 contact only {m3.rsquared:.4f} | contact + distance {m4.rsquared:.4f} | AIC distance {m1.aic:.0f} vs contact {m3.aic:.0f}")
H.drop(columns='dbin').to_csv(OUTDIR + 'hic_coupling_IMR90.csv.gz', index=False)
POSLAYERS_FILE_044_END
echo "  updated  scripts/hic_test_imr90.py"
mkdir -p "scripts"
cat > "scripts/isochore_law.py" << 'POSLAYERS_FILE_045_END'
"""Isochore law, re-done after review point 6.

For each tissue and lag L we report three things, all on the SAME normalisation (raw per-gene SD):
  old     : corr(raw) - corr(quadratic-GC-corrected)      [the published metric, kept for comparison]
  insample: [cov(raw) - cov(raw minus per-sample LINEAR GC fit)] / (sd_i sd_j)  versus  Var(b) g_i g_j / (sd_i sd_j)
            -> nearly an identity; the gap is the cross term cov(b g, residual), i.e. how far b is from independent of biology
  heldout : b_s estimated on odd chromosomes only; prediction Var(b_odd) g_i g_j / (sd_i sd_j) for gene pairs on EVEN chromosomes,
            observation = GC covariance on even chromosomes from their own fit. This is the non-trivial test: a single
            per-sample scalar measured elsewhere in the genome must predict the positional covariance in held-out chromosomes.
Also: corr(b_odd, b_even) across samples, and the share of the quadratic correction carried by the linear term.
"""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyannotables as pa
CHR = [str(i) for i in range(1, 23)]
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()][['Chromosome', 'Start']]; G.columns = ['chr', 'start']; G['chr'] = G.chr.astype(str)
G = G[G.chr.isin(CHR)].join(BM[['gc']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN']).set_index('SAMPID')
LAGS = [1, 2, 5, 10, 20, 30]
OUT = OUTDIR + 'isochore_v2.csv'
done = set(pd.read_csv(OUT).tissue) if os.path.exists(OUT) else set()

def lag_pairs(chrs, keep, L):
    """index pairs (i, i+L) within the same chromosome, restricted to genes where keep is True"""
    I = []
    for c in np.unique(chrs[keep]):
        idx = np.where((chrs == c) & keep)[0]
        if len(idx) > L + 5: I.append(np.column_stack([idx[:-L], idx[L:]]))
    return np.vstack(I) if I else np.empty((0, 2), int)

def fit_slopes(D, g, rows):
    """per-sample OLS of D[rows, s] on [1, g[rows]]; returns intercepts a_s and slopes b_s"""
    X = np.column_stack([np.ones(rows.sum()), g[rows]]); B = np.linalg.lstsq(X, D[rows], rcond=None)[0]
    return B[0], B[1]

for f in sorted(glob.glob(DATA + 'gtex/*_gct.gz')):
    t = os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', '').lower().replace('-', '_')
    if t == 'kidney_medulla' or t in done: continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
    C = C.loc[C.index.intersection(G.index)]; C = C[C.median(axis=1) >= 10]
    g = G.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)   # genes x samples
    gz = (g.gc.values - g.gc.mean()) / g.gc.std(); S = D.shape[1]
    sd = D.std(1)                                                                  # common normalisation for every quantity
    allg = np.ones(len(gz), bool); odd = np.isin(chrs, [str(i) for i in range(1, 23, 2)]); even = ~odd
    a_all, b_all = fit_slopes(D, gz, allg); a_o, b_o = fit_slopes(D, gz, odd); a_e, b_e = fit_slopes(D, gz, even)
    R_lin = D - (a_all[None, :] + np.outer(gz, b_all))                             # linear GC removed
    Xq = np.column_stack([np.ones(len(gz)), gz, gz ** 2 - (gz ** 2).mean()]); Bq = np.linalg.lstsq(Xq, D, rcond=None)[0]; R_q = D - Xq @ Bq
    R_e = D - (a_e[None, :] + np.outer(gz, b_e))                                   # even chromosomes, their own fit
    cov = lambda M, P: np.mean((M[P[:, 0]] - M[P[:, 0]].mean(1, keepdims=True)) * (M[P[:, 1]] - M[P[:, 1]].mean(1, keepdims=True)), 1)
    def corr(M, P):
        s = M.std(1) + 1e-12; return cov(M, P) / (s[P[:, 0]] * s[P[:, 1]])
    row = {'tissue': t, 'samples': S, 'genes': len(gz), 'var_b': b_all.var(), 'var_b_odd': b_o.var(), 'var_b_even': b_e.var(),
           'corr_b_odd_even': np.corrcoef(b_o, b_e)[0, 1]}
    for L in LAGS:
        P = lag_pairs(chrs, allg, L); nn = sd[P[:, 0]] * sd[P[:, 1]]
        row[f'old_obs_L{L}'] = np.mean(corr(D, P) - corr(R_q, P))
        row[f'old_pred_L{L}'] = np.mean(b_all.var() * gz[P[:, 0]] * gz[P[:, 1]] / nn)
        row[f'in_obs_L{L}'] = np.mean((cov(D, P) - cov(R_lin, P)) / nn)
        row[f'in_pred_L{L}'] = row[f'old_pred_L{L}']
        row[f'quad_obs_L{L}'] = np.mean((cov(D, P) - cov(R_q, P)) / nn)               # full basis, same normalisation
        Pe = lag_pairs(chrs, even, L); ne = sd[Pe[:, 0]] * sd[Pe[:, 1]]
        row[f'held_pred_L{L}'] = np.mean(b_o.var() * gz[Pe[:, 0]] * gz[Pe[:, 1]] / ne)
        row[f'held_obs_L{L}'] = np.mean((cov(D, Pe) - cov(R_e, Pe)) / ne)
    pd.DataFrame([row]).to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
    print(t, S, f"corr(b_odd,b_even)={row['corr_b_odd_even']:.3f}  L10 old {row['old_obs_L10']:.4f}/{row['old_pred_L10']:.4f}  held {row['held_obs_L10']:.4f}/{row['held_pred_L10']:.4f}", flush=True)
POSLAYERS_FILE_045_END
echo "  new      scripts/isochore_law.py"
mkdir -p "scripts"
cat > "scripts/lib_tad.py" << 'POSLAYERS_FILE_046_END'
"""TAD partitions and gene-to-TAD assignment, shared by tad_test.py, fric1.py and tad_bootstrap.py.
TADs: McArthur & Capra 20-bin TAD landscape (hg19); genes placed by GRCh37 midpoint."""
import glob, numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA

CHR = [str(i) for i in range(1, 23)] + ['X']
_found = glob.glob(DATA + 'TAD-full/*/data/20binsTADlandscape')
TD = _found[0] if _found else DATA + 'TAD-full/MISSING/data/20binsTADlandscape'
MATCH = {'liver': 'Liver_leung2015', 'adrenal_gland': 'adrenal_schmitt2016', 'artery_aorta': 'aorta_leung2015', 'bladder': 'bladder_schmitt2016',
         'brain_cortex': 'cortex_DLPFC_schmitt2016', 'brain_frontal_cortex_ba9': 'cortex_DLPFC_schmitt2016', 'heart_left_ventricle': 'leftVentricle_leung2015',
         'lung': 'lung_schmitt2016', 'pancreas': 'pancreas_schmitt2016', 'muscle_skeletal': 'psoasMuscle_schmitt2016',
         'small_intestine_terminal_ileum': 'smallBowel_schmitt2016', 'spleen': 'spleen_schmitt2016',
         'cells_ebv_transformed_lymphocytes': 'GM12878_lymphoblastoid_Lieberman', 'cells_cultured_fibroblasts': 'IMR90_fetalLungFibroblast_Lieberman'}

G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()][['Chromosome', 'Start', 'End']]
G37.columns = ['chr', 'start', 'end']; G37['chr'] = G37.chr.astype(str); G37 = G37[G37.chr.isin(CHR)].copy(); G37['mid'] = (G37.start + G37.end) / 2


def tads(name):
    """TAD intervals for one cell type or tissue from the 20-bin landscape (bin 6 starts closed by a bin-15 end)."""
    d = f'{TD}/{name}/'
    b6 = pd.read_csv(d + f'bin_6_{name}.bed', sep='\t', header=None, names=['chr', 's', 'e'])
    b15 = pd.read_csv(d + f'bin_15_{name}.bed', sep='\t', header=None, names=['chr', 's', 'e'])
    b6['bin'] = b6.e - b6.s + 1; b6['exp_end'] = b6.s + 10 * b6.bin - 1; ends = {c: set(g.e) for c, g in b15.groupby('chr')}
    b6['end'] = [next((e for e in range(r.exp_end - 3, r.exp_end + 4) if e in ends.get(r.chr, set())), np.nan) for _, r in b6.iterrows()]
    T = b6.dropna(subset=['end']).rename(columns={'s': 'start'}); T['chr'] = T.chr.str.replace('chr', '')
    T = T[T.chr.isin(CHR)].sort_values(['chr', 'start']).reset_index(drop=True); T['tad_id'] = np.arange(len(T))
    return T


def assign(T):
    """Series gene_id -> TAD id (-1 when the gene midpoint is outside every TAD)."""
    out = pd.Series(-1, index=G37.index)
    for c, g in T.groupby('chr'):
        idx = G37.index[G37.chr == c]; mm = G37.loc[idx, 'mid'].values; st_ = g.start.values; en = g.end.values.astype(int); ids = g.tad_id.values
        j = np.searchsorted(st_, mm, side='right') - 1; ok = (j >= 0) & (mm <= en[np.clip(j, 0, len(en) - 1)]); out.loc[idx[ok]] = ids[j[ok]]
    return out
POSLAYERS_FILE_046_END
echo "  new      scripts/lib_tad.py"
mkdir -p "scripts"
cat > "scripts/lib_tumour.py" << 'POSLAYERS_FILE_047_END'
"""Shared tumour preprocessing for the TCGA analyses (cohesin_v2, cont_tests, aneuploidy, replicate, tad_tumour, cohesin_freedman_lane).
Previously each script read cohesin_v2.py as text, edited it and ran it with exec; the variants are now options of prepare().

prepare(c, cn='gistic' | 'continuous', codes=('01',), return_raw=False, return_gids=False)
    -> (D, [D_raw,] chrs, samp, cov[, gids])
    D       GC-corrected deviations with each gene's own copy-number effect removed
    D_raw   the same before the copy-number correction
    cn      'gistic': remove the mean of each GISTIC level per gene; 'continuous': regress each gene on its own segment log2 ratio (linear + quadratic)
"""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA, OUTDIR
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene name': 'sym'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc', 'sym']], how='inner')
mu = pd.read_csv(DATA + 'tcga/GDC-PANCAN.mutect2_snv.tsv', sep='\t', usecols=['Sample_ID', 'gene', 'effect', 'filter']); mu = mu[mu['filter'] == 'PASS']
NONCOD = ['synonymous_variant', 'intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant', 'non_coding_transcript_exon_variant']
coding = mu[~mu.effect.isin(NONCOD)]; trunc = coding[coding.effect.str.contains('stop_gained|frameshift|splice_acceptor|splice_donor', regex=True)]
burden = mu.groupby('Sample_ID').size()
IMM = ['PTPRC', 'CD2', 'CD3E', 'CD3D', 'CD247', 'LCK', 'CD48', 'CD53', 'CD52', 'CORO1A', 'LAPTM5', 'CD37', 'FCER1G', 'TYROBP', 'CD74', 'HLA-DRA', 'CCL5', 'IL2RG', 'CXCR4', 'SELL']
STR = ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A3', 'DCN', 'LUM', 'FBN1', 'FAP', 'PDGFRB', 'THY1', 'SPARC', 'VCAN', 'MMP2', 'CDH11', 'FN1', 'POSTN', 'COL11A1', 'SFRP2', 'ACTA2']
SUB = {'LAML': (['MPO', 'ELANE', 'AZU1'], ['CD14', 'LYZ', 'CSF1R']), 'GBM': (['OLIG2', 'SOX2', 'PDGFRA'], ['CHI3L1', 'CD44', 'MET']), 'BLCA': (['GATA3', 'KRT20', 'UPK1B', 'UPK2', 'FOXA1', 'PPARG', 'FGFR3'], ['KRT5', 'KRT6A', 'KRT14', 'CD44', 'CDH3']),
       'UCEC': (['PGR', 'ESR1', 'PAX8', 'MSX1'], ['CDKN2A', 'L1CAM', 'WT1']), 'STAD': (['CDX2', 'MUC2', 'TFF3'], ['VIM', 'ZEB1']), 'COAD': (['CDX2', 'VIL1', 'LGALS4'], ['VIM', 'ZEB1'])}
def prepare(c, cn='gistic', codes=('01',), return_raw=False, return_gids=False):
    z = np.load(OUTDIR + f'expr_{c}.npz', allow_pickle=True); X = z['X']; genes = z['genes']; samp = z['samples']
    keep = np.array([s[13:15] in codes for s in samp]) & pd.Index(samp).isin(burden.index); X = X[:, keep]; samp = samp[keep]
    C = np.clip(2 ** X - 1, 0, None); ok = pd.Index(genes).isin(G38.index); C = C[ok]; genes = genes[ok]
    good = np.median(C, 1) >= 10; C = C[good]; genes = genes[good]
    g = G38.loc[genes].reset_index(); o = np.lexsort((g.start.values, g.chr.values)); C = C[o]; g = g.iloc[o].reset_index(drop=True)
    Y = np.log2(C / C.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
    if cn == 'continuous':
        Zc = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CNd = pd.DataFrame(Zc['CN'], index=Zc['genes'], columns=Zc['samples'])
        have = pd.Index(samp).isin(CNd.columns); D = D[:, have]; Y = Y[:, have]; samp = samp[have]
        CNv = CNd.reindex(index=g.iloc[:, 0].values, columns=samp).fillna(0).values; D_raw = D.copy()
        for i in range(D.shape[0]):                          # remove each gene's own continuous dosage effect (linear + quadratic)
            x = CNv[i]
            if np.std(x) < 1e-6: continue
            X1 = np.column_stack([np.ones(len(x)), x, x ** 2]); D[i] -= X1 @ np.linalg.lstsq(X1, D[i], rcond=None)[0]
        cna_burden = np.mean(np.abs(CNv) > 0.2, 0)
    else:
        hdr = pd.read_csv(DATA + 'tcga/GDC-PANCAN.gistic.tsv', sep='\t', nrows=0).columns
        CN = pd.read_csv(DATA + 'tcga/GDC-PANCAN.gistic.tsv', sep='\t', usecols=[hdr[0]] + [x for x in samp if x in set(hdr)], index_col=0); CN.index = CN.index.str.split('.').str[0]
        CN = CN.reindex(index=g.iloc[:, 0].values, columns=samp); have = (CN.notna().mean(axis=0) > 0.5).values
        D = D[:, have]; Y = Y[:, have]; samp = samp[have]; CNv = np.nan_to_num(CN.values[:, have].astype(float), nan=0.0); D_raw = D.copy()
        for i in range(D.shape[0]):
            cnv = CNv[i]
            for lv in np.unique(cnv): D[i, cnv == lv] -= D[i, cnv == lv].mean()
        cna_burden = np.mean(CNv != 0, 0)
    sym = g.sym.astype(str).values; sc = lambda gs: np.nanmean([(Y[sym == s][0] - Y[sym == s][0].mean()) / (Y[sym == s][0].std() + 1e-9) for s in gs if (sym == s).any()], 0)
    up, dn = SUB[c]
    cov = pd.DataFrame({'immune': sc(IMM), 'stromal': sc(STR), 'subtype': sc(up) - sc(dn), 'cna_burden': cna_burden, 'log_tmb': np.log10(burden.reindex(samp).values + 1)}, index=samp)
    out = [D] + ([D_raw] if return_raw else []) + [g.chr.values, samp, cov] + ([g.iloc[:, 0].values] if return_gids else [])
    return tuple(out)
def scores(D, chrs, lags_cis=(1, 2, 3), lags_far=(20, 30)):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    def lagmean(L):
        parts = [Z[idx[:-L]] * Z[idx[L:]] for c in CHR for idx in [np.where(chrs == c)[0]] if len(idx) > L + 5]
        return np.concatenate(parts, 0).mean(0)
    return np.mean([lagmean(L) for L in lags_cis], 0) - np.mean([lagmean(L) for L in lags_far], 0)
POSLAYERS_FILE_047_END
echo "  new      scripts/lib_tumour.py"
mkdir -p "scripts"
cat > "scripts/make_tss_table.py" << 'POSLAYERS_FILE_048_END'
"""Gene TSS table (hg19) used by hic_extract.py: one row per Ensembl gene on chromosomes 1-22 and X, strand-aware TSS from Ensembl GRCh37."""
import pandas as pd, pyannotables as pa
import os
from poslayers.config import DATA
os.makedirs(DATA + 'hic', exist_ok=True)
G = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G = G[~G.index.duplicated()]
G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)] + ['X'])]
tss = G.Start.where(G.Strand.astype(str).isin(['1', '+']), G.End)
pd.DataFrame({'gene_id': G.index, 'chr': G.Chromosome.astype(str).values, 'tss': tss.astype(int).values}).to_csv(DATA + 'hic/genes_hg19_tss.tsv', sep='\t', index=False)
POSLAYERS_FILE_048_END
echo "  new      scripts/make_tss_table.py"
mkdir -p "scripts"
cat > "scripts/meta_stag2.py" << 'POSLAYERS_FILE_049_END'
"""Review point 4: inverse-variance combination of the STAG2 effect in bladder cancer (TCGA) and AML (BeatAML2), computed from the
exported estimates rather than typed into the figure. Primary specification in both cohorts: STAG2, any coding mutation.
Percent effects are combined on the percent scale; SE is recovered from the HC3 confidence interval (BLCA) or reported directly (AML)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd
from scipy import stats
b = pd.read_csv(OUTDIR + 'cohesin_freedman_lane.csv'); b = b[(b.cohort == 'BLCA') & b.primary].iloc[0]
a = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv'); a = a[a.group == 'STAG2, any coding'].iloc[0]
est = np.array([b.adj_pct, a.adj_pct]); se = np.array([(b.ci_high - b.ci_low) / (2 * 1.959964), a.se_pct]); w = 1 / se ** 2
m = np.sum(w * est) / w.sum(); s = 1 / np.sqrt(w.sum()); z = m / s; Q = np.sum(w * (est - m) ** 2)
row = {'combination': 'STAG2 BLCA + AML, any coding', 'blca_pct': b.adj_pct, 'blca_se': se[0], 'aml_pct': a.adj_pct, 'aml_se': se[1],
       'meta_pct': m, 'meta_se': s, 'ci_low': m - 1.959964 * s, 'ci_high': m + 1.959964 * s, 'p': 2 * stats.norm.sf(abs(z)),
       'cochran_Q': Q, 'p_heterogeneity': stats.chi2.sf(Q, 1)}
pd.DataFrame([row]).to_csv(OUTDIR + 'stag2_meta.csv', index=False); print({k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()})
POSLAYERS_FILE_049_END
echo "  new      scripts/meta_stag2.py"
mkdir -p "scripts"
cat > "scripts/orient.py" << 'POSLAYERS_FILE_050_END'
"""Cis layer by orientation and intergenic distance of adjacent gene pairs, per GTEx tissue (after GC and technical correction)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, pyannotables as pa
OUT = OUTDIR + 'orientation_by_tissue.csv'; PAIRS = OUTDIR + 'pair_correlations.parquet'
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start', 'End', 'Strand']]
G38.columns = ['chr', 'start', 'end', 'strand']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc', 'type']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
BINS = [-np.inf, 0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]; LAB = ['overlap', '0-1kb', '1-5kb', '5-20kb', '20-100kb', '100-500kb', '>500kb']
done = set(pd.read_csv(OUT).tissue) if os.path.exists(OUT) else set(); allpairs = []
for f in sys.argv[1:]:
    t = os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', '')
    if t in done or t == 'kidney_medulla': continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
    C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
    g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); D = resid_genes(D, np.column_stack([gcz, gcz ** 2]))
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    if tech.shape[1] < D.shape[1] - 10: D = resid_samples(D, tech)
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    same = g.chr.values[:-1] == g.chr.values[1:]; i = np.where(same)[0]
    P = pd.DataFrame({'g1': g.index[i], 'g2': g.index[i + 1], 'r': (Z[i] * Z[i + 1]).mean(1), 'dist': g.start.values[i + 1] - g.end.values[i],
                      's1': g.strand.values[i], 's2': g.strand.values[i + 1], 'pc': (g.type.values[i] == 'protein_coding') & (g.type.values[i + 1] == 'protein_coding')})
    P['orientation'] = np.where(P.s1 == P.s2, 'tandem', np.where((P.s1 == '-') & (P.s2 == '+'), 'divergent', 'convergent'))
    P['dist_bin'] = pd.cut(P.dist, BINS, labels=LAB)
    rng = np.random.default_rng(1); j = rng.permutation(len(Z)); floor = float(np.mean((Z[j[:-1]] * Z[j[1:]]).mean(1)))
    S = P.groupby(['orientation', 'dist_bin'], observed=True).r.agg(['size', 'mean']).reset_index().rename(columns={'size': 'n_pairs', 'mean': 'mean_r'})
    S['tissue'] = t; S['random_pair_floor'] = floor
    S.to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
    P['tissue'] = t; P[['tissue', 'g1', 'g2', 'r', 'dist', 'orientation', 'pc']].to_csv(OUTDIR + f'pairs_{t}.csv.gz', index=False)
    print(t, len(P), 'pairs; floor', round(floor, 4), flush=True)
POSLAYERS_FILE_050_END
echo "  updated  scripts/orient.py"
mkdir -p "scripts"
cat > "scripts/pilot_domain_tests.py" << 'POSLAYERS_FILE_051_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, numpy as np, pandas as pd, pyannotables as pa
sys.path.insert(0, os.path.join(os.environ.get('LIVER_SPECTRA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Liver_Spectra')), 'scripts'))  # liver pipeline: github.com/Danpc11/Liver_Spectra
from common import *
X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl')); comp = pd.read_pickle(inter('comp.pkl'))
s = M.index[M.estadio != 'Control']
G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()][['Chromosome', 'Start', 'End']]
k = keep.sort_values(['chr', 'grid_index']).set_index('gene_id').join(G37, how='left')
k['length'] = (k.End - k.Start).clip(lower=200); k['mid'] = (k.Start + k.End) / 2
# local gene density: retained genes within +/-1 Mb (isochore / GC proxy)
dens = []
for c, g in k.groupby('chr'):
    m = g.mid.values; dens += list(np.array([np.sum(np.abs(m - x) <= 1e6) for x in m]))
k['density'] = np.array(dens, dtype=float)
genes = k.index[k.Start.notna() & k.index.isin(A.index)]; k = k.loc[genes]
Y = A.loc[genes, s].values; Yc = Y - Y.mean(0, keepdims=True); D0 = Yc - Yc.mean(1, keepdims=True)         # per-gene deviations across biopsies
mu = Yc.mean(1)
def resid_samples(D, F):   # regress each gene's profile across samples on sample-level factors F (samples x p)
    F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G):     # within each sample, remove trends with gene properties G (genes x q)
    G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
C = comp.loc[s].values
props = np.column_stack([np.log(k.length), np.log(k.length) ** 2, mu, mu ** 2, k.density, k.density ** 2])
D1 = resid_samples(D0, C)                       # composition removed
D2 = resid_genes(D1, (props - props.mean(0)) / props.std(0))   # + gene-property-linked technical trends
U_, S_, Vt = np.linalg.svd(D2 - D2.mean(1, keepdims=True), full_matrices=False); ve = S_ ** 2 / np.sum(S_ ** 2)
def drop_pcs(D, q): return D - (U_[:, :q] * S_[:q]) @ Vt[:q]
variants = {'raw deviations': D0, '- composition': D1, '- composition - gene properties': D2,
            '... - 5 PCs': drop_pcs(D2, 5), '... - 10 PCs': drop_pcs(D2, 10), '... - 20 PCs': drop_pcs(D2, 20)}
print('variance across biopsies explained by PC1-5: ' + ', '.join(f'{100*v:.1f}%' for v in ve[:5]) + f'; top 20: {100*ve[:20].sum():.1f}%')
chrs = k.chr.values; lags = [1, 2, 3, 5, 10, 15, 20, 30, 45, 60, 100, 150, 200, 300]
def lagprof(D, perm=False):
    out = {L: [] for L in lags}
    for c in CHR:
        idx = np.where(chrs == c)[0]
        if perm: idx = np.random.default_rng(3).permutation(idx)
        Z = D[idx]; Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
        for L in lags:
            if len(Z) > L + 5: out[L] += list((Z[:-L] * Z[L:]).mean(1))
    return np.array([np.mean(out[L]) for L in lags])
rows = []
for name, D in variants.items():
    r = lagprof(D); p = lagprof(D, perm=True); ex = r - p
    m = (np.array(lags) <= 100) & (ex > 0.003); b = np.polyfit(np.array(lags)[m], np.log(ex[m]), 1) if m.sum() >= 3 else [np.nan]
    rows.append({'variant': name, 'cis_L1': ex[0], 'L2': ex[1], 'domain_L10_30': ex[[4, 5, 6, 7]].mean(), 'L60_100': ex[[9, 10]].mean(), 'long_L150_300': ex[[11, 12, 13]].mean(), 'decay_length_genes': -1 / b[0]})
R = pd.DataFrame(rows); print(R.round(4).to_string(index=False)); R.to_csv(OUTDIR + 'p4_domain_scale_tests.csv', index=False)
POSLAYERS_FILE_051_END
echo "  new      scripts/pilot_domain_tests.py"
mkdir -p "scripts"
cat > "scripts/pilot_liver_gc.py" << 'POSLAYERS_FILE_052_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.environ.get('LIVER_SPECTRA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Liver_Spectra')), 'scripts'))  # liver pipeline: github.com/Danpc11/Liver_Spectra
from common import *
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False)
BM = BM.rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene start (bp)': 'start', 'Gene end (bp)': 'end', 'Chromosome/scaffold name': 'chr'}).drop_duplicates('gid').set_index('gid')
X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl')); comp = pd.read_pickle(inter('comp.pkl'))
s = M.index[M.estadio != 'Control']
k = keep.sort_values(['chr', 'grid_index']).set_index('gene_id').join(BM[['gc', 'start', 'end']], how='inner'); k = k[k.index.isin(A.index)]
k['length'] = (k.end - k.start).clip(lower=200); k['mid'] = (k.start + k.end) / 2
dens = []
for c, g in k.groupby('chr', sort=False):
    m = g.mid.values; dens += list(np.array([np.sum(np.abs(m - x) <= 1e6) for x in m]))
k['density'] = np.array(dens, float)
genes = k.index; chrs = k.chr.astype(str).values
print(f'liver biopsies: {len(genes)} genes with GC (of {len(A)})')
Y = A.loc[genes, s].values; Yc = Y - Y.mean(0, keepdims=True); mu = Yc.mean(1); D0 = Yc - mu[:, None]
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
zs = lambda a: (a - a.mean(0)) / a.std(0)
gc = k.gc.values; L = np.log(k.length.values)
D1 = resid_samples(D0, comp.loc[s].values)
lags = [1, 2, 3, 5, 10, 15, 20, 30, 150, 200, 300]
def prof(D, perm=False):
    out = {l: [] for l in lags}
    for c in CHR:
        idx = np.where(chrs == c)[0]
        if perm: idx = np.random.default_rng(3).permutation(idx)
        Z = D[idx]; Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
        for l in lags:
            if len(Z) > l + 5: out[l] += list((Z[:-l] * Z[l:]).mean(1))
    return np.array([np.mean(out[l]) for l in lags])
def summ(name, D):
    ex = prof(D) - prof(D, perm=True); return {'adjustment': name, 'cis_L1': ex[0], 'L2': ex[1], 'domain_L10_30': ex[4:8].mean(), 'long_L150_300': ex[8:].mean()}
rows = [summ('- composition', D1),
        summ('- composition - GC', resid_genes(D1, zs(np.column_stack([gc, gc ** 2])))),
        summ('- composition - length, expression, density', resid_genes(D1, zs(np.column_stack([L, L ** 2, mu, mu ** 2, k.density, k.density ** 2])))),
        summ('- composition - all gene properties incl. GC', resid_genes(D1, zs(np.column_stack([gc, gc ** 2, L, L ** 2, mu, mu ** 2, k.density, k.density ** 2]))))]
R = pd.DataFrame(rows); print(R.round(4).to_string(index=False))
# GC is itself clustered along the genome (isochores): neighbour correlation of gene GC
gz = []; 
for c in CHR:
    v = k.gc.values[chrs == c]; v = (v - v.mean()) / v.std(); gz.append([np.mean(v[:-l] * v[l:]) for l in (1, 10, 30)])
print('autocorrelation of gene GC content along the genome at 1, 10, 30 genes:', np.round(np.mean(gz, 0), 3))
# per-biopsy GC slope: how strongly each sample's deviations track gene GC
slope = np.array([np.polyfit(zs(gc), D1[:, i], 1)[0] for i in range(D1.shape[1])])
print(f'per-biopsy GC slope: SD {slope.std():.3f}; differs by cohort (ANOVA F = {__import__("scipy").stats.f_oneway(*[slope[M.loc[s, "cohorte"].values == c] for c in M.loc[s].cohorte.unique()]).statistic:.1f})')
R.to_csv(OUTDIR + 'p5a_liver_gc.csv', index=False)
POSLAYERS_FILE_052_END
echo "  new      scripts/pilot_liver_gc.py"
mkdir -p "scripts"
cat > "scripts/pilot_spectra.py" << 'POSLAYERS_FILE_053_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.environ.get('LIVER_SPECTRA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Liver_Spectra')), 'scripts'))  # liver pipeline: github.com/Danpc11/Liver_Spectra
from common import *
X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl'))
s = M.index[M.estadio != 'Control']; grid = load_grid(); N = chrom_lengths(grid)
U = pd.read_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv')); U['chr'] = U.chr.astype(str)
def per_chrom_signal(vals, g, n):
    x = np.zeros(n); x[g.grid_index.values - 1] = vals; return x
def pgram(x): n = len(x); return np.abs(np.fft.rfft(x)[1:n // 2 + 1]) ** 2
def whiten(P, n):
    k = np.arange(1, len(P) + 1); lx = np.log10(k / n); b = np.polyfit(lx, np.log10(P + 1e-12), 1); w = P / 10 ** (b[0] * lx + b[1]); return w / w.mean()
EUL = 0.5772156649
rows = []; stat_tot = dyn_tot = 0.0
for c in CHR:
    g = keep[keep.chr == c].sort_values('grid_index'); n = N[c]
    Y = A.loc[g.gene_id, s].values                      # genes x samples, batch-corrected log expression
    Yc = Y - Y.mean(0, keepdims=True)                   # centred across genes within each sample (as in the pipeline)
    mu = Yc.mean(1); D = Yc - mu[:, None]               # static landscape and per-sample deviations
    stat_tot += np.sum(mu ** 2) * Y.shape[1]; dyn_tot += np.sum(D ** 2)
    Wfull = np.array([whiten(pgram(per_chrom_signal(Yc[:, i], g, n)), n) for i in range(Y.shape[1])])
    Wmu = whiten(pgram(per_chrom_signal(mu, g, n)), n)
    Wdyn = np.array([whiten(pgram(per_chrom_signal(D[:, i], g, n)), n) for i in range(Y.shape[1])])
    cons = np.exp(np.log(Wfull).mean(0) + EUL); cons_dyn = np.exp(np.log(Wdyn).mean(0) + EUL)
    for k in range(len(cons)):
        rows.append({'chr': c, 'k': k + 1, 'period': n / (k + 1), 'cons_full': cons[k], 'w_static': Wmu[k], 'cons_dynamic': cons_dyn[k],
                     'frac_full_gt3': np.mean(Wfull[:, k] > 3), 'frac_dyn_gt3': np.mean(Wdyn[:, k] > 3)})
R = pd.DataFrame(rows).merge(U[['chr', 'k', 'universal']], on=['chr', 'k'], how='left')
R.to_csv(OUTDIR + 'p1_spectra.csv', index=False)
u = R.universal.fillna(False).astype(bool)
print('frequencies:', len(R), '| universal peaks (pipeline):', int(u.sum()))
print('Parseval share of total positional variance: static landscape %.1f%%, per-sample deviations %.1f%%' % (100 * stat_tot / (stat_tot + dyn_tot), 100 * dyn_tot / (stat_tot + dyn_tot)))
print('r(log consensus spectrum, log spectrum of the static landscape) = %.3f' % np.corrcoef(np.log(R.cons_full), np.log(R.w_static))[0, 1])
print('universal peaks with static-landscape power > 3x: %d / %d' % ((R.w_static[u] > 3).sum(), u.sum()))
print('frequencies universal in the DYNAMIC spectra (>3x in >=90%% of biopsies): %d' % (R.frac_dyn_gt3 >= 0.9).sum())
print('r(log consensus, log dynamic consensus) = %.3f' % np.corrcoef(np.log(R.cons_full), np.log(R.cons_dynamic))[0, 1])
POSLAYERS_FILE_053_END
echo "  new      scripts/pilot_spectra.py"
mkdir -p "scripts"
cat > "scripts/predict.py" << 'POSLAYERS_FILE_054_END'
"""Quantitative predictions per GTEx tissue.
(1) Isochore law: GC-induced correlation between genes i,j = Var(b) * gc_i * gc_j / (sd_i * sd_j), b = per-sample GC slope.
(2) eQTL law: correlation induced by a shared causal variant = sum over shared credible sets of
    P(same variant) * 2p(1-p) * (afc_1/2) * (afc_2/2) / (sd_1 * sd_2)   [afc: log2 allelic fold change; afc/2 = per-allele log2 effect]."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyannotables as pa, pyarrow.parquet as pq
CHR = [str(i) for i in range(1, 23)] + ['X']; norm = lambda s: s.lower().replace('-', '_')
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start', 'End']]; G38.columns = ['chr', 'start', 'end']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
SUS = {norm(os.path.basename(f).replace('_v11_eQTLs_SuSiE_summary.parquet', '')): f for f in glob.glob(DATA + 'gtex_eqtl/*_v11_eQTLs_SuSiE_summary.parquet')}
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
LAGS = [1, 2, 5, 10, 20, 30]
for f in sys.argv[1:]:
    t = norm(os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', ''))
    if t == 'kidney_medulla' or os.path.exists(OUTDIR + f'pred_pairs_{t}.csv.gz'): continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
    g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std()
    # ---- (1) isochore law
    b = np.array([np.polyfit(gcz, D[:, i], 1)[0] for i in range(D.shape[1])]); varb = b.var(); sd_raw = D.std(1)
    Dg = resid_genes(D, np.column_stack([gcz, gcz ** 2]))
    def lagcorr(M, L):
        Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9); return np.mean(np.concatenate([np.mean(Z[i[:-L]] * Z[i[L:]], 1) for c in CHR for i in [np.where(chrs == c)[0]] if len(i) > L + 5]))
    law = {'tissue': t, 'var_GC_slope': varb}
    for L in LAGS:
        pr = np.concatenate([varb * gcz[i[:-L]] * gcz[i[L:]] / (sd_raw[i[:-L]] * sd_raw[i[L:]]) for c in CHR for i in [np.where(chrs == c)[0]] if len(i) > L + 5]).mean()
        law[f'pred_L{L}'] = pr; law[f'obs_gc_component_L{L}'] = lagcorr(D, L) - lagcorr(Dg, L)
    pd.DataFrame([law]).to_csv(OUTDIR + 'isochore_law.csv', mode='a', header=not os.path.exists(OUTDIR + 'isochore_law.csv'), index=False)
    # ---- (2) eQTL law, on GC- and technically corrected expression (as for the cis layer)
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    R = resid_samples(Dg, tech) if tech.shape[1] < Dg.shape[1] - 10 else Dg; sd = R.std(1); Z = (R - R.mean(1, keepdims=True)) / (sd[:, None] + 1e-9)
    same = chrs[:-1] == chrs[1:]; ii = np.where(same)[0]; pos = {gid: k for k, gid in enumerate(g.index)}
    if t in SUS:
        E = pq.read_table(SUS[t], columns=['phenotype_id', 'variant_id', 'pip', 'cs_id', 'af', 'afc']).to_pandas(); E['gid'] = E.phenotype_id.str.split('.').str[0]; E = E[E.gid.isin(pos)]
        CS = {gg: [(dict(zip(c.variant_id, c.pip)), dict(zip(c.variant_id, c.afc)), dict(zip(c.variant_id, c.af))) for _, c in d.groupby('cs_id')] for gg, d in E.groupby('gid')}
        rows = []
        for k in ii:
            g1, g2 = g.index[k], g.index[k + 1]; A, B = CS.get(g1), CS.get(g2)
            if A is None or B is None: continue
            cov = 0.0; pmax = 0.0
            for p1, f1, af1 in A:
                for p2, f2, af2 in B:
                    sh = [v for v in p1.keys() & p2.keys()]
                    if not sh: continue
                    ps = sum(p1[v] * p2[v] for v in sh); pmax = max(pmax, ps)
                    vs = [v for v in sorted(sh, key=lambda v: -p1[v] * p2[v]) if np.isfinite(f1.get(v, np.nan)) and np.isfinite(f2.get(v, np.nan))]
                    if not vs: continue
                    v = vs[0]; p = af1[v]; cov += ps * 2 * p * (1 - p) * (f1[v] / 2) * (f2[v] / 2)
            if pmax < 0.1: continue
            rows.append({'g1': g1, 'g2': g2, 'p_coloc': pmax, 'pred_r': cov / (sd[k] * sd[k + 1]), 'obs_r': float(np.mean(Z[k] * Z[k + 1])), 'dist': g.start.values[k + 1] - g.end.values[k]})
        P = pd.DataFrame(rows); P['tissue'] = t; P.to_csv(OUTDIR + f'pred_pairs_{t}.csv.gz', index=False)
        print(t, f'varb {varb:.4f} | pred/obs GC L20 {law["pred_L20"]:.4f}/{law["obs_gc_component_L20"]:.4f} | coloc pairs {len(P)} | r(pred,obs) {P[["pred_r", "obs_r"]].corr().iloc[0, 1]:.2f}', flush=True)
POSLAYERS_FILE_054_END
echo "  updated  scripts/predict.py"
mkdir -p "scripts"
cat > "scripts/predict_first.py" << 'POSLAYERS_FILE_055_END'
"""STEP 1 — predictions fixed from GTEx whole blood only (baseline cis coupling), blind to the edited-cell data."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import pandas as pd, numpy as np, pyannotables as pa
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
P = pd.read_csv(OUTDIR + 'pairs_whole_blood.csv.gz')
P['s1'] = sym.reindex(P.g1).values; P['s2'] = sym.reindex(P.g2).values
def neighbours(target, k=6):
    """genes within k positions of the target on the GTEx-blood gene order, with the coupling of the intervening pairs"""
    i1 = P.index[(P.s1 == target) | (P.s2 == target)]
    if not len(i1): return pd.DataFrame()
    rows = []
    for i in i1:
        r = P.loc[i]; nb = r.s2 if r.s1 == target else r.s1
        rows.append({'target': target, 'neighbour': nb, 'steps': 1, 'coupling_r': r.r, 'dist_kb': r.dist / 1e3})
    # two steps away: product of consecutive couplings (expected transmitted coupling)
    for i in i1:
        r = P.loc[i]; mid = r.s2 if r.s1 == target else r.s1
        j = P.index[((P.s1 == mid) | (P.s2 == mid))]
        for jj in j:
            r2 = P.loc[jj]; nb2 = r2.s2 if r2.s1 == mid else r2.s1
            if nb2 in (target, mid): continue
            rows.append({'target': target, 'neighbour': nb2, 'steps': 2, 'coupling_r': r.r * r2.r, 'dist_kb': np.nan})
    return pd.DataFrame(rows)
pred = pd.concat([neighbours(t) for t in ['BCL11A', 'HBG1', 'HBG2', 'HBB', 'HBD']]).drop_duplicates(['target', 'neighbour'])
pred['predicted_direction'] = np.where(pred.coupling_r > 0, 'same as target', 'opposite to target')
pred = pred.sort_values(['target', 'coupling_r'], ascending=[True, False])
pred.to_csv(OUTDIR + 'predictions_from_GTEx_blood.csv', index=False)
pd.set_option('display.width', 200); print(pred.round(3).to_string(index=False))
print('\nblood-wide reference: median |r| %.3f, 90th percentile %.3f' % (P.r.abs().median(), P.r.quantile(0.9)))
POSLAYERS_FILE_055_END
echo "  updated  scripts/predict_first.py"
mkdir -p "scripts"
cat > "scripts/replicate.py" << 'POSLAYERS_FILE_056_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm
from lib_tumour import coding, prepare, scores, trunc

COH = ['STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']; rows = []
for c in ['LAML', 'GBM']:
    D, chrs, samp, cov = prepare(c, codes=('01', '03')); e = scores(D, chrs)
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(COH)].Sample_ID)).astype(float)
        Xc = cov[['cna_burden', 'log_tmb', 'immune', 'stromal', 'subtype']]; Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join((Xc - Xc.mean()) / Xc.std()))
        base = e[m == 0].mean(); row = {'cohort': c, 'class': mclass, 'n_tumours': len(e), 'n_mutant': int(m.sum()), 'raw_pct': 100 * (e[m == 1].mean() - base) / base if m.sum() else np.nan}
        if m.sum() >= 3:
            fit = sm.OLS(e, Xd).fit(cov_type='HC3'); r = e - sm.OLS(e, sm.add_constant(((Xc - Xc.mean()) / Xc.std()).values)).fit().fittedvalues
            obs = r[m == 1].mean() - r[m == 0].mean(); rng = np.random.default_rng(0); null = np.array([(lambda p: r[p == 1].mean() - r[p == 0].mean())(rng.permutation(m)) for _ in range(10000)])
            row.update({'adj_pct': 100 * fit.params['mutant'] / base, 'ci_low': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high': 100 * fit.conf_int().loc['mutant', 1] / base, 'p_perm': (np.sum(np.abs(null) >= abs(obs)) + 1) / 10001})
        rows.append(row)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'replication_LAML_GBM.csv', index=False); print(R.round(3).to_string(index=False))
POSLAYERS_FILE_056_END
echo "  updated  scripts/replicate.py"
mkdir -p "scripts"
cat > "scripts/robust_families.py" << 'POSLAYERS_FILE_057_END'
"""Gene families, classic clusters and read-through transcripts (Fig. 2d). Input: pairs_<tissue>.csv.gz from atlas.py."""
import glob, os, re, numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
stem = lambda s: (re.match(r'^([A-Z]+)\d', s).group(1) if re.match(r'^([A-Z]+)\d', s) else None)
CLUST = ('HIST', 'H2A', 'H2B', 'H3C', 'H4C', 'H1-', 'OR', 'KRT', 'KRTAP', 'PCDH', 'ZNF', 'UGT', 'CYP', 'HLA', 'IGH', 'IGK', 'IGL', 'TRB', 'TRA', 'KIR', 'LCE', 'SPRR', 'MT1', 'MT2',
         'DEF', 'SERPIN', 'APO', 'GST', 'CLEC', 'LILR', 'SIGLEC', 'TAS2R', 'PRAME', 'GAGE', 'MAGE', 'SPANX', 'USP17', 'NBPF', 'TBC1D3', 'GOLGA6', 'FAM90')
def readthrough(s): return ('-' in s) and not re.search(r'-(AS|DT|IT|OT|OS)\d*$', s) and not s.startswith(('HLA-', 'H1-', 'H2A', 'H2B', 'H3-', 'H4-'))
rows = []
for f in sorted(glob.glob(OUTDIR + 'pairs_*.csv.gz')):
    t = os.path.basename(f)[6:-7]; P = pd.read_csv(f); P = P[P.dist > 0].copy()
    s1 = sym.reindex(P.g1).fillna('').values; s2 = sym.reindex(P.g2).fillna('').values
    st1 = np.array([stem(x) for x in s1], dtype=object); st2 = np.array([stem(x) for x in s2], dtype=object)
    same_stem = (st1 == st2) & pd.notna(st1)
    clus = np.array([a.startswith(CLUST) or b.startswith(CLUST) for a, b in zip(s1, s2)])
    rt = np.array([readthrough(a) or readthrough(b) for a, b in zip(s1, s2)])
    keep = ~(same_stem | clus | rt); base = P.r.mean()
    rows.append({'tissue': t, 'pairs': len(P), 'excluded_pct': 100 * (1 - keep.mean()), 'r_all': base, 'r_same_stem': P.r[same_stem].mean(), 'r_clusters': P.r[clus].mean(),
                 'r_readthrough': P.r[rt].mean(), 'r_clean': P.r[keep].mean(), 'clean_vs_all_pct': 100 * (P.r[keep].mean() / base - 1),
                 'r_clean_le20kb': P.r[keep & (P.dist <= 2e4)].mean(), 'r_clean_gt100kb': P.r[keep & (P.dist > 1e5)].mean()})
pd.DataFrame(rows).to_csv(OUTDIR + 'robust_families.csv', index=False); print(len(rows), 'tissues')
POSLAYERS_FILE_057_END
echo "  new      scripts/robust_families.py"
mkdir -p "scripts"
cat > "scripts/spectral_form.py" << 'POSLAYERS_FILE_058_END'
"""Spectral form of the cis layer. Residual expression (GC + technical covariates removed), standardised per gene.
Route A: autocorrelation rho(L) across samples -> fit exponential(s) -> lambda_acf.
Route B: mean periodogram of the residual profiles along gene order -> fit discrete Lorentzian(s) + floor -> lambda_spec;
compare with a power law + floor. Theory: exponential ACF <=> Lorentzian spectrum with the same lambda, no discrete peaks."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, pyannotables as pa
from scipy.optimize import curve_fit
CHR = [str(i) for i in range(1, 23)] + ['X']; norm = lambda s: s.lower().replace('-', '_')
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
def lor1(f, a, lam, c):     # discrete Lorentzian (spectrum of an exponential ACF with length lam) + white floor
    q = np.exp(-1 / lam); return c + a * (1 - q ** 2) / (1 - 2 * q * np.cos(2 * np.pi * f) + q ** 2)
def lor2(f, a1, l1, a2, l2, c): return lor1(f, a1, l1, 0) + lor1(f, a2, l2, 0) + c
def plaw(f, k, al, c): return c + k * f ** (-al)
OUT = OUTDIR + 'lorentz_v2.csv'
for t in sys.argv[1:]:
    f = DATA + f'gtex/gene_reads_adult_gtex_v11_{t}_gct.gz'
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
    g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); D = resid_genes(D, np.column_stack([gcz, gcz ** 2]))
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    if tech.shape[1] < D.shape[1] - 10: D = resid_samples(D, tech)
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    # route A: autocorrelation
    lags = np.arange(1, 61); rho = np.array([np.mean(np.concatenate([np.mean(Z[i[:-L]] * Z[i[L:]], 1) for c in CHR for i in [np.where(chrs == c)[0]] if len(i) > L + 5])) for L in lags])
    fl = rho[40:].mean(); y = rho - fl
    e1 = curve_fit(lambda L, a, lam: a * np.exp(-L / lam), lags[:30], y[:30], p0=[0.15, 5], maxfev=20000)[0]
    e2 = curve_fit(lambda L, a1, l1, a2, l2: a1 * np.exp(-L / l1) + a2 * np.exp(-L / l2), lags[:30], y[:30], p0=[0.1, 1, 0.05, 10], bounds=([0, 0.2, 0, 2], [1, 5, 1, 100]), maxfev=20000)[0]
    rss1 = np.sum((y[:30] - e1[0] * np.exp(-lags[:30] / e1[1])) ** 2); rss2 = np.sum((y[:30] - (e2[0] * np.exp(-lags[:30] / e2[1]) + e2[2] * np.exp(-lags[:30] / e2[3]))) ** 2)
    # route B: mean periodogram of residual profiles (per chromosome, normalised by length), on a common log-frequency grid
    bins = np.logspace(np.log10(1 / 1500), np.log10(0.5), 61); acc = np.zeros(60); cnt = np.zeros(60)
    for c in CHR:
        i = np.where(chrs == c)[0]; n = len(i)
        if n < 200: continue
        P = np.mean(np.abs(np.fft.rfft(Z[i], axis=0)[1:n // 2 + 1]) ** 2, axis=1) / n; fr = np.arange(1, n // 2 + 1) / n
        k = np.digitize(fr, bins) - 1; ok = (k >= 0) & (k < 60); np.add.at(acc, k[ok], P[ok]); np.add.at(cnt, k[ok], 1)
    m = cnt > 0; fx = np.sqrt(bins[:-1] * bins[1:])[m]; S = acc[m] / cnt[m]
    # Review point 7: fit and evaluate on the SAME scale. The log of a bin-averaged periodogram has variance ~ 1/n_bin
    # (n_bin = number of Fourier frequencies averaged in the bin), so fit log S with weights sqrt(n_bin) and compute the
    # Gaussian AIC from the same weighted residual sum of squares.
    wcnt = cnt[m]
    def fit(fun, p0, bnds):
        lf = lambda f, *q: np.log(np.maximum(fun(f, *q), 1e-12))
        p = curve_fit(lf, fx, np.log(S), p0=p0, bounds=bnds, sigma=1 / np.sqrt(wcnt), maxfev=200000)[0]
        r = np.log(S) - lf(fx, *p); return p, np.sum(wcnt * r ** 2)
    pL1, rL1 = fit(lor1, [1, 5, 0.8], ([0, 0.3, 0], [100, 500, 10])); pL2, rL2 = fit(lor2, [0.3, 1, 0.5, 15, 0.8], ([0, 0.2, 0, 2, 0], [100, 5, 100, 500, 10])); pP, rP = fit(plaw, [0.01, 0.5, 0.8], ([0, 0, 0], [10, 3, 10]))
    nf = len(fx); aic = lambda rss, k: nf * np.log(rss / nf) + 2 * k + 2 * k * (k + 1) / (nf - k - 1)   # AICc
    resid = S / lor2(fx, *pL2); row = {'tissue': t, 'samples': Z.shape[1], 'lam_acf_1exp': e1[1], 'lam_acf_short': e2[1], 'lam_acf_long': e2[3], 'acf_rss_1exp': rss1, 'acf_rss_2exp': rss2,
        'lam_spec_1lor': pL1[1], 'lam_spec_short': pL2[1], 'lam_spec_long': pL2[3], 'AIC_lor1': aic(rL1, 3), 'AIC_lor2': aic(rL2, 5), 'AIC_powerlaw': aic(rP, 3), 'powerlaw_alpha': pP[1],
        'max_peak_over_fit': resid.max(), 'n_bins_over_1.5x_fit': int((resid > 1.5).sum())}
    pd.DataFrame([row]).to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
    np.save(OUTDIR + f'spec_v2_{t}.npy', np.vstack([fx, S, lor2(fx, *pL2), plaw(fx, *pP)]))
    print(t, {k: round(v, 3) if isinstance(v, float) else v for k, v in row.items() if k != 'tissue'}, flush=True)
POSLAYERS_FILE_058_END
echo "  new      scripts/spectral_form.py"
mkdir -p "scripts"
cat > "scripts/tad_bootstrap.py" << 'POSLAYERS_FILE_059_END'
"""Review point 9: genomic-block bootstrap (10-Mb blocks) of the same-TAD effect in each tissue, with the model of Fig. 4b:
r ~ same_tad + distance bin + orientation, on adjacent pairs with both genes in a TAD."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import os, glob, numpy as np, pandas as pd, pyannotables as pa
from lib_tad import MATCH, tads, assign
TADA = {n: assign(tads(n)) for n in sorted(set(MATCH.values()))}
BINS = [0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]; G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
rng = np.random.default_rng(5); rows = []
for t, tn in MATCH.items():
    f = [x for x in glob.glob(OUTDIR + 'pairs_*.csv.gz') if os.path.basename(x)[6:-7].lower().replace('-', '_') == t][0]
    d = pd.read_csv(f); d = d[d.dist > 0].copy(); a = TADA[tn]; d['t1'] = a.reindex(d.g1).values; d['t2'] = a.reindex(d.g2).values
    d = d[(d.t1 >= 0) & (d.t2 >= 0)].reset_index(drop=True); d['same'] = (d.t1 == d.t2).astype(float)
    X = np.column_stack([d.same.values, pd.get_dummies(pd.cut(d.dist, BINS).astype(str)).values.astype(float), pd.get_dummies(d.orientation, drop_first=True).values.astype(float)])
    y = d.r.values; est = np.linalg.lstsq(X, y, rcond=None)[0][0]
    blk = (G.Chromosome.reindex(d.g1).astype(str).values + ':' + (G.Start.reindex(d.g1).fillna(0).values // 1e7).astype(int).astype(str))
    ub, inv = np.unique(blk, return_inverse=True); grp = [np.where(inv == k)[0] for k in range(len(ub))]
    bs = np.array([np.linalg.lstsq(X[ix], y[ix], rcond=None)[0][0] for ix in (np.concatenate([grp[k] for k in rng.integers(0, len(grp), len(grp))]) for _ in range(500))])
    rows.append({'tissue': t, 'pairs': len(d), 'blocks': len(ub), 'b_same_tad': est, 'ci_low': np.percentile(bs, 2.5), 'ci_high': np.percentile(bs, 97.5), 'boot_se': bs.std(),
                 'p_boot_normal': 2 * __import__('scipy.stats', fromlist=['norm']).norm.sf(abs(est / bs.std()))})
    print(rows[-1]['tissue'], round(est, 4), round(rows[-1]['ci_low'], 4), round(rows[-1]['ci_high'], 4), flush=True)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_block_bootstrap.csv', index=False)
print(f"CI excludes 0 in {(R.ci_low > 0).sum()}/{len(R)} tissues; median effect {R.b_same_tad.median():.3f}")
POSLAYERS_FILE_059_END
echo "  new      scripts/tad_bootstrap.py"
mkdir -p "scripts"
cat > "scripts/tad_test.py" << 'POSLAYERS_FILE_060_END'
"""Does the cis layer respect TAD boundaries? Adjacent gene pairs within one TAD vs across a boundary, at matched distance,
using the tissue's own TAD partition (McArthur & Capra 20-bin landscape, hg19; genes placed with GRCh37 coordinates)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import glob, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
from lib_tad import CHR, MATCH, tads, assign
TADA = {n: assign(tads(n)) for n in sorted(set(MATCH.values()))}
norm = lambda s: s.lower().replace('-', '_'); BINS = [0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
rows, spec = [], []
for t, tn in MATCH.items():
    d = P[t].merge(S[t][['g1', 'g2', 'egene1', 'egene2', 'share_any', 'share_same']], on=['g1', 'g2']); d = d[d.dist > 0].copy(); d['bin'] = pd.cut(d.dist, BINS).astype(str)
    def fit(tadname, sub=None):
        a = TADA[tadname]; x = d.copy(); x['t1'] = a.reindex(x.g1).values; x['t2'] = a.reindex(x.g2).values; x = x[(x.t1 >= 0) & (x.t2 >= 0)]
        x['same_tad'] = (x.t1 == x.t2).astype(int)
        if sub is not None: x = x[sub(x)]
        m = smf.ols('r ~ same_tad + C(bin) + C(orientation)', data=x).fit(); return m.params['same_tad'], m.pvalues['same_tad'], len(x), x.same_tad.mean()
    b_all, p_all, n_all, f_same = fit(tn)
    b_ng, p_ng, n_ng, _ = fit(tn, sub=lambda x: ~x.share_any)                      # no shared eQTL: non-genetic component
    others = [fit(o)[0] for o in TADA if o != tn]
    rows.append({'tissue': t, 'tad_map': tn, 'pairs': n_all, 'frac_same_tad': f_same, 'b_same_tad': b_all, 'p': p_all, 'b_same_tad_no_shared_eQTL': b_ng, 'p_no_shared_eQTL': p_ng,
                 'b_same_tad_other_maps_median': np.median(others), 'own_map_rank_among_maps': 1 + sum(o > b_all for o in others), 'n_maps': 1 + len(others)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_cis_test.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
print(f"\nsame-TAD effect > 0 in {(R.b_same_tad > 0).sum()}/{len(R)} (P<0.05 in {(R.p < 0.05).sum()}); median {R.b_same_tad.median():.3f}; without shared eQTL median {R.b_same_tad_no_shared_eQTL.median():.3f}; own map > median of other maps in {(R.b_same_tad > R.b_same_tad_other_maps_median).sum()}/{len(R)}")
POSLAYERS_FILE_060_END
echo "  updated  scripts/tad_test.py"
mkdir -p "scripts"
cat > "scripts/tad_tumour.py" << 'POSLAYERS_FILE_061_END'
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyannotables as pa, statsmodels.api as sm
from lib_tumour import CHR, coding, prepare, trunc

B = pd.read_csv(glob.glob(DATA + 'TAD-full/*/data/boundariesByStability/100kbBookendBoundaries_mainText/100kbBookendBoundaries_byStability.bed')[0], sep='\t')
B['chr'] = B.chr.str.replace('chr', ''); B['mid'] = (B['loc'] + B['loc2']) / 2
G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()]; mid37 = ((G37.Start + G37.End) / 2); chr37 = G37.Chromosome.astype(str)
def pair_classes(gids, chrs):
    i = np.where(chrs[:-1] == chrs[1:])[0]; a = gids[i]; b = gids[i + 1]
    m1 = mid37.reindex(a).values; m2 = mid37.reindex(b).values; c1 = chr37.reindex(a).values
    cls = np.full(len(i), 'na', dtype=object)
    for c, bb in B.groupby('chr'):
        sel = np.where((c1 == c) & np.isfinite(m1) & np.isfinite(m2))[0]
        if not len(sel): continue
        lo = np.minimum(m1[sel], m2[sel]); hi = np.maximum(m1[sel], m2[sel]); bm = bb['mid'].values; st = bb.stability_percentile.values
        for k, (l, h) in zip(sel, zip(lo, hi)):
            inside = (bm > l) & (bm < h)
            cls[k] = 'within' if not inside.any() else ('stable' if st[inside].max() >= 0.75 else 'weak')
    return i, cls
def score(D, chrs, i, mask):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9); near = (Z[i[mask]] * Z[i[mask] + 1]).mean(0)
    far = np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in (20, 30)], 0)
    return near - far
rows = []
for c, gm in [('BLCA', ['STAG2']), ('UCEC', ['CTCF'])]:
    D, chrs, samp, cov, gids = prepare(c, return_gids=True); i, cls = pair_classes(gids, chrs)
    print(c, pd.Series(cls).value_counts().to_dict(), flush=True)
    S = {k: score(D, chrs, i, cls == k) for k in ('within', 'stable')}
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, hq); m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(gm)].Sample_ID)).astype(float)[keep]
            Xc = cov[keep]; Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join((Xc - Xc.mean()) / Xc.std()))
            res = {'cohort': c, 'gene': gm[0], 'mutation_class': mclass, 'tmb_q': hq, 'n_mut': int(m.sum())}
            for k in ('within', 'stable'):
                y = S[k][keep]; f = sm.OLS(y, Xd).fit(cov_type='HC3'); base = y[m == 0].mean()
                res[f'{k}_wt'] = base; res[f'{k}_effect_pct'] = 100 * f.params['mutant'] / abs(base); res[f'{k}_p'] = f.pvalues['mutant']
            dd = S['within'][keep] - S['stable'][keep]; f = sm.OLS(dd, Xd).fit(cov_type='HC3'); res['within_minus_stable_effect'] = f.params['mutant']; res['p_difference'] = f.pvalues['mutant']
            rows.append(res)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_tumour_results.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
POSLAYERS_FILE_061_END
echo "  updated  scripts/tad_tumour.py"
mkdir -p "scripts"
cat > "scripts/test_edits.py" << 'POSLAYERS_FILE_062_END'
"""STEP 2 — test the fixed predictions in the edited erythroblasts (GSE264491)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import pandas as pd, numpy as np, pyannotables as pa
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
C = pd.read_csv(DATA + 'edit/GSE264491_merged_counts.csv.gz', index_col=0)
grp = {'unedited': C.columns[:3], 'BCL11A': C.columns[3:6], 'HBG': C.columns[6:9]}
C = C[C.median(axis=1) >= 10]; cpm = np.log2(C / C.sum(0) * 1e6 + 1)
def lfc(g): return cpm[grp[g]].mean(1) - cpm[grp['unedited']].mean(1)
def tstat(g):
    a = cpm[grp[g]].values; b = cpm[grp['unedited']].values
    return (a.mean(1) - b.mean(1)) / np.sqrt(a.var(1, ddof=1) / 3 + b.var(1, ddof=1) / 3 + 1e-9)
R = pd.DataFrame({'sym': sym.reindex(C.index).values, 'lfc_BCL11A': lfc('BCL11A'), 'lfc_HBG': lfc('HBG'), 't_BCL11A': tstat('BCL11A'), 't_HBG': tstat('HBG')}, index=C.index)
print('== on-target effects'); print(R[R.sym.isin(['BCL11A', 'HBG1', 'HBG2', 'HBB', 'HBD', 'HBE1', 'HBBP1'])].round(2).to_string(index=False))
pred = pd.read_csv(OUTDIR + 'predictions_from_GTEx_blood.csv'); s2i = {v: k for k, v in sym.items() if v in set(pred.neighbour)}
pred['gid'] = pred.neighbour.map(s2i); pred = pred[pred.gid.isin(R.index)]
pred['lfc_BCL11A'] = R.lfc_BCL11A.reindex(pred.gid).values; pred['lfc_HBG'] = R.lfc_HBG.reindex(pred.gid).values
print('\n== predicted neighbours, observed log2 fold change')
print(pred[['target', 'neighbour', 'steps', 'coupling_r', 'lfc_BCL11A', 'lfc_HBG']].round(3).to_string(index=False))
# genome-wide: do neighbours of the edited locus move with their baseline coupling?
P = pd.read_csv(OUTDIR + 'pairs_whole_blood.csv.gz'); P = P[P.g1.isin(R.index) & P.g2.isin(R.index)]
for edit, tcol in [('BCL11A', 't_BCL11A'), ('HBG', 't_HBG')]:
    d = P.assign(t1=R[tcol].reindex(P.g1).values, t2=R[tcol].reindex(P.g2).values).dropna(subset=['t1', 't2'])
    d['bin'] = pd.qcut(d.r, 5)
    g = d.groupby('bin', observed=True).apply(lambda x: pd.Series({'n': len(x), 'baseline_r': x.r.mean(), 'corr_of_edit_effects': stats.pearsonr(x.t1, x.t2)[0]}), include_groups=False)
    print(f'\n== {edit} edit: concordance of neighbour responses by baseline coupling'); print(g.round(3).to_string())
    print('  slope of concordance on baseline coupling: %.2f' % np.polyfit(d.r, (d.t1 - d.t1.mean()) / d.t1.std() * (d.t2 - d.t2.mean()) / d.t2.std(), 1)[0])
R.to_csv(OUTDIR + 'edit_effects.csv'); pred.to_csv(OUTDIR + 'predictions_tested.csv', index=False)
# cis vs trans: distance of responding genes from the edited locus
loc = {'BCL11A': ('2', 60450000), 'HBG': ('11', 5250000)}
pos = pd.DataFrame({'chr': G.Chromosome.astype(str).reindex(R.index), 'mid': ((G.Start + G.End) / 2).reindex(R.index)})
for edit, tcol in [('BCL11A', 't_BCL11A'), ('HBG', 't_HBG')]:
    c, p0 = loc[edit]; near = (pos.chr == c) & ((pos.mid - p0).abs() < 2e6)
    big = R[tcol].abs() > 4
    print(f'\n{edit} edit: responding genes (|t|>4) = {int(big.sum())}; of these within 2 Mb of the edited site = {int((big & near).sum())}; genes within 2 Mb = {int(near.sum())}')
    if (big & near).sum(): print('  ', R.sym[big & near].tolist())
POSLAYERS_FILE_062_END
echo "  updated  scripts/test_edits.py"
mkdir -p "scripts"
cat > "scripts/theory_sim.py" << 'POSLAYERS_FILE_063_END'
"""Simulations for the positional-layer theory.
x[s,j] = mu[j] + beta[j]*z[s] + c[s,j] + b[s]*gc[j] + e[s,j]
  mu: tissue landscape (no positional structure); c: cis co-regulation (moving average of local regulators, decay length lam);
  b*gc: per-sample technical GC bias acting on isochore-clustered GC; beta*z: stage effect, locally smoothed.
Checks: (1) spectral decomposition identity; (2) 'universal peaks' with real vs random gene order; (3) GC creates a false
domain scale; (4) recovery of the cis decay length by naive vs corrected estimators."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd
rng = np.random.default_rng(42)
def ar1(n, phi, rng): x = np.zeros(n); x[0] = rng.normal()
for _ in range(0): pass
def ar1(n, phi, rng):
    x = np.empty(n); x[0] = rng.normal()
    for i in range(1, n): x[i] = phi * x[i - 1] + np.sqrt(1 - phi ** 2) * rng.normal()
    return x
def simulate(n_chr=4, n=1000, S=200, lam=3.0, cis_sd=0.6, gc_sd=0.35, stage_sd=0.15, land_sd=2.0, seed=0, cis=True, gc=True, stage=True):
    r = np.random.default_rng(seed); chrs = np.repeat(np.arange(n_chr), n); N = n_chr * n
    mu = r.normal(0, land_sd, N)                                    # landscape: independent across positions
    GC = np.concatenate([ar1(n, 0.97, r) for _ in range(n_chr)])     # isochore-clustered GC
    z = r.integers(0, 6, S).astype(float)
    X = np.tile(mu, (S, 1)) + r.normal(0, 1, (S, N))
    if cis:
        k = np.arange(-int(6 * lam), int(6 * lam) + 1); w = np.exp(-np.abs(k) / lam); w /= np.sqrt(np.sum(w ** 2))
        for c in range(n_chr):
            sl = slice(c * n, (c + 1) * n); U = r.normal(0, 1, (S, n + 2 * len(k)))
            X[:, sl] += cis_sd * np.array([np.convolve(u, w, mode='same')[len(k):len(k) + n] for u in U])
    if gc: X += np.outer(r.normal(0, gc_sd, S), GC)
    if stage:
        beta = np.concatenate([np.convolve(r.normal(0, 1, n + 20), np.ones(3) / np.sqrt(3), mode='same')[10:10 + n] for _ in range(n_chr)]) * stage_sd
        X += np.outer(z, beta)
    return X, chrs, GC, mu, z
def pgram(x): n = len(x); return np.abs(np.fft.rfft(x)[1:n // 2 + 1]) ** 2
def whiten(P, n):
    k = np.arange(1, len(P) + 1); lx = np.log10(k / n); b = np.polyfit(lx, np.log10(P + 1e-12), 1); w = P / 10 ** (b[0] * lx + b[1]); return w / w.mean()
def universal(X, chrs, perm_seed=None):
    tot = 0
    for c in np.unique(chrs):
        idx = np.where(chrs == c)[0]
        if perm_seed is not None: idx = np.random.default_rng(perm_seed).permutation(idx)
        Xc = X[:, idx] - X[:, idx].mean(1, keepdims=True); W = np.array([whiten(pgram(x), len(idx)) for x in Xc]); tot += int((np.mean(W > 3, 0) >= 0.9).sum())
    return tot
LAGS = np.array([1, 2, 3, 5, 7, 10, 15, 20, 30, 50])
def lagprof(D, chrs):
    out = []
    for L in LAGS:
        v = []
        for c in np.unique(chrs):
            Z = D[:, chrs == c]; Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9); v.append(np.mean(Z[:, :-L] * Z[:, L:]))
        out.append(np.mean(v))
    return np.array(out)
def deviations(X, GC=None, z=None):
    D = X - X.mean(0)                                                 # remove the landscape (per-gene mean across samples)
    if z is not None: Zd = np.column_stack([np.ones(len(z)), z]); D = D - Zd @ np.linalg.lstsq(Zd, D, rcond=None)[0]
    if GC is not None: G = np.column_stack([np.ones(len(GC)), GC, GC ** 2]); D = D - (G @ np.linalg.lstsq(G, D.T, rcond=None)[0]).T
    return D
def decay_length(prof):
    m = prof > 0.01; b = np.polyfit(LAGS[m], np.log(prof[m]), 1) if m.sum() >= 3 else [np.nan]; return -1 / b[0]
if __name__ == '__main__':
    out = {}
    # (1) identity: mean periodogram = |FT(mu_hat)|^2 + FT of the lag-summed (circular) sample covariance
    X, chrs, GC, mu, z = simulate(seed=1); c0 = chrs == 0; Xc = X[:, c0]; n = Xc.shape[1]
    mbar = Xc.mean(0); Dv = Xc - mbar
    lhs = np.mean([np.abs(np.fft.fft(x)) ** 2 for x in Xc], 0)
    C = np.array([np.mean(np.sum(Dv * np.roll(Dv, -L, axis=1), 1)) for L in range(n)])
    rhs = np.abs(np.fft.fft(mbar)) ** 2 + np.real(np.fft.fft(C))
    out['identity_max_rel_error'] = float(np.max(np.abs(lhs - rhs) / lhs))
    out['landscape_share_of_power'] = float(np.sum(np.abs(np.fft.fft(mbar)) ** 2) / np.sum(lhs))
    # (2) universal peaks with real vs random order, with and without positional structure
    for name, kw in [('landscape only', dict(cis=False, gc=False, stage=False)), ('full model', {})]:
        Xs, ch, *_ = simulate(seed=2, **kw); out[f'universal_{name}_real'] = universal(Xs, ch); out[f'universal_{name}_random'] = universal(Xs, ch, perm_seed=5)
    # (3) GC creates a false domain scale; correction recovers the cis-only truth
    Xg, ch, GCg, _, zg = simulate(seed=3, cis=True, gc=True); Xt, cht, *_ = simulate(seed=3, cis=True, gc=False)
    prof_raw = lagprof(deviations(Xg, z=zg), ch); prof_corr = lagprof(deviations(Xg, GC=GCg, z=zg), ch); prof_truth = lagprof(deviations(Xt, z=zg), cht)
    pd.DataFrame({'lag': LAGS, 'raw': prof_raw, 'GC_corrected': prof_corr, 'truth_no_GC': prof_truth}).to_csv(OUTDIR + 'sim_lag_profiles.csv', index=False)
    # (4) recovery of the true cis decay length
    rec = []
    for lam in [1, 2, 3, 5, 8]:
        for rep in range(3):
            Xr, chr_, GCr, _, zr = simulate(seed=100 + 10 * lam + rep, lam=lam)
            rec.append({'true_lambda': lam, 'rep': rep, 'naive': decay_length(lagprof(deviations(Xr, z=zr), chr_)), 'corrected': decay_length(lagprof(deviations(Xr, GC=GCr, z=zr), chr_))})
    REC = pd.DataFrame(rec); REC.to_csv(OUTDIR + 'sim_decay_recovery.csv', index=False)
    pd.Series(out).to_csv(OUTDIR + 'sim_summary.csv')
    print(pd.Series(out).round(4).to_string()); print(pd.DataFrame({'lag': LAGS, 'raw': prof_raw, 'GC_corr': prof_corr, 'truth': prof_truth}).round(3).to_string(index=False))
    print(REC.groupby('true_lambda')[['naive', 'corrected']].median().round(2).to_string())
POSLAYERS_FILE_063_END
echo "  updated  scripts/theory_sim.py"
mkdir -p "src/poslayers"
cat > "src/poslayers/config.py" << 'POSLAYERS_FILE_064_END'
"""Where inputs, intermediate tables and figures live. Every script imports these three strings.

    POSLAYERS_DATA     raw public inputs (default: data/)        -- see data/MANIFEST.md
    POSLAYERS_RESULTS  every table the scripts produce (default: results/)
    POSLAYERS_FIGS     figure files (default: figures/output/)

Each value ends with '/', so scripts write DATA + 'gtex/...' or OUTDIR + f'pairs_{t}.csv.gz'.
"""
import os

def _dir(var, default):
    p = os.path.abspath(os.environ.get(var, default)).rstrip('/') + '/'
    return p

DATA = _dir('POSLAYERS_DATA', 'data')
OUTDIR = _dir('POSLAYERS_RESULTS', 'results')
FIGDIR = _dir('POSLAYERS_FIGS', 'figures/output')
for _d in (OUTDIR, FIGDIR):
    os.makedirs(_d, exist_ok=True)
POSLAYERS_FILE_064_END
echo "  new      src/poslayers/config.py"
mkdir -p "src/poslayers"
cat > "src/poslayers/decompose.py" << 'POSLAYERS_FILE_065_END'
"""Exact decomposition of a positional expression matrix (Supplementary Note 1, sections 1-2)."""
from __future__ import annotations
import numpy as np


def lag_covariance(Y: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Circular lag covariance C(L) of deviations Y (samples x genes), averaged over samples.

    C(L) = (1/S) sum_s sum_j y(s, j) y(s, j + L), the quantity whose Fourier transform
    is the covariance spectrum S(f).

    The sum is CIRCULAR: index j + L wraps around the end of the array, which is what
    makes the Wiener-Khinchin identity in periodogram_identity() exact. Pass one
    chromosome at a time. For a linear (non-wrapping) lag profile at a few short lags,
    use lag_profile_linear() instead; the two agree to O(L / n_genes).
    """
    Y = np.asarray(Y, float)
    n = Y.shape[1]
    max_lag = n if max_lag is None else min(max_lag, n)
    F = np.fft.rfft(Y, axis=1)
    C = np.fft.irfft((F * np.conj(F)).real, n=n, axis=1).mean(0)
    return C[:max_lag]


def periodogram_identity(X: np.ndarray) -> dict:
    """Verify  mean_s |FT(x_s)|^2  ==  |FT(xbar)|^2 + FT(C).

    X is samples x genes, genes ordered along a chromosome. Returns the three spectra and
    the maximum relative error of the identity, which is an algebraic identity and so
    should be at machine precision.
    """
    X = np.asarray(X, float)
    xbar = X.mean(0)
    Y = X - xbar
    lhs = (np.abs(np.fft.rfft(X, axis=1)) ** 2).mean(0)
    land = np.abs(np.fft.rfft(xbar)) ** 2
    cov = (np.abs(np.fft.rfft(Y, axis=1)) ** 2).mean(0)
    err = np.max(np.abs(lhs - land - cov) / (np.abs(lhs) + 1e-12))
    return {"mean_periodogram": lhs, "landscape": land, "covariance": cov, "max_relative_error": float(err)}


def landscape_share(X: np.ndarray) -> float:
    """Fraction of the spectral power of a single sample carried by the tissue landscape."""
    X = np.asarray(X, float)
    xbar = X.mean(0)
    land = (xbar ** 2).sum() * X.shape[0]
    resid = ((X - xbar) ** 2).sum()
    return float(land / (land + resid))



def _gc_basis(gc, degree):
    """Orthonormal basis of [1, g, ..., g^degree] with its EFFECTIVE rank.

    Review point 5: with constant GC (or fewer distinct GC values than the degree), the polynomial columns are collinear.
    A plain QR still returns `degree + 1` columns, and projecting on the spurious ones removes signal unrelated to GC.
    We keep only directions whose singular value is non-negligible."""
    g = np.asarray(gc, float)
    sd = g.std()
    g = (g - g.mean()) / sd if sd > 0 else np.zeros_like(g)
    D = np.vstack([g ** k for k in range(degree + 1)]).T
    U, s, _ = np.linalg.svd(D, full_matrices=False)
    keep = s > s[0] * 1e-10 * max(D.shape)
    return U[:, keep]


def gc_correct(X: np.ndarray, gc: np.ndarray, degree: int = 2) -> np.ndarray:
    """Remove a per-sample polynomial trend on gene GC content (the technical isochore layer).

    X is samples x genes. The basis includes the constant column, so each sample's mean across genes is also removed: a
    per-sample offset is not positional information, and rows come out centred. When GC is constant the basis has rank 1
    and only the per-sample mean is removed. On simulated data with no GC bias the correction costs 0.3% of the
    adjacent-gene correlation; when real regulation tracks GC it can remove much more at domain scale (see
    scripts/gc_correlated_sim.py), so it cannot separate technical bias from GC-associated biology.
    """
    X = np.asarray(X, float)
    if X.ndim != 2 or X.shape[1] != len(gc):
        raise ValueError(f'X must be samples x genes with {len(gc)} genes; got shape {X.shape}')
    Q = _gc_basis(gc, degree)
    return X - (X @ Q) @ Q.T


def gc_slopes(X: np.ndarray, gc: np.ndarray) -> np.ndarray:
    """Per-sample GC slope b_s (OLS on standardised GC with an intercept), whose variance drives the isochore law.
    Returns zeros when GC is constant, where the slope is not identifiable."""
    X = np.asarray(X, float)
    g = np.asarray(gc, float)
    if X.shape[1] != len(g):
        raise ValueError(f'X must be samples x genes with {len(g)} genes; got shape {X.shape}')
    if g.std() == 0:
        return np.zeros(X.shape[0])
    g = (g - g.mean()) / g.std()
    Y = X - X.mean(1, keepdims=True)
    return (Y @ g) / (g @ g)


def lag_profile_linear(Y: np.ndarray, n_genes: int, n_chrom: int, max_lag: int = 60) -> np.ndarray:
    """Non-circular lag correlation profile, averaged over chromosomes and samples.

    Y is samples x genes with the genes of n_chrom chromosomes of n_genes each, concatenated. This is the estimator used for
    the coupling atlas: genes standardised across samples, pairs taken within a chromosome without wrapping. Lags with no
    pairs (max_lag >= n_genes) are returned as NaN instead of failing.
    """
    Y = np.asarray(Y, float)
    if Y.ndim != 2 or Y.shape[1] != n_genes * n_chrom:
        raise ValueError(f'Y must be samples x (n_genes * n_chrom) = {n_genes * n_chrom} columns; got shape {Y.shape}')
    Z = (Y - Y.mean(0)) / (Y.std(0) + 1e-12)
    out = np.full(max_lag + 1, np.nan)
    for L in range(min(max_lag, n_genes - 1) + 1):
        acc, cnt = 0.0, 0
        for c in range(n_chrom):
            s = slice(c * n_genes, (c + 1) * n_genes)
            A, B = Z[:, s][:, :n_genes - L], Z[:, s][:, L:]
            acc += float(np.mean(A * B)) * A.shape[1]; cnt += A.shape[1]
        out[L] = acc / cnt
    return out
POSLAYERS_FILE_065_END
echo "  updated  src/poslayers/decompose.py"
mkdir -p "tests"
cat > "tests/test_audit.py" << 'POSLAYERS_FILE_066_END'
"""Regression tests for the problems found in the audit of the repository."""
import numpy as np
from poslayers import simulate_genome, lag_covariance, lag_profile_linear
from poslayers.decompose import gc_correct


def _fit_decay(prof):
    """Log-linear fit over the lags where the profile is still above 5% of its lag-1 value.

    A fixed lag window lets the noisy tail dominate when the decay is short: with a fixed
    window of 30 lags, a true length of 2 genes came back as 6. The relative window
    recovers lambda to within 3% for lambda >= 2, and returns nan below that, which is the
    honest answer for a decay of about one gene.
    """
    lags = np.arange(1, len(prof))
    v = prof[lags]
    ok = (v > 0.05 * v[0]) & (v > 1e-4)
    if ok.sum() < 4:
        return np.nan
    slope = np.polyfit(lags[ok], np.log(v[ok]), 1)[0]
    return -1 / slope if slope < 0 else np.nan


def test_kernel_autocorrelation_decays_with_lambda():
    """No calibration constant: the decay length of the autocorrelation IS lambda.

    This is the bug the audit found. The code used to divide by 1.56, which made the
    GC-corrected estimator undershoot by about a third.
    """
    for lam in (2.0, 3.0, 5.0, 8.0, 12.0):
        s = simulate_genome(n_genes=800, n_chrom=4, n_samples=300, decay_len=lam,
                            gc_bias_sd=0.0, noise_sd=0.0, cis_sd=1.0, seed=2)
        Y = s["X"] - s["X"].mean(0)
        fitted = _fit_decay(lag_profile_linear(Y, s["n_genes_per_chrom"], s["n_chrom"], 40))
        assert 0.85 * lam < fitted < 1.15 * lam, (lam, fitted)


def test_gc_correction_does_not_eat_real_signal():
    """With no GC bias at all, correcting costs almost nothing."""
    s = simulate_genome(n_genes=800, n_chrom=4, n_samples=200, decay_len=5.0,
                        gc_bias_sd=0.0, seed=6)
    Y = s["X"] - s["X"].mean(0)
    raw = lag_profile_linear(Y, 800, 4, 5)[1]
    cor = lag_profile_linear(gc_correct(Y, s["gc"]), 800, 4, 5)[1]
    assert cor > 0.95 * raw


def test_gc_correction_removes_a_pure_artefact():
    """With a GC bias and no cis layer, correction should flatten the lag profile."""
    s = simulate_genome(n_genes=800, n_chrom=4, n_samples=200, decay_len=3.0,
                        gc_bias_sd=0.6, cis_sd=0.0, noise_sd=0.3, seed=11)
    Y = s["X"] - s["X"].mean(0)
    raw = lag_profile_linear(Y, 800, 4, 20)
    cor = lag_profile_linear(gc_correct(Y, s["gc"]), 800, 4, 20)
    assert abs(cor[10]) < 0.5 * abs(raw[10])


def test_lag_covariance_is_circular_as_documented():
    """lag_covariance wraps; lag_profile_linear does not. Both are intended."""
    rng = np.random.default_rng(0)
    Y = rng.normal(size=(40, 120))
    n = Y.shape[1]
    circ = [float(np.mean([np.sum(Y[s] * np.roll(Y[s], -L)) for s in range(40)])) for L in range(4)]
    assert np.allclose(lag_covariance(Y, max_lag=4), circ)


def test_gc_correct_centres_each_sample():
    """The polynomial basis includes the constant column, so rows come out centred."""
    rng = np.random.default_rng(1)
    g = rng.normal(size=200)
    X = rng.normal(size=(20, 200)) + 7.0
    assert np.allclose(gc_correct(X, g).mean(1), 0.0, atol=1e-9)


def test_gc_correct_with_constant_gc_removes_only_the_sample_mean():
    """Review point 5: a constant GC track has rank 1; nothing beyond the per-sample mean may be removed."""
    rng = np.random.default_rng(0); X = rng.normal(size=(100, 10)); g = np.full(10, 0.4)
    assert np.allclose(gc_correct(X, g), X - X.mean(1, keepdims=True), atol=1e-10)


def test_gc_slopes_with_constant_gc_are_zero_not_nan():
    from poslayers.decompose import gc_slopes
    s = gc_slopes(np.random.default_rng(1).normal(size=(5, 8)), np.ones(8))
    assert np.all(np.isfinite(s)) and np.allclose(s, 0)


def test_lag_profile_handles_lags_beyond_the_chromosome():
    """Review point 5: max_lag >= n_genes used to raise; lags without pairs are now NaN."""
    p = lag_profile_linear(np.random.default_rng(2).normal(size=(20, 10)), n_genes=10, n_chrom=1, max_lag=12)
    assert np.isfinite(p[:10]).all() and np.isnan(p[10:]).all()


def test_shape_mismatch_is_an_error():
    import pytest
    with pytest.raises(ValueError):
        gc_correct(np.zeros((4, 5)), np.arange(6.0))


def test_gc_correction_removes_gc_tracking_biology_at_domain_scale():
    """Review point 12, documented limitation: when real regulation tracks GC, GC correction removes much of it at
    domain scale but little between neighbours."""
    s = simulate_genome(n_genes=800, n_chrom=4, n_samples=200, decay_len=5.0, gc_bias_sd=0.0, seed=50)
    rng = np.random.default_rng(90); gc = s["gc"]; cis = s["cis"]
    bio = np.outer(rng.normal(0, 1, cis.shape[0]), np.convolve(gc, np.ones(9) / 9, "same")); bio *= cis.std() / bio.std()
    new = np.sqrt(0.5) * cis + np.sqrt(0.5) * bio; Y = new - new.mean(0)
    truth = lag_profile_linear(Y, 800, 4, 10); kept = lag_profile_linear(gc_correct(Y, gc), 800, 4, 10)
    assert kept[10] / truth[10] < 0.7 and kept[1] / truth[1] > 0.8
POSLAYERS_FILE_066_END
echo "  updated  tests/test_audit.py"
mkdir -p "tests"
cat > "tests/test_laws.py" << 'POSLAYERS_FILE_067_END'
"""Tests of the decomposition identity and of the three laws."""
import numpy as np
import pytest
from poslayers import (periodogram_identity, landscape_share, lag_covariance,
                       isochore_law, eqtl_law, saturation_exponent, fit_contact_law,
                       simulate_genome)
from poslayers.decompose import gc_slopes, gc_correct


@pytest.fixture(scope="module")
def sim():
    return simulate_genome(n_genes=300, n_chrom=2, n_samples=120, decay_len=5.0, seed=7)


def test_identity_is_exact(sim):
    """|M(f)|^2 + S(f) equals the mean single-sample periodogram to machine precision."""
    assert periodogram_identity(sim["X"])["max_relative_error"] < 1e-8


def test_landscape_dominates(sim):
    assert 0.8 < landscape_share(sim["X"]) < 1.0


def test_lag_covariance_decays(sim):
    X = sim["X"]
    C = lag_covariance(X - X.mean(0), max_lag=40)
    assert C[0] > C[5] > C[20]


def test_isochore_law_recovers_a_simulated_bias():
    """With a known per-sample GC slope, the law predicts the induced correlation."""
    s = simulate_genome(n_genes=400, n_chrom=2, n_samples=200, decay_len=3.0,
                        gc_bias_sd=0.5, cis_sd=0.0, noise_sd=0.3, seed=3)
    X, gc = s["X"], s["gc"]
    Y = X - X.mean(0)
    var_b = float(np.var(gc_slopes(X, gc)))
    sd = Y.std(0)
    i, j = 0, 1
    pred = isochore_law(var_b, gc[i], gc[j], sd[i], sd[j])
    obs = float(np.mean(Y[:, i] * Y[:, j]) / (sd[i] * sd[j]))
    assert np.sign(pred) == np.sign(obs) or abs(obs) < 0.02
    assert abs(pred - obs) < 0.25


def test_gc_correction_removes_the_artefact():
    """With no cis layer, GC correction should leave no lag-1 correlation."""
    s = simulate_genome(n_genes=400, n_chrom=2, n_samples=200, decay_len=3.0,
                        gc_bias_sd=0.6, cis_sd=0.0, noise_sd=0.3, seed=11)
    X, gc = s["X"], s["gc"]
    Y = X - X.mean(0)
    lag1 = lambda M: float(np.mean((M[:, :-1] / (M.std(0)[:-1] + 1e-12)) *
                                   (M[:, 1:] / (M.std(0)[1:] + 1e-12))))
    assert abs(lag1(gc_correct(Y, gc))) < abs(lag1(Y))


def test_eqtl_law_sign_and_scale():
    """Opposite effects give negative coupling; the scale is 2p(1-p) b1 b2."""
    assert eqtl_law(0.3, 0.4, -0.5) < 0
    assert eqtl_law(0.3, 0.4, 0.5) > 0
    assert np.isclose(eqtl_law(0.5, 1.0, 1.0), 0.5)
    assert np.isclose(eqtl_law(0.3, 0.4, 0.5, p_same_causal=0.5),
                      0.5 * eqtl_law(0.3, 0.4, 0.5))


def test_saturation_exponent_is_one_minus_occupancy():
    """k -> 1 far from saturation, k -> 0 when saturated, k = 0.5 at c = K."""
    assert saturation_exponent(1e-3, 1.0) > 0.99
    assert saturation_exponent(1e3, 1.0) < 0.01
    assert np.isclose(saturation_exponent(1.0, 1.0), 0.5)


def test_power_law_parameter_recovery():
    """Parameter recovery only: data generated by a power law return its exponent.
    (Review point 11: this does NOT test any mechanism; an earlier version of this test claimed it did.)"""
    c = np.logspace(0, 3, 20)
    fit = fit_contact_law(c, 0.002 + 0.02 * c ** 0.5)
    assert abs(fit["power_law"]["k"] - 0.5) < 0.05


def test_hub_model_is_preferred_when_it_generated_the_data():
    """Model comparison works in the direction it can: noisy data from a fixed-size hub favour the hub form."""
    from poslayers.laws import _hub
    rng = np.random.default_rng(3); c = np.logspace(0, 3, 20); se = 0.003 * np.ones_like(c)
    wins = [fit_contact_law(c, _hub(c, 0.002, 0.15, 60.0) + rng.normal(0, se), sigma=se)
            for _ in range(30)]
    assert np.mean([w["hub"]["chi2"] < w["power_law"]["chi2"] for w in wins]) > 0.9


def test_power_law_vs_hub_does_not_identify_saturation():
    """Characterisation test, kept so the limitation stays visible.

    Data generated by a heterogeneous saturating response (each pair's occupancy c / (c + K), K log-normal) are ALSO fitted
    better by the fixed-hub form than by a power law. So 'the power law beats the hub' in real data argues against a single
    fixed-size hub but does not support saturation; the evidence for saturation must come from elsewhere (the expression
    tertile gradient)."""
    rng = np.random.default_rng(4); c = np.logspace(0, 3, 20); se = 0.003 * np.ones_like(c)
    def sat(): K = np.exp(rng.normal(np.log(30), 1.5, 4000)); return np.array([np.mean(0.15 * ci / (ci + K)) for ci in c])
    wins = [fit_contact_law(c, sat() + rng.normal(0, se), sigma=se) for _ in range(20)]
    assert np.mean([w["hub"]["chi2"] < w["power_law"]["chi2"] for w in wins]) > 0.5
POSLAYERS_FILE_067_END
echo "  updated  tests/test_laws.py"
chmod +x run_pipeline.sh

echo; echo "Checking..."
PY=""
for c in python3.13 python3.12 python3.11 python3.10 python3 python; do
  if command -v $c >/dev/null 2>&1 && $c -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" 2>/dev/null; then PY=$c; break; fi
done
left=$(grep -rln "/home/claude\|/mnt/user-data" --include=*.py --include=*.sh scripts figures src 2>/dev/null || true)
[ -z "$left" ] && echo "  ok  no absolute paths" || { echo "  FAILED: absolute paths in: $left"; exit 1; }
if [ -z "$PY" ]; then
  echo "  skipped tests: no Python >= 3.10 found. GitHub Actions will run them on push."
else
  if $PY -c "import pyflakes" 2>/dev/null; then
    out=$($PY -m pyflakes scripts/*.py figures/code/*.py src/poslayers/*.py | grep -E "undefined name|may be undefined" | grep -v "common" || true)
    [ -z "$out" ] && echo "  ok  no undefined names" || { echo "$out"; exit 1; }
  else echo "  skipped undefined-name check: pip install pyflakes"; fi
  if $PY -c "import pytest, numpy, scipy" 2>/dev/null; then
    PYTHONPATH="$PWD/src:${PYTHONPATH:-}" $PY -m pytest -q | tail -1
  else echo "  skipped tests: pip install -e \".[dev]\" in a Python >= 3.10 environment"; fi
fi
cat << 'MSG'

Done. Review with "git status" and "git diff --stat", then:

  git add -A
  git commit -m "v1.1.0: answer the code and methods review (see REVIEW_RESPONSE.md)"
  git push

GitHub Actions then runs the tests on Python 3.10 and 3.12 and the undefined-name check.
MSG
