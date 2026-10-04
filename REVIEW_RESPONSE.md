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
