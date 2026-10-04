# Which script produces what

Every script reads its inputs from `$POSLAYERS_DATA` (default `data/`). See the data table in
the top-level README for where to download each input.

| Script | Produces | Figure |
| --- | --- | --- |
| `atlas.py` | 36-tissue atlas: landscape share, GC batch η², cis and domain coupling, decay length | Fig. 1c,e,f · Fig. 2b · ED Fig. 4a |
| `theory_sim.py` | Simulated genomes; verifies the decomposition identity and the decay-length estimators | Fig. 1b,h |
| `predict.py` | Isochore law, predicted versus observed GC component | Fig. 1g · ED Fig. 4b |
| `predict_cal.py` | Calibration of the isochore law across distances | ED Fig. 4b |
| `orient.py` | Coupling by intergenic distance and orientation | Fig. 2c |
| `lorentz.py` | Covariance spectrum; one versus two Lorentzians versus power law | Fig. 2e,f · ED Fig. 1b |
| `eqtl_share.py` | eQTL-sharing classes per tissue | Fig. 3b |
| `eqtl_test.py` | Adjusted effect of sharing a same- or opposite-direction variant | Fig. 3c |
| `coloc.py` | SuSiE colocalisation probability per adjacent pair | Fig. 3d |
| `coloc_test.py` | Dose-response of coupling on colocalisation; tissue specificity | Fig. 3d,f |
| `tad_test.py` | Same-TAD increase in coupling at equal distance | Fig. 4b |
| `hic_test.py` | Hi-C contact versus coupling in GM12878; the contact exponent | Fig. 4c,d |
| `hic_test_imr90.py` | The same in IMR-90 | Fig. 4c,d |
| `derive.py` | Hub model versus power law; rejects the fixed-size hub | Supplementary Note 1.5 |
| `derive2.py` | Contact exponent by expression tertile; the saturation law | Fig. 4f |
| `tad_tumour.py` | Cis excess within versus across conserved boundaries in tumours | Fig. 4e |
| `extract.py` | Extracts the TCGA cohorts from the Pan-Cancer release | — |
| `cn_continuous.py` | Per-gene continuous copy number from the segment file | — |
| `cont_tests.py` | STAG2 and CTCF effects on cis excess, copy-number corrected | Fig. 4g |
| `cohesin_v2.py` | The same using GISTIC calls (sensitivity analysis) | Fig. 4g |
| `aneuploidy.py` | Long-range floor versus copy-number burden | ED Fig. 2a,b |
| `beat_test.py` | STAG2 replication in BeatAML2 | Fig. 4g |
| `replicate.py` | Inverse-variance combination of the bladder and AML estimates | Fig. 4g |
| `edit_cis.py` | Neighbour responses to the HBG1/2 and BCL11A edits | Fig. 5c,d |
| `predict_first.py` | The blind neighbour predictions, made before opening the edited-cell data | Fig. 5c |
| `test_edits.py` | Scores those predictions | Fig. 5c |
| `drug_test.py` | Propagation to coupled neighbours in SLAM-seq | Fig. 5f |
| `crispri_test.py` | Neighbour response to CRISPRi by distance and coupling | Fig. 5e |
| `crispri_analyse.py` | The coupling × knockdown interaction test | Fig. 5e |
| `fric1.py` | TAD effect across distance, with and without contact | ED Fig. 3 |
| `fric2.py` | TAD clustering of active genes in tumours, with and without copy-number correction | ED Fig. 2c |
| `simulate_demo.py` | The web simulator, in Python | — |

## Superseded

`cohesin_test.py` was the first pass at the cohesin analysis and is kept only for the record;
`cont_tests.py` supersedes it and is what the paper reports.
