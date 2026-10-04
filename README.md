# Positional Layers

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
pip install -e .
```

Python 3.10 or newer. Runtime dependencies are numpy, scipy, pandas, statsmodels and matplotlib;
`hicstraw` is needed only for the Hi-C extraction.

## The laws in thirty seconds

```python
from poslayers import simulate_genome, periodogram_identity, isochore_law, eqtl_law, saturation_exponent

sim = simulate_genome(decay_len=5, gc_bias_sd=0.25, seed=1)

# 1. The decomposition is an identity, not an approximation
periodogram_identity(sim["X"])["max_relative_error"]      # ~1e-15

# 2. Isochore law: how much correlation a per-sample GC bias creates, with no fitted parameter
isochore_law(var_b=0.06, gc_i=1.2, gc_j=0.9, sd_i=1.0, sd_j=1.0)

# 3. eQTL law: how much a shared causal variant couples two genes, sign included
eqtl_law(freq=0.3, beta1=0.4, beta2=-0.5)                 # negative: opposite effects

# 4. Saturation law: the contact exponent is one minus regulatory occupancy
saturation_exponent(contact=1.0, K=1.0)                   # 0.5
```

## Layout

| Path | What it holds |
| --- | --- |
| `src/poslayers/` | Reference implementation: `decompose.py`, `laws.py`, `simulate.py` |
| `scripts/` | Analysis scripts, one per result: atlas, eQTL, Hi-C, tumours, perturbations |
| `figures/code/` | One script per main and Extended Data figure |
| `docs/` | The interactive simulator served by GitHub Pages |
| `tests/` | Tests of the identity and of each law |
| `data/` | Empty. Inputs are public; see below |

## Data

Every input is public and none is redistributed here.

| Source | Where |
| --- | --- |
| GTEx v10/v11 counts, attributes, cis-eQTLs, SuSiE sets | https://gtexportal.org |
| Liver biopsy cohorts | GEO GSE130970, GSE135251, GSE162694 |
| GDC Pan-Cancer expression, mutations, copy number | https://xenabrowser.net |
| BeatAML2 | https://biodev.github.io/BeatAML2 |
| TAD stability landscapes | https://github.com/emcarthur/TAD-stability-heritability |
| In situ Hi-C (GM12878, IMR-90) | GEO GSE63525 |
| Edited erythroblasts | GEO GSE264491 |
| SLAM-seq after BET inhibition | GEO GSE100708 |
| Genome-wide Perturb-seq | https://gwps.wi.mit.edu |

Put the downloads under `data/` and point `POSLAYERS_DATA` at it:

```bash
export POSLAYERS_DATA=$PWD/data
python scripts/atlas.py           # the 36-tissue atlas
python scripts/predict.py         # the isochore law
python scripts/hic_test.py        # the contact exponent
python scripts/cont_tests.py      # STAG2 and the dosage term
```

## Reproducing the figures

```bash
python figures/code/fig1.py       # ... through fig6.py
python figures/code/ed4.py        # heatmap summary
```

Figures are written to `figures/output/`. The shared style, colour tokens and schematic
primitives live in `figures/code/style.py` and `figures/code/schem.py`.

See `scripts/SCRIPTS.md` for which script produces which figure panel.

## Tests

```bash
pytest -q
```

The tests check that the decomposition identity holds to machine precision, that the isochore law
recovers a simulated bias, that the eQTL law carries the right sign, and that the saturation
exponent is one minus occupancy.

## Citation

Pérez-Calixto, D. *et al.* Shared regulation of neighbouring genes underlies the positional
transcriptome (2026).

## Licence

MIT for the code. The public datasets keep the licences of their original sources.
