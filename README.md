# Positional layers

Code for **"Common cis-regulatory inputs shape local gene co-expression across human tissues"**.

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
| `figures/code/` | One script per main and Extended Data figure; `style.py` places panel letters and keeps legends off the data |
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

## Citation

Pérez-Calixto, D. *et al.* Shared cis-regulation drives local gene co-expression across human tissues (2026).

## Licence

MIT for the code. The public datasets keep the licences of their original sources.
