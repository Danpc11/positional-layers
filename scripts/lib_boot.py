"""Genomic-block bootstrap for regression coefficients estimated on gene pairs.

Pairs that share genes or neighbourhoods are not independent. block_bootstrap() refits a linear model on resamples of
10-Mb genomic blocks (defined on the first gene of each pair) and returns a percentile interval and a bootstrap SE for
the requested coefficients. The design matrix is built once with patsy, so categorical levels stay fixed across replicates.

Blocks are assigned from the GRCh38 start of the first gene of each pair (Ensembl 100). A gene without a GRCh38 position
gets a block of its own, so it is resampled independently rather than merged into another block. The block size is
10 Mb by default and can be changed with POSLAYERS_BLOCK_MB for sensitivity analyses. p_block is a normal approximation
(estimate / bootstrap SE), not an empirical bootstrap or permutation P value under the null.
"""
import os, warnings
import numpy as np
import pyannotables as pa
from patsy import dmatrices
from scipy import stats

_G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; _G = _G[~_G.index.duplicated()]


BLOCK = float(os.environ.get('POSLAYERS_BLOCK_MB', '10')) * 1e6


def blocks_for(g1, size=None):
    size = BLOCK if size is None else size
    g1 = np.asarray(g1); start = _G.Start.reindex(g1).values; chrom = _G.Chromosome.reindex(g1).astype(str).values
    out = np.array([f'{c}:{int(s // size)}' if np.isfinite(s) else f'unplaced:{g}' for c, s, g in zip(chrom, start, g1)], dtype=object)
    n = int(np.sum(~np.isfinite(start)))
    if n: warnings.warn(f'{n} genes without a GRCh38 position were given their own block')
    return out


def block_bootstrap(formula, data, terms, B=300, seed=0):
    y, X = dmatrices(formula, data, return_type='dataframe')
    blk = blocks_for(data.loc[X.index, 'g1'].values); groups = [np.where(blk == b)[0] for b in np.unique(blk)]
    cols = [list(X.columns).index(t) for t in terms]; Xv, yv = X.values, y.values.ravel()
    est = np.linalg.lstsq(Xv, yv, rcond=None)[0][cols]
    rng = np.random.default_rng(seed); draws = []
    for _ in range(B):
        ix = np.concatenate([groups[k] for k in rng.integers(0, len(groups), len(groups))])
        draws.append(np.linalg.lstsq(Xv[ix], yv[ix], rcond=None)[0][cols])
    draws = np.array(draws); se = draws.std(0)
    return {t: {'est': est[i], 'ci_low': np.percentile(draws[:, i], 2.5), 'ci_high': np.percentile(draws[:, i], 97.5),
                'boot_se': se[i], 'p_block': 2 * stats.norm.sf(abs(est[i] / se[i])) if se[i] > 0 else np.nan, 'blocks': len(groups)}
            for i, t in enumerate(terms)}
