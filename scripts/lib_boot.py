"""Genomic-block bootstrap for regression coefficients estimated on gene pairs.

Pairs that share genes or neighbourhoods are not independent. block_bootstrap() refits a linear model on resamples of
10-Mb genomic blocks (defined on the first gene of each pair) and returns a percentile interval and a bootstrap SE for
the requested coefficients. The design matrix is built once with patsy, so categorical levels stay fixed across replicates.
"""
import numpy as np, pandas as pd, pyannotables as pa
from patsy import dmatrices
from scipy import stats

_G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; _G = _G[~_G.index.duplicated()]


def blocks_for(g1, size=1e7):
    return (_G.Chromosome.reindex(g1).astype(str).values + ':' + (_G.Start.reindex(g1).fillna(0).values // size).astype(int).astype(str))


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
