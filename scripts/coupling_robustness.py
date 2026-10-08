"""Robustness of adjacent-gene coupling as a property of gene pairs (not only a Pearson average).

For each tissue (GTEx, GC- and technically-corrected expression, as in the atlas):
  1. Pearson and Spearman coupling of adjacent pairs, and their agreement across pairs;
  2. reproducibility: donors split at random into two halves; correlation between the coupling of each pair in half A and
     in half B, for adjacent pairs and for distant pairs (> 50 genes apart on the same chromosome);
  3. a matched null: distant pairs matched to each adjacent pair by deciles of mean expression and of GC content;
  4. pair-level tests: fraction of pairs with Benjamini-Hochberg q < 0.05 (Pearson test), adjacent versus matched distant.
Each analysis is run with 0 and with 15 expression principal components removed (genome-wide, non-positional co-expression).
Usage: python coupling_robustness.py TISSUE [TISSUE ...]
Outputs: OUTDIR/coupling_robustness.csv (rows replaced per tissue) and, with 15 PCs removed,
         OUTDIR/coupling_robustness_pairs_<tissue>.csv.gz (pair-level values; Extended Data Fig. 6)
"""
import sys
import numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
from poslayers.config import DATA, OUTDIR, replace_groups_csv

BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])]

def residuals(tis):
    C = pd.read_csv(DATA + f'gtex/gene_reads_adult_gtex_v11_{tis}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [s for s in C.columns if s in SA.index and pd.notna(SA.loc[s, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(BM.index).intersection(G.index)]; C = C[C.median(axis=1) >= 10]
    g = G.loc[C.index].sort_values(['Chromosome', 'Start']); C = C.loc[g.index]
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); mean_expr = Y.mean(1); D = Y - Y.mean(1, keepdims=True)
    gz = BM.loc[C.index, 'gc'].values; gz = (gz - gz.mean()) / gz.std()
    G1 = np.column_stack([np.ones(len(gz)), gz, gz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    if T.shape[1] < D.shape[1] - 10:                       # as in atlas.py: technical covariates only when samples clearly outnumber them
        D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T
    donors = np.array(['-'.join(s.split('-')[:2]) for s in samp])
    return D, g.Chromosome.astype(str).values, mean_expr, gz, donors

def pair_corr(D, i, j, rank=False):
    X = D[:, :]
    if rank: X = np.apply_along_axis(stats.rankdata, 1, X)
    Z = (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)
    return (Z[i] * Z[j]).mean(1)

def matched_distant(chrs, me, gz, ai, aj, rng):
    """Distant partners (> 50 genes apart, same chromosome) matched by deciles of mean expression and GC content."""
    key = pd.qcut(me, 10, labels=False) * 10 + pd.qcut(gz, 10, labels=False); pool = {}
    for k, c in enumerate(chrs): pool.setdefault((c, key[k]), []).append(k)
    di, dj = [], []
    for i, j in zip(ai, aj):
        for _ in range(20):
            ci, cj = pool[(chrs[i], key[i])], pool[(chrs[j], key[j])]; u, v = ci[rng.integers(len(ci))], cj[rng.integers(len(cj))]
            if abs(u - v) > 50: di.append(u); dj.append(v); break
    return np.array(di), np.array(dj)

def frac_q(r, n_s):
    t_ = r * np.sqrt((n_s - 2) / np.clip(1 - r ** 2, 1e-12, None)); pv = 2 * stats.t.sf(np.abs(t_), n_s - 2)
    o = np.argsort(pv); q = np.empty_like(pv); q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    return float(np.mean((q < 0.05) & (r > 0))), float(np.mean((q < 0.05) & (r < 0)))

for tis in sys.argv[1:]:
    D0, chrs, me, gz, donors = residuals(tis); out = []
    idx = np.arange(len(chrs)); same = chrs[:-1] == chrs[1:]; ai, aj = idx[:-1][same], idx[1:][same]
    di, dj = matched_distant(chrs, me, gz, ai, aj, np.random.default_rng(7))
    u = np.unique(donors); hA = np.isin(donors, np.random.default_rng(8).choice(u, len(u) // 2, replace=False))
    for K in (0, 15):
        D = D0.copy()
        if K:
            U, sv, Vt = np.linalg.svd(D - D.mean(1, keepdims=True), full_matrices=False); D = D - (U[:, :K] * sv[:K]) @ Vt[:K]
        n_s = D.shape[1]
        r_adj, r_adj_s = pair_corr(D, ai, aj), pair_corr(D, ai, aj, rank=True); r_dis, r_dis_s = pair_corr(D, di, dj), pair_corr(D, di, dj, rank=True)
        rep = lambda i, j: stats.pearsonr(pair_corr(D[:, hA], i, j), pair_corr(D[:, ~hA], i, j))[0]
        pa_pos, pa_neg = frac_q(r_adj, n_s); pd_pos, pd_neg = frac_q(r_dis, n_s)
        row = {'tissue': tis, 'pcs_removed': K, 'samples': n_s, 'adjacent_pairs': len(ai), 'matched_distant_pairs': len(di),
               'pearson_adjacent_median': np.median(r_adj), 'spearman_adjacent_median': np.median(r_adj_s), 'pearson_spearman_r_across_pairs': stats.pearsonr(r_adj, r_adj_s)[0],
               'pearson_distant_median': np.median(r_dis), 'spearman_distant_median': np.median(r_dis_s),
               'split_half_reproducibility_adjacent': rep(ai, aj), 'split_half_reproducibility_distant': rep(di, dj),
               'frac_adjacent_q05_positive': pa_pos, 'frac_adjacent_q05_negative': pa_neg, 'frac_distant_q05_positive': pd_pos, 'frac_distant_q05_negative': pd_neg}
        if K == 15:
            half = lambda i, j: (pair_corr(D[:, hA], i, j), pair_corr(D[:, ~hA], i, j))
            (aA, aB), (dA, dB) = half(ai, aj), half(di, dj)
            pd.concat([pd.DataFrame({'type': 'adjacent', 'pearson': r_adj, 'spearman': r_adj_s, 'half_a': aA, 'half_b': aB}),
                       pd.DataFrame({'type': 'matched distant', 'pearson': r_dis, 'spearman': r_dis_s, 'half_a': dA, 'half_b': dB})]).round(5) \
              .to_csv(OUTDIR + f'coupling_robustness_pairs_{tis}.csv.gz', index=False)
        out.append(row); print(tis, K, 'done', flush=True)
    replace_groups_csv(pd.DataFrame(out), OUTDIR + 'coupling_robustness.csv', ['tissue'], [tis], key=['tissue', 'pcs_removed'])
