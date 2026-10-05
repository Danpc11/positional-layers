"""Friction 2: is the elevated clustering of active genes within TADs in tumours (Liang et al. 2026) a copy-number effect?
Clustering z per sample: among the top-3,000 expressed protein-coding genes, number of TADs with >= 2 active genes versus
200 random gene sets of the same size (as in their bootstrap test). TADs: segments between conserved boundaries
(stability percentile >= 0.5, McArthur & Capra), genes placed with GRCh37 coordinates."""
import sys
from poslayers.config import DATA, OUTDIR, replace_groups_csv
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
COHORTS = sys.argv[1].split(',') if len(sys.argv) > 1 else ['BLCA', 'UCEC', 'STAD', 'COAD']
from poslayers.config import upsert_csv
SEED = {'BLCA': 11, 'UCEC': 12, 'STAD': 13, 'COAD': 14}
for c in COHORTS:
    rng = np.random.default_rng(SEED[c])   # per-cohort seed: results do not depend on which cohorts are run together
    z_ = np.load(OUTDIR + f'expr_{c}.npz', allow_pickle=True); X = z_['X']; genes = z_['genes']; samp = z_['samples']
    keep = pd.Index(genes).isin(tad.index[tad >= 0]); X = X[keep]; genes = genes[keep]; TI = tad.reindex(genes).values.astype(int)
    Y = np.log2(np.clip(2 ** X - 1, 0, None) / np.clip(2 ** X - 1, 0, None).sum(0) * 1e6 + 1)
    typ = np.array([s[13:15] for s in samp]); T = typ == '01'; N = typ == '11'
    CNz = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CN = pd.DataFrame(CNz['CN'], index=CNz['genes'], columns=CNz['samples'])
    tum = np.where(T & pd.Index(samp).isin(CN.columns))[0]; nor = np.where(N)[0]
    CNt = CN.reindex(index=genes, columns=samp[tum]).fillna(0).values; burden = np.mean(np.abs(CNt) > 0.2, 0)
    zN = cluster_z(Y[:, nor]); zT = cluster_z(Y[:, tum]); M = np.abs(CNt) <= 0.2; zT_clean = cluster_z(Y[:, tum], ok_mask=M)
    zN_masked = cluster_z(Y[:, nor], ok_mask=M[:, rng.choice(M.shape[1], len(nor), replace=True)])   # same gene exclusion as random tumours
    per = pd.concat([pd.DataFrame({'cohort': c, 'sample': samp[nor], 'group': 'normal', 'z': zN}), pd.DataFrame({'cohort': c, 'sample': samp[tum], 'group': 'tumour', 'z': zT}),
                     pd.DataFrame({'cohort': c, 'sample': samp[tum], 'group': 'tumour_CNA_removed', 'z': zT_clean}), pd.DataFrame({'cohort': c, 'sample': samp[nor], 'group': 'normal_matched_exclusion', 'z': zN_masked})])
    replace_groups_csv(per, OUTDIR + 'friction2_per_sample_z.csv', ['cohort'], [c], key=['cohort', 'sample', 'group'])
    rows.append({'cohort': c, 'normals': len(nor), 'tumours': len(tum), 'z_normal': np.median(zN), 'z_tumour': np.median(zT), 'P_tumour_vs_normal': stats.mannwhitneyu(zT, zN).pvalue,
                 'rho_z_vs_CNA_burden': stats.spearmanr(zT, burden)[0], 'P_rho': stats.spearmanr(zT, burden)[1],
                 'z_tumour_CNA_genes_removed': np.nanmedian(zT_clean), 'z_normal_same_exclusion': np.nanmedian(zN_masked), 'P_clean_vs_normal_matched': stats.mannwhitneyu(zT_clean[np.isfinite(zT_clean)], zN_masked[np.isfinite(zN_masked)]).pvalue})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); upsert_csv(R, OUTDIR + 'friction2_tumour_tad_clustering.csv', ['cohort']); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
