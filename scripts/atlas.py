"""Positional-layer atlas across GTEx tissues: landscape, technical isochore (GC) and cis covariance."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR, upsert_csv
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
    upsert_csv(pd.DataFrame([row]), OUT, ['tissue']); print(t, D.shape, 'ok', flush=True)
