import os, sys
from poslayers.config import OUTDIR
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
