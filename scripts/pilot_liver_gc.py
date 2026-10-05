import os, sys
from poslayers.config import DATA, OUTDIR
import sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.environ.get('LIVER_SPECTRA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Liver_Spectra')), 'scripts'))  # liver pipeline: github.com/Danpc11/Liver_Spectra
from common import *
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False)
BM = BM.rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene start (bp)': 'start', 'Gene end (bp)': 'end', 'Chromosome/scaffold name': 'chr'}).drop_duplicates('gid').set_index('gid')
X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl')); comp = pd.read_pickle(inter('comp.pkl'))
s = M.index[M.estadio != 'Control']
k = keep.sort_values(['chr', 'grid_index']).set_index('gene_id').join(BM[['gc', 'start', 'end']], how='inner'); k = k[k.index.isin(A.index)]
k['length'] = (k.end - k.start).clip(lower=200); k['mid'] = (k.start + k.end) / 2
dens = []
for c, g in k.groupby('chr', sort=False):
    m = g.mid.values; dens += list(np.array([np.sum(np.abs(m - x) <= 1e6) for x in m]))
k['density'] = np.array(dens, float)
genes = k.index; chrs = k.chr.astype(str).values
print(f'liver biopsies: {len(genes)} genes with GC (of {len(A)})')
Y = A.loc[genes, s].values; Yc = Y - Y.mean(0, keepdims=True); mu = Yc.mean(1); D0 = Yc - mu[:, None]
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
zs = lambda a: (a - a.mean(0)) / a.std(0)
gc = k.gc.values; L = np.log(k.length.values)
D1 = resid_samples(D0, comp.loc[s].values)
lags = [1, 2, 3, 5, 10, 15, 20, 30, 150, 200, 300]
def prof(D, perm=False):
    out = {l: [] for l in lags}
    for c in CHR:
        idx = np.where(chrs == c)[0]
        if perm: idx = np.random.default_rng(3).permutation(idx)
        Z = D[idx]; Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
        for l in lags:
            if len(Z) > l + 5: out[l] += list((Z[:-l] * Z[l:]).mean(1))
    return np.array([np.mean(out[l]) for l in lags])
def summ(name, D):
    ex = prof(D) - prof(D, perm=True); return {'adjustment': name, 'cis_L1': ex[0], 'L2': ex[1], 'domain_L10_30': ex[4:8].mean(), 'long_L150_300': ex[8:].mean()}
rows = [summ('- composition', D1),
        summ('- composition - GC', resid_genes(D1, zs(np.column_stack([gc, gc ** 2])))),
        summ('- composition - length, expression, density', resid_genes(D1, zs(np.column_stack([L, L ** 2, mu, mu ** 2, k.density, k.density ** 2])))),
        summ('- composition - all gene properties incl. GC', resid_genes(D1, zs(np.column_stack([gc, gc ** 2, L, L ** 2, mu, mu ** 2, k.density, k.density ** 2]))))]
R = pd.DataFrame(rows); print(R.round(4).to_string(index=False))
# GC is itself clustered along the genome (isochores): neighbour correlation of gene GC
gz = []; 
for c in CHR:
    v = k.gc.values[chrs == c]; v = (v - v.mean()) / v.std(); gz.append([np.mean(v[:-l] * v[l:]) for l in (1, 10, 30)])
print('autocorrelation of gene GC content along the genome at 1, 10, 30 genes:', np.round(np.mean(gz, 0), 3))
# per-biopsy GC slope: how strongly each sample's deviations track gene GC
slope = np.array([np.polyfit(zs(gc), D1[:, i], 1)[0] for i in range(D1.shape[1])])
print(f'per-biopsy GC slope: SD {slope.std():.3f}; differs by cohort (ANOVA F = {__import__("scipy").stats.f_oneway(*[slope[M.loc[s, "cohorte"].values == c] for c in M.loc[s].cohorte.unique()]).statistic:.1f})')
R.to_csv(OUTDIR + 'p5a_liver_gc.csv', index=False)
