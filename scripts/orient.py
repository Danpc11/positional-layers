"""Cis layer by orientation and intergenic distance of adjacent gene pairs, per GTEx tissue (after GC and technical correction)."""
import os, sys
from poslayers.config import DATA, OUTDIR, upsert_csv, RESUME
import sys, os, numpy as np, pandas as pd, pyannotables as pa
OUT = OUTDIR + 'orientation_by_tissue.csv'; PAIRS = OUTDIR + 'pair_correlations.parquet'
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start', 'End', 'Strand']]
G38.columns = ['chr', 'start', 'end', 'strand']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc', 'type']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
BINS = [-np.inf, 0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]; LAB = ['overlap', '0-1kb', '1-5kb', '5-20kb', '20-100kb', '100-500kb', '>500kb']
done = set(pd.read_csv(OUT).tissue) if (RESUME and os.path.exists(OUT)) else set(); allpairs = []
for f in sys.argv[1:]:
    t = os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', '')
    if t in done or t == 'kidney_medulla': continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
    C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
    g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); D = resid_genes(D, np.column_stack([gcz, gcz ** 2]))
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    if tech.shape[1] < D.shape[1] - 10: D = resid_samples(D, tech)
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    same = g.chr.values[:-1] == g.chr.values[1:]; i = np.where(same)[0]
    P = pd.DataFrame({'g1': g.index[i], 'g2': g.index[i + 1], 'r': (Z[i] * Z[i + 1]).mean(1), 'dist': g.start.values[i + 1] - g.end.values[i],
                      's1': g.strand.values[i], 's2': g.strand.values[i + 1], 'pc': (g.type.values[i] == 'protein_coding') & (g.type.values[i + 1] == 'protein_coding')})
    P['orientation'] = np.where(P.s1 == P.s2, 'tandem', np.where((P.s1 == '-') & (P.s2 == '+'), 'divergent', 'convergent'))
    P['dist_bin'] = pd.cut(P.dist, BINS, labels=LAB)
    rng = np.random.default_rng(1); j = rng.permutation(len(Z)); floor = float(np.mean((Z[j[:-1]] * Z[j[1:]]).mean(1)))
    S = P.groupby(['orientation', 'dist_bin'], observed=True).r.agg(['size', 'mean']).reset_index().rename(columns={'size': 'n_pairs', 'mean': 'mean_r'})
    S['tissue'] = t; S['random_pair_floor'] = floor
    upsert_csv(S, OUT, ['tissue', 'orientation', 'dist_bin'])
    P['tissue'] = t; P[['tissue', 'g1', 'g2', 'r', 'dist', 'orientation', 'pc']].to_csv(OUTDIR + f'pairs_{t}.csv.gz', index=False)
    print(t, len(P), 'pairs; floor', round(floor, 4), flush=True)
