"""Quantitative predictions per GTEx tissue.
(1) Isochore law: GC-induced correlation between genes i,j = Var(b) * gc_i * gc_j / (sd_i * sd_j), b = per-sample GC slope.
(2) eQTL law: correlation induced by a shared causal variant = sum over shared credible sets of
    P(same variant) * 2p(1-p) * (afc_1/2) * (afc_2/2) / (sd_1 * sd_2)   [afc: log2 allelic fold change; afc/2 = per-allele log2 effect]."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, glob, numpy as np, pandas as pd, pyannotables as pa, pyarrow.parquet as pq
CHR = [str(i) for i in range(1, 23)] + ['X']; norm = lambda s: s.lower().replace('-', '_')
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start', 'End']]; G38.columns = ['chr', 'start', 'end']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
SUS = {norm(os.path.basename(f).replace('_v11_eQTLs_SuSiE_summary.parquet', '')): f for f in glob.glob(DATA + 'gtex_eqtl/*_v11_eQTLs_SuSiE_summary.parquet')}
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
LAGS = [1, 2, 5, 10, 20, 30]
for f in sys.argv[1:]:
    t = norm(os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', ''))
    if t == 'kidney_medulla' or os.path.exists(OUTDIR + f'pred_pairs_{t}.csv.gz'): continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
    g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std()
    # ---- (1) isochore law
    b = np.array([np.polyfit(gcz, D[:, i], 1)[0] for i in range(D.shape[1])]); varb = b.var(); sd_raw = D.std(1)
    Dg = resid_genes(D, np.column_stack([gcz, gcz ** 2]))
    def lagcorr(M, L):
        Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9); return np.mean(np.concatenate([np.mean(Z[i[:-L]] * Z[i[L:]], 1) for c in CHR for i in [np.where(chrs == c)[0]] if len(i) > L + 5]))
    law = {'tissue': t, 'var_GC_slope': varb}
    for L in LAGS:
        pr = np.concatenate([varb * gcz[i[:-L]] * gcz[i[L:]] / (sd_raw[i[:-L]] * sd_raw[i[L:]]) for c in CHR for i in [np.where(chrs == c)[0]] if len(i) > L + 5]).mean()
        law[f'pred_L{L}'] = pr; law[f'obs_gc_component_L{L}'] = lagcorr(D, L) - lagcorr(Dg, L)
    pd.DataFrame([law]).to_csv(OUTDIR + 'isochore_law.csv', mode='a', header=not os.path.exists(OUTDIR + 'isochore_law.csv'), index=False)
    # ---- (2) eQTL law, on GC- and technically corrected expression (as for the cis layer)
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    R = resid_samples(Dg, tech) if tech.shape[1] < Dg.shape[1] - 10 else Dg; sd = R.std(1); Z = (R - R.mean(1, keepdims=True)) / (sd[:, None] + 1e-9)
    same = chrs[:-1] == chrs[1:]; ii = np.where(same)[0]; pos = {gid: k for k, gid in enumerate(g.index)}
    if t in SUS:
        E = pq.read_table(SUS[t], columns=['phenotype_id', 'variant_id', 'pip', 'cs_id', 'af', 'afc']).to_pandas(); E['gid'] = E.phenotype_id.str.split('.').str[0]; E = E[E.gid.isin(pos)]
        CS = {gg: [(dict(zip(c.variant_id, c.pip)), dict(zip(c.variant_id, c.afc)), dict(zip(c.variant_id, c.af))) for _, c in d.groupby('cs_id')] for gg, d in E.groupby('gid')}
        rows = []
        for k in ii:
            g1, g2 = g.index[k], g.index[k + 1]; A, B = CS.get(g1), CS.get(g2)
            if A is None or B is None: continue
            cov = 0.0; pmax = 0.0
            for p1, f1, af1 in A:
                for p2, f2, af2 in B:
                    sh = [v for v in p1.keys() & p2.keys()]
                    if not sh: continue
                    ps = sum(p1[v] * p2[v] for v in sh); pmax = max(pmax, ps)
                    vs = [v for v in sorted(sh, key=lambda v: -p1[v] * p2[v]) if np.isfinite(f1.get(v, np.nan)) and np.isfinite(f2.get(v, np.nan))]
                    if not vs: continue
                    v = vs[0]; p = af1[v]; cov += ps * 2 * p * (1 - p) * (f1[v] / 2) * (f2[v] / 2)
            if pmax < 0.1: continue
            rows.append({'g1': g1, 'g2': g2, 'p_coloc': pmax, 'pred_r': cov / (sd[k] * sd[k + 1]), 'obs_r': float(np.mean(Z[k] * Z[k + 1])), 'dist': g.start.values[k + 1] - g.end.values[k]})
        P = pd.DataFrame(rows); P['tissue'] = t; P.to_csv(OUTDIR + f'pred_pairs_{t}.csv.gz', index=False)
        print(t, f'varb {varb:.4f} | pred/obs GC L20 {law["pred_L20"]:.4f}/{law["obs_gc_component_L20"]:.4f} | coloc pairs {len(P)} | r(pred,obs) {P[["pred_r", "obs_r"]].corr().iloc[0, 1]:.2f}', flush=True)
