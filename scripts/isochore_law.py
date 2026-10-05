"""Isochore law: in-sample and out-of-sample tests of the GC component of covariance.

For each tissue and lag L we report three things, all on the SAME normalisation (raw per-gene SD):
  old     : corr(raw) - corr(quadratic-GC-corrected)      [difference of correlations, for comparison]
  insample: [cov(raw) - cov(raw minus per-sample LINEAR GC fit)] / (sd_i sd_j)  versus  Var(b) g_i g_j / (sd_i sd_j)
            -> nearly an identity; the gap is the cross term cov(b g, residual), i.e. how far b is from independent of biology
  heldout : b_s estimated on odd chromosomes only; prediction Var(b_odd) g_i g_j / (sd_i sd_j) for gene pairs on EVEN chromosomes,
            observation = GC covariance on even chromosomes from their own fit. This is the non-trivial test: a single
            per-sample scalar measured elsewhere in the genome must predict the positional covariance in held-out chromosomes.
Also: corr(b_odd, b_even) across samples, and the share of the quadratic correction carried by the linear term.
"""
import os
from poslayers.config import DATA, OUTDIR, upsert_csv, RESUME
import os
import glob
import numpy as np
import pandas as pd
import pyannotables as pa
CHR = [str(i) for i in range(1, 23)]
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()][['Chromosome', 'Start']]; G.columns = ['chr', 'start']; G['chr'] = G.chr.astype(str)
G = G[G.chr.isin(CHR)].join(BM[['gc']], how='inner')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN']).set_index('SAMPID')
LAGS = [1, 2, 5, 10, 20, 30]
OUT = OUTDIR + 'isochore_v2.csv'
done = set(pd.read_csv(OUT).tissue) if (RESUME and os.path.exists(OUT)) else set()

def lag_pairs(chrs, keep, L):
    """index pairs (i, i+L) within the same chromosome, restricted to genes where keep is True"""
    I = []
    for c in np.unique(chrs[keep]):
        idx = np.where((chrs == c) & keep)[0]
        if len(idx) > L + 5: I.append(np.column_stack([idx[:-L], idx[L:]]))
    return np.vstack(I) if I else np.empty((0, 2), int)

def fit_slopes(D, g, rows):
    """per-sample OLS of D[rows, s] on [1, g[rows]]; returns intercepts a_s and slopes b_s"""
    X = np.column_stack([np.ones(rows.sum()), g[rows]]); B = np.linalg.lstsq(X, D[rows], rcond=None)[0]
    return B[0], B[1]

for f in sorted(glob.glob(DATA + 'gtex/*_gct.gz')):
    t = os.path.basename(f).replace('gene_reads_adult_gtex_v11_', '').replace('gene_reads_v10_', '').replace('_gct.gz', '').lower().replace('-', '_')
    if t == 'kidney_medulla' or t in done: continue
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
    C = C.loc[C.index.intersection(G.index)]; C = C[C.median(axis=1) >= 10]
    g = G.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)   # genes x samples
    gz = (g.gc.values - g.gc.mean()) / g.gc.std(); S = D.shape[1]
    sd = D.std(1)                                                                  # common normalisation for every quantity
    allg = np.ones(len(gz), bool); odd = np.isin(chrs, [str(i) for i in range(1, 23, 2)]); even = ~odd
    a_all, b_all = fit_slopes(D, gz, allg); a_o, b_o = fit_slopes(D, gz, odd); a_e, b_e = fit_slopes(D, gz, even)
    R_lin = D - (a_all[None, :] + np.outer(gz, b_all))                             # linear GC removed
    Xq = np.column_stack([np.ones(len(gz)), gz, gz ** 2 - (gz ** 2).mean()]); Bq = np.linalg.lstsq(Xq, D, rcond=None)[0]; R_q = D - Xq @ Bq
    R_e = D - (a_e[None, :] + np.outer(gz, b_e))                                   # even chromosomes, their own fit
    cov = lambda M, P: np.mean((M[P[:, 0]] - M[P[:, 0]].mean(1, keepdims=True)) * (M[P[:, 1]] - M[P[:, 1]].mean(1, keepdims=True)), 1)
    def corr(M, P):
        s = M.std(1) + 1e-12; return cov(M, P) / (s[P[:, 0]] * s[P[:, 1]])
    row = {'tissue': t, 'samples': S, 'genes': len(gz), 'var_b': b_all.var(), 'var_b_odd': b_o.var(), 'var_b_even': b_e.var(),
           'corr_b_odd_even': np.corrcoef(b_o, b_e)[0, 1]}
    for L in LAGS:
        P = lag_pairs(chrs, allg, L); nn = sd[P[:, 0]] * sd[P[:, 1]]
        row[f'old_obs_L{L}'] = np.mean(corr(D, P) - corr(R_q, P))
        row[f'old_pred_L{L}'] = np.mean(b_all.var() * gz[P[:, 0]] * gz[P[:, 1]] / nn)
        row[f'in_obs_L{L}'] = np.mean((cov(D, P) - cov(R_lin, P)) / nn)
        row[f'in_pred_L{L}'] = row[f'old_pred_L{L}']
        row[f'quad_obs_L{L}'] = np.mean((cov(D, P) - cov(R_q, P)) / nn)               # full basis, same normalisation
        Pe = lag_pairs(chrs, even, L); ne = sd[Pe[:, 0]] * sd[Pe[:, 1]]
        row[f'held_pred_L{L}'] = np.mean(b_o.var() * gz[Pe[:, 0]] * gz[Pe[:, 1]] / ne)
        row[f'held_obs_L{L}'] = np.mean((cov(D, Pe) - cov(R_e, Pe)) / ne)
    upsert_csv(pd.DataFrame([row]), OUT, ['tissue'])
    print(t, S, f"corr(b_odd,b_even)={row['corr_b_odd_even']:.3f}  L10 old {row['old_obs_L10']:.4f}/{row['old_pred_L10']:.4f}  held {row['held_obs_L10']:.4f}/{row['held_pred_L10']:.4f}", flush=True)
