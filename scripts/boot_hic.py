"""Donor x genomic-block bootstrap for the Hi-C results (Fig. 4c,d,f).
Each replicate resamples GTEx donors with replacement (re-standardising expression, then recomputing every pair's coupling) AND
resamples 10-Mb genomic blocks with replacement (blocks defined on the first gene's position). Statistics per replicate:
  contact_coef : coefficient of log2 O/E contact in r ~ B-spline(log10 distance, 5 df) + log_oe
  k_distance   : log-log slope of mean coupling on mean contact across the six distance bins (as in Fig. 4d)
  k_contact    : power-law exponent fitted to mean coupling over 20 contact quantiles (as in Supplementary Note 1)
  k_t1..k_t3   : the same exponent within expression tertiles (Fig. 4f)
Resumable: appends one row per replicate to boot_hic_<cell>.csv."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, pyannotables as pa
from scipy.optimize import curve_fit
from patsy import dmatrix
cell, tissue, B = sys.argv[1], sys.argv[2], int(sys.argv[3])
OUT = OUTDIR + f'boot_hic_{cell}.csv'; done = len(pd.read_csv(OUT)) if os.path.exists(OUT) else 0
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
C = pd.read_csv(DATA + f'gtex/gene_reads_adult_gtex_v11_{tissue}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description')
C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]; samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); expr = Y.mean(1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[C.index, 'gc'].values; gz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gz)), gz, gz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T                  # technical residuals, computed once on all donors
pos = {g: i for i, g in enumerate(C.index)}
H = pd.read_csv(OUTDIR + f'hic_coupling_{cell}.csv.gz', dtype={'chr': str}); H = H[H.g1.isin(pos) & H.g2.isin(pos) & (H.contact_KR > 0)].reset_index(drop=True)
i1 = H.g1.map(pos).values; i2 = H.g2.map(pos).values
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
start = G.Start.reindex(H.g1).values; block = (H.chr.astype(str) + ':' + (np.nan_to_num(start) // 1e7).astype(int).astype(str)).values
ub, binv = np.unique(block, return_inverse=True); members = [np.where(binv == k)[0] for k in range(len(ub))]
X0 = np.asarray(dmatrix('bs(x, df=5)', {'x': H.log_d.values}, return_type='dataframe')); Xf = np.column_stack([X0, H.log_oe.values])
dbin = pd.cut(H.tss_distance, [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6], labels=False).values
emin = np.minimum(expr[i1], expr[i2]); tert = pd.qcut(emin, 3, labels=False)
pw = lambda c, r0, A, k: r0 + A * c ** k
def k_fit(c, r, nq):
    q = pd.qcut(np.log10(c), nq, labels=False, duplicates='drop'); g = pd.DataFrame({'c': c, 'r': r, 'q': q}).groupby('q').agg(c=('c', 'median'), r=('r', 'mean'), s=('r', 'sem'))
    try: return curve_fit(pw, g.c, g.r, p0=[0.0, 0.002, 0.5], sigma=g.s, bounds=([-0.05, 0, 0.01], [0.05, 5, 3]), maxfev=20000)[0][2]
    except Exception: return np.nan
def stats_for(r, idx):
    out = {}
    beta = np.linalg.lstsq(Xf[idx], r[idx], rcond=None)[0]; out['contact_coef'] = beta[-1]
    dd = pd.DataFrame({'b': dbin[idx], 'r': r[idx], 'c': H.contact_KR.values[idx], 'd': H.tss_distance.values[idx]}).groupby('b').agg(r=('r', 'mean'), c=('c', 'mean'))
    out['k_distance'] = np.polyfit(np.log10(dd.c), np.log10(dd.r.clip(lower=1e-4)), 1)[0]
    out['k_contact'] = k_fit(H.contact_KR.values[idx], r[idx], 20)
    for t in range(3):
        j = idx[tert[idx] == t]; out[f'k_t{t + 1}'] = k_fit(H.contact_KR.values[j], r[j], 15)
    return out
def coupling(cols):
    M = D[:, cols]; Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9); r = np.empty(len(H))
    for s in range(0, len(H), 150000): r[s:s + 150000] = (Z[i1[s:s + 150000]] * Z[i2[s:s + 150000]]).mean(1)
    return r
rng = np.random.default_rng(2024 + done)
if done == 0:
    est = stats_for(coupling(np.arange(D.shape[1])), np.arange(len(H))); est['rep'] = -1
    pd.DataFrame([est]).to_csv(OUT, index=False); print('point estimate', {k: round(v, 4) for k, v in est.items()}, flush=True); done = 1
for b in range(done, B + 1):
    cols = rng.integers(0, D.shape[1], D.shape[1]); idx = np.concatenate([members[k] for k in rng.integers(0, len(members), len(members))])
    st = stats_for(coupling(cols), idx); st['rep'] = b
    pd.DataFrame([st]).to_csv(OUT, mode='a', header=False, index=False)
    if b % 10 == 0: print(b, {k: round(v, 3) for k, v in st.items()}, flush=True)
