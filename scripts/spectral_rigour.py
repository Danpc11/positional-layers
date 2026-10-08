"""Spectral checks requested in review: autocorrelation-preserving nulls, Whittle likelihood, and a base-pair decay scale.

For one GTEx tissue (same preprocessing as atlas.py):
  1. Peaks of the mean profile (landscape) of each chromosome tested against nulls that keep local autocorrelation:
     block permutation (blocks of 10, 30 and 100 genes, 200 replicates) and AR(1) red noise fitted to the profile.
     A frequency is called a peak when its power exceeds the 95th percentile of the null MAXIMUM over frequencies
     (family-wise error 5% per chromosome).
  2. Covariance spectra of the corrected residuals fitted by Whittle likelihood (one Lorentzian, two Lorentzians, power
     law, each with a constant floor), summed over chromosomes; AIC compared.
  3. Coupling against physical distance: mean correlation of gene pairs in 20 log-spaced bins from 5 kb to 5 Mb, minus
     the mean for pairs 20-40 Mb apart; a two-exponential model fitted on the base-pair axis.
Usage: python spectral_rigour.py TISSUE
Outputs: OUTDIR/spectral_rigour.csv (row replaced per tissue), OUTDIR/coupling_by_bp_<tissue>.csv and
         OUTDIR/landscape_null_<tissue>_chr1.csv (Fig. 2e and Extended Data Fig. 2f)
"""
import sys
import numpy as np, pandas as pd, pyannotables as pa
from scipy.optimize import minimize, curve_fit
from poslayers.config import DATA, OUTDIR, replace_groups_csv

BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])]
tis = sys.argv[1]
C = pd.read_csv(DATA + f'gtex/gene_reads_adult_gtex_v11_{tis}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
samp = [s for s in C.columns if s in SA.index and pd.notna(SA.loc[s, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(BM.index).intersection(G.index)]; C = C[C.median(axis=1) >= 10]
g = G.loc[C.index].sort_values(['Chromosome', 'Start']); C = C.loc[g.index]; chrs = g.Chromosome.astype(str).values
tss = np.where(g.Strand.astype(str).isin(['1', '+']), g.Start, g.End).astype(float)
X = np.log2(C.values / C.values.sum(0) * 1e6 + 1); mu = X.mean(1); D = X - mu[:, None]
gz = BM.loc[C.index, 'gc'].values; gz = (gz - gz.mean()) / gz.std(); G1 = np.column_stack([np.ones(len(gz)), gz, gz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
if T.shape[1] < D.shape[1] - 10:                           # as in atlas.py: technical covariates only when samples clearly outnumber them
    D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T
rng = np.random.default_rng(11); row = {'tissue': tis, 'samples': len(samp), 'genes': len(g)}
pw = lambda v: np.abs(np.fft.rfft(v - v.mean())[1:len(v) // 2]) ** 2 / len(v)
# ---- 1. peaks of the landscape against autocorrelation-preserving nulls
counts = {f'block{b}': 0 for b in (10, 30, 100)}; counts['ar1'] = 0; nfreq = 0
for c in np.unique(chrs):
    m = mu[chrs == c]; n = len(m)
    if n < 200: continue
    P = pw(m); nfreq += len(P)
    for b in (10, 30, 100):
        blocks = [m[i:i + b] for i in range(0, n, b)]; mx = []
        for _ in range(200):
            perm = np.concatenate([blocks[k] for k in rng.permutation(len(blocks))]); mx.append(pw(perm).max())
        counts[f'block{b}'] += int((P > np.percentile(mx, 95)).sum())
    z = m - m.mean(); phi = np.corrcoef(z[:-1], z[1:])[0, 1]; s2 = z.var() * (1 - phi ** 2)
    f = np.arange(1, len(P) + 1) / n; Sar = s2 / (1 + phi ** 2 - 2 * phi * np.cos(2 * np.pi * f))
    counts['ar1'] += int((P / Sar > -np.log(0.05 / len(P))).sum())                 # exponential ordinates, Bonferroni over frequencies
    if c == '1':                                                                    # spectrum and thresholds of one chromosome (Extended Data Fig. 2f)
        bmax = []; rng_f = np.random.default_rng(99)                               # own generator: leaves the main nulls unchanged
        for _ in range(200):
            blocks = [m[i:i + 30] for i in range(0, n, 30)]; bmax.append(pw(np.concatenate([blocks[k] for k in rng_f.permutation(len(blocks))])).max())
        pd.DataFrame({'tissue': tis, 'chromosome': c, 'frequency': f, 'landscape_power': P, 'red_noise_threshold': Sar * -np.log(0.05 / len(P)),
                      'block30_threshold': np.percentile(bmax, 95)}).to_csv(OUTDIR + f'landscape_null_{tis}_chr{c}.csv', index=False)
row.update({f'landscape_peaks_above_{k}_null': v for k, v in counts.items()}); row['landscape_frequencies_tested'] = nfreq
# ---- 2. Whittle fits of the covariance spectrum (per-chromosome periodograms averaged over samples)
spec = []
for c in np.unique(chrs):
    Y = D[chrs == c]; n = Y.shape[0]
    if n < 200: continue
    I = (np.abs(np.fft.rfft(Y, axis=0)[1:n // 2]) ** 2 / n).mean(1); spec.append((np.arange(1, n // 2) / n, I))
f_all = np.concatenate([s[0] for s in spec]); I_all = np.concatenate([s[1] for s in spec]); m_s = D.shape[1]
lor = lambda f, A, lam: A * (1 - np.exp(-2 / lam)) / (1 - 2 * np.exp(-1 / lam) * np.cos(2 * np.pi * f) + np.exp(-2 / lam))
models = {'one Lorentzian': (lambda f, p: lor(f, np.exp(p[0]), np.exp(p[1])) + np.exp(p[2]), [np.log(0.3), 0.0, np.log(0.5)]),
          'two Lorentzians': (lambda f, p: lor(f, np.exp(p[0]), np.exp(p[1])) + lor(f, np.exp(p[2]), np.exp(p[3])) + np.exp(p[4]), [np.log(0.3), 0.0, np.log(0.1), np.log(15), np.log(0.5)]),
          'power law': (lambda f, p: np.exp(p[0]) * f ** (-np.exp(p[1])) + np.exp(p[2]), [np.log(0.01), np.log(0.8), np.log(0.5)])}
fits = {}
for name, (S, p0) in models.items():
    nll = lambda p: m_s * np.sum(np.log(S(f_all, p)) + I_all / S(f_all, p))
    best = min((minimize(nll, np.array(p0) + d, method='Nelder-Mead', options={'maxiter': 20000, 'xatol': 1e-6, 'fatol': 1e-6}) for d in (0, 0.3, -0.3)), key=lambda r: r.fun)
    fits[name] = best; row[f'whittle_AIC_{name}'] = 2 * best.fun + 2 * len(p0)
lams = sorted([np.exp(fits['two Lorentzians'].x[1]), np.exp(fits['two Lorentzians'].x[3])])
row.update({'whittle_short_length_genes': lams[0], 'whittle_long_length_genes': lams[1],
            'whittle_dAIC_one_minus_two': row['whittle_AIC_one Lorentzian'] - row['whittle_AIC_two Lorentzians'],
            'whittle_dAIC_power_minus_two': row['whittle_AIC_power law'] - row['whittle_AIC_two Lorentzians']})
# ---- 3. coupling against base-pair distance
Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-12); edges = np.logspace(np.log10(5e3), np.log10(5e6), 21)
sums = np.zeros(20); cnt = np.zeros(20); far = []
for c in np.unique(chrs):
    idx = np.where(chrs == c)[0]; pos = tss[idx]; Zc = Z[idx]
    for i in range(len(idx)):
        d = np.abs(pos[i + 1:] - pos[i]); near = np.where(d <= 5e6)[0]
        if len(near):
            r = Zc[i + 1 + near] @ Zc[i] / Z.shape[1]; b = np.digitize(d[near], edges) - 1; ok = (b >= 0) & (b < 20)
            np.add.at(sums, b[ok], r[ok]); np.add.at(cnt, b[ok], 1)
        farj = np.where((d > 2e7) & (d < 4e7))[0]
        if len(farj): j = rng.choice(farj, min(3, len(farj)), replace=False); far.extend(Zc[i + 1 + j] @ Zc[i] / Z.shape[1])
base = float(np.mean(far)); mid = np.sqrt(edges[:-1] * edges[1:]); prof = sums / np.maximum(cnt, 1) - base
two_exp = lambda d, a1, l1, a2, l2: a1 * np.exp(-d / l1) + a2 * np.exp(-d / l2)
try:
    p, _ = curve_fit(two_exp, mid, prof, p0=[0.15, 3e4, 0.02, 1e6], sigma=1 / np.sqrt(cnt), bounds=([0, 1e3, 0, 1e5], [1, 5e5, 1, 2e7]), maxfev=20000)
    row.update({'bp_short_length_kb': p[1] / 1e3, 'bp_long_length_Mb': p[3] / 1e6, 'bp_short_amplitude': p[0], 'bp_long_amplitude': p[2]})
except Exception as e: row['bp_fit_error'] = str(e)[:80]
row.update({'bp_far_baseline': base, 'median_gene_spacing_kb': float(np.median(np.diff(tss)[np.diff(tss) > 0]) / 1e3)})
pd.DataFrame({'tissue': tis, 'distance_bp': mid, 'pairs': cnt, 'coupling_minus_far': prof}).to_csv(OUTDIR + f'coupling_by_bp_{tis}.csv', index=False)
replace_groups_csv(pd.DataFrame([row]), OUTDIR + 'spectral_rigour.csv', ['tissue'], [tis], key=['tissue'])
print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})
