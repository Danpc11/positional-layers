"""Spectral form of the cis layer. Residual expression (GC + technical covariates removed), standardised per gene.
Route A: autocorrelation rho(L) across samples -> fit exponential(s) -> lambda_acf.
Route B: mean periodogram of the residual profiles along gene order -> fit discrete Lorentzian(s) + floor -> lambda_spec;
compare with a power law + floor. Theory: exponential ACF <=> Lorentzian spectrum with the same lambda, no discrete peaks."""
import sys, os, numpy as np, pandas as pd, pyannotables as pa
from scipy.optimize import curve_fit
CHR = [str(i) for i in range(1, 23)] + ['X']; norm = lambda s: s.lower().replace('-', '_')
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
SA = pd.read_csv('/home/claude/repo/data/external/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
def resid_samples(D, F): F1 = np.column_stack([np.ones(F.shape[0]), F]); B = np.linalg.lstsq(F1, D.T, rcond=None)[0]; return (D.T - F1 @ B).T
def resid_genes(D, G): G1 = np.column_stack([np.ones(G.shape[0]), G]); B = np.linalg.lstsq(G1, D, rcond=None)[0]; return D - G1 @ B
def lor1(f, a, lam, c):     # discrete Lorentzian (spectrum of an exponential ACF with length lam) + white floor
    q = np.exp(-1 / lam); return c + a * (1 - q ** 2) / (1 - 2 * q * np.cos(2 * np.pi * f) + q ** 2)
def lor2(f, a1, l1, a2, l2, c): return lor1(f, a1, l1, 0) + lor1(f, a2, l2, 0) + c
def plaw(f, k, al, c): return c + k * f ** (-al)
OUT = '/home/claude/atlas/lorentz_results.csv'
for t in sys.argv[1:]:
    f = f'/home/claude/repo/data/external/gtex/gene_reads_adult_gtex_v11_{t}_gct.gz'
    C = pd.read_csv(f, sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
    g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); Yc = Y - Y.mean(0, keepdims=True); D = Yc - Yc.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); D = resid_genes(D, np.column_stack([gcz, gcz ** 2]))
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    tech = np.column_stack([rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    if tech.shape[1] < D.shape[1] - 10: D = resid_samples(D, tech)
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    # route A: autocorrelation
    lags = np.arange(1, 61); rho = np.array([np.mean(np.concatenate([np.mean(Z[i[:-L]] * Z[i[L:]], 1) for c in CHR for i in [np.where(chrs == c)[0]] if len(i) > L + 5])) for L in lags])
    fl = rho[40:].mean(); y = rho - fl
    e1 = curve_fit(lambda L, a, lam: a * np.exp(-L / lam), lags[:30], y[:30], p0=[0.15, 5], maxfev=20000)[0]
    e2 = curve_fit(lambda L, a1, l1, a2, l2: a1 * np.exp(-L / l1) + a2 * np.exp(-L / l2), lags[:30], y[:30], p0=[0.1, 1, 0.05, 10], bounds=([0, 0.2, 0, 2], [1, 5, 1, 100]), maxfev=20000)[0]
    rss1 = np.sum((y[:30] - e1[0] * np.exp(-lags[:30] / e1[1])) ** 2); rss2 = np.sum((y[:30] - (e2[0] * np.exp(-lags[:30] / e2[1]) + e2[2] * np.exp(-lags[:30] / e2[3]))) ** 2)
    # route B: mean periodogram of residual profiles (per chromosome, normalised by length), on a common log-frequency grid
    bins = np.logspace(np.log10(1 / 1500), np.log10(0.5), 61); acc = np.zeros(60); cnt = np.zeros(60)
    for c in CHR:
        i = np.where(chrs == c)[0]; n = len(i)
        if n < 200: continue
        P = np.mean(np.abs(np.fft.rfft(Z[i], axis=0)[1:n // 2 + 1]) ** 2, axis=1) / n; fr = np.arange(1, n // 2 + 1) / n
        k = np.digitize(fr, bins) - 1; ok = (k >= 0) & (k < 60); np.add.at(acc, k[ok], P[ok]); np.add.at(cnt, k[ok], 1)
    m = cnt > 0; fx = np.sqrt(bins[:-1] * bins[1:])[m]; S = acc[m] / cnt[m]
    def fit(fun, p0, bnds):
        p = curve_fit(fun, fx, S, p0=p0, bounds=bnds, maxfev=50000)[0]; r = np.log(S) - np.log(fun(fx, *p)); return p, np.sum(r ** 2)
    pL1, rL1 = fit(lor1, [1, 5, 0.8], ([0, 0.3, 0], [100, 500, 10])); pL2, rL2 = fit(lor2, [0.3, 1, 0.5, 15, 0.8], ([0, 0.2, 0, 2, 0], [100, 5, 100, 500, 10])); pP, rP = fit(plaw, [0.01, 0.5, 0.8], ([0, 0, 0], [10, 3, 10]))
    nf = len(fx); aic = lambda rss, k: nf * np.log(rss / nf) + 2 * k
    resid = S / lor2(fx, *pL2); row = {'tissue': t, 'samples': Z.shape[1], 'lam_acf_1exp': e1[1], 'lam_acf_short': e2[1], 'lam_acf_long': e2[3], 'acf_rss_1exp': rss1, 'acf_rss_2exp': rss2,
        'lam_spec_1lor': pL1[1], 'lam_spec_short': pL2[1], 'lam_spec_long': pL2[3], 'AIC_lor1': aic(rL1, 3), 'AIC_lor2': aic(rL2, 5), 'AIC_powerlaw': aic(rP, 3), 'powerlaw_alpha': pP[1],
        'max_peak_over_fit': resid.max(), 'n_bins_over_1.5x_fit': int((resid > 1.5).sum())}
    pd.DataFrame([row]).to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
    np.save(f'/home/claude/atlas/spec_{t}.npy', np.vstack([fx, S, lor2(fx, *pL2), plaw(fx, *pP)]))
    print(t, {k: round(v, 3) if isinstance(v, float) else v for k, v in row.items() if k != 'tissue'}, flush=True)
