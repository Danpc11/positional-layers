import numpy as np, pandas as pd, statsmodels.api as sm
src = open('cohesin_v2.py').read(); src = src[:src.index('rows = []')]
i0 = src.index("    hdr = pd.read_csv('GDC-PANCAN.gistic.tsv'"); i1 = src.index("    sym = g.sym.astype(str).values")
src = src[:i0] + '''    Zc = np.load(f'cn_cont_{c}.npz', allow_pickle=True); CNd = pd.DataFrame(Zc['CN'], index=Zc['genes'], columns=Zc['samples'])
    have = pd.Index(samp).isin(CNd.columns); D = D[:, have]; Y = Y[:, have]; samp = samp[have]
    CNv = CNd.reindex(index=g.iloc[:, 0].values, columns=samp).fillna(0).values
    global D_RAW; D_RAW = D.copy()
    for i in range(D.shape[0]):                          # remove each gene's own continuous dosage effect (linear + quadratic)
        x = CNv[i]
        if np.std(x) < 1e-6: continue
        X1 = np.column_stack([np.ones(len(x)), x, x ** 2]); D[i] -= X1 @ np.linalg.lstsq(X1, D[i], rcond=None)[0]
''' + src[i1:]
src = src.replace("'cna_burden': np.mean(CNv != 0, 0)", "'cna_burden': np.mean(np.abs(CNv) > 0.2, 0)")
src = src.replace("SUB = {", "SUB = {'GBM': (['OLIG2', 'SOX2', 'PDGFRA'], ['CHI3L1', 'CD44', 'MET']), ")
exec(src)
import sys, os
WANT = sys.argv[1].split(',')
def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
an, te = [], []
for c, gm in [x for x in [('BLCA', ['STAG2']), ('UCEC', ['CTCF']), ('STAD', None), ('COAD', None), ('GBM', None)] if x[0] in WANT]:
    D, chrs, samp, cov = prepare(c); cna = cov.cna_burden.values
    for nm, M in [('raw', D_RAW), ('continuous CN corrected', D)]:
        far = lagscore(M, chrs, (20, 30)); near = lagscore(M, chrs, (1, 2, 3))
        an.append({'cohort': c, 'expression': nm, 'mean_far': far.mean(), 'rho_far_cna': pd.Series(far).corr(pd.Series(cna), method='spearman'), 'rho_cisexcess_cna': pd.Series(near - far).corr(pd.Series(cna), method='spearman')})
    if gm is None: continue
    e = lagscore(D, chrs, (1, 2, 3)) - lagscore(D, chrs, (20, 30))
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, hq); y = e[keep]; m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(gm)].Sample_ID)).astype(float)[keep]
            Xc = cov[keep]; Xs = (Xc - Xc.mean()) / Xc.std(); Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join(Xs))
            fit = sm.OLS(y, Xd).fit(cov_type='HC3'); base = y[m == 0].mean(); r = y - sm.OLS(y, sm.add_constant(Xs.values)).fit().fittedvalues
            obs = r[m == 1].mean() - r[m == 0].mean(); rng = np.random.default_rng(0); null = np.array([(lambda p: r[p == 1].mean() - r[p == 0].mean())(rng.permutation(m)) for _ in range(10000)])
            te.append({'cohort': c, 'gene': gm[0], 'class': mclass, 'tmb_q': hq, 'n_mut': int(m.sum()), 'adj_pct': 100 * fit.params['mutant'] / base,
                       'ci_low': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high': 100 * fit.conf_int().loc['mutant', 1] / base, 'p_perm': (np.sum(np.abs(null) >= abs(obs)) + 1) / 10001})
    print(c, 'done', flush=True)
A = pd.DataFrame(an); T = pd.DataFrame(te)
A.to_csv('aneuploidy_continuousCN.csv', mode='a', header=not os.path.exists('aneuploidy_continuousCN.csv'), index=False)
if len(T): T.to_csv('cohesin_continuousCN.csv', mode='a', header=not os.path.exists('cohesin_continuousCN.csv'), index=False)
pd.set_option('display.width', 200); print(A.round(3).to_string(index=False)); print(T.round(3).to_string(index=False))
