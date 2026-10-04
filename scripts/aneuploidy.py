import numpy as np, pandas as pd, statsmodels.api as sm
src = open('cohesin_v2.py').read(); src = src[:src.index('rows = []')]
src = src.replace("keep = np.array([s[13:15] == '01' for s in samp])", "keep = np.array([s[13:15] in ('01', '03') for s in samp])")
src = src.replace("    for i in range(D.shape[0]):\n        cn = CNv[i]", "    D_raw = D.copy()\n    for i in range(D.shape[0]):\n        cn = CNv[i]")
src = src.replace("    return D, g.chr.values, samp, cov", "    return D, D_raw, g.chr.values, samp, cov")
src = src.replace("SUB = {", "SUB = {'LAML': (['MPO', 'ELANE', 'AZU1'], ['CD14', 'LYZ', 'CSF1R']), 'GBM': (['OLIG2', 'SOX2', 'PDGFRA'], ['CHI3L1', 'CD44', 'MET']), ")
exec(src)
def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
rows = []
for c in ['BLCA', 'UCEC', 'STAD', 'COAD', 'GBM']:
    D, Draw, chrs, samp, cov = prepare(c); cna = cov.cna_burden.values
    for name, M in [('raw', Draw), ('copy-number corrected', D)]:
        far = lagscore(M, chrs, (20, 30)); near = lagscore(M, chrs, (1, 2, 3))
        rows.append({'cohort': c, 'expression': name, 'tumours': len(far), 'mean_far': far.mean(), 'mean_near': near.mean(),
                     'rho_far_cna': pd.Series(far).corr(pd.Series(cna), method='spearman'), 'rho_cis_excess_cna': pd.Series(near - far).corr(pd.Series(cna), method='spearman'), 'cohort_mean_cna': cna.mean()})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv('aneuploidy_scaling.csv', index=False); pd.set_option('display.width', 220); print(R.round(3).to_string(index=False))
raw = R[R.expression == 'raw']; print('\nacross cohorts (raw): Spearman(mean far floor, mean CNA burden) = %.2f' % raw.mean_far.corr(raw.cohort_mean_cna, method='spearman'))
