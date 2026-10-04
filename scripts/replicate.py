import numpy as np, pandas as pd, statsmodels.api as sm
src = open('cohesin_v2.py').read(); src = src[:src.index('rows = []')]
src = src.replace("keep = np.array([s[13:15] == '01' for s in samp])", "keep = np.array([s[13:15] in ('01', '03') for s in samp])")
src = src.replace("SUB = {", "SUB = {'LAML': (['MPO', 'ELANE', 'AZU1'], ['CD14', 'LYZ', 'CSF1R']), 'GBM': (['OLIG2', 'SOX2', 'PDGFRA'], ['CHI3L1', 'CD44', 'MET']), ")
exec(src)
COH = ['STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']; rows = []
for c in ['LAML', 'GBM']:
    D, chrs, samp, cov = prepare(c); e = scores(D, chrs)
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(COH)].Sample_ID)).astype(float)
        Xc = cov[['cna_burden', 'log_tmb', 'immune', 'stromal', 'subtype']]; Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join((Xc - Xc.mean()) / Xc.std()))
        base = e[m == 0].mean(); row = {'cohort': c, 'class': mclass, 'n_tumours': len(e), 'n_mutant': int(m.sum()), 'raw_pct': 100 * (e[m == 1].mean() - base) / base if m.sum() else np.nan}
        if m.sum() >= 3:
            fit = sm.OLS(e, Xd).fit(cov_type='HC3'); r = e - sm.OLS(e, sm.add_constant(((Xc - Xc.mean()) / Xc.std()).values)).fit().fittedvalues
            obs = r[m == 1].mean() - r[m == 0].mean(); rng = np.random.default_rng(0); null = np.array([(lambda p: r[p == 1].mean() - r[p == 0].mean())(rng.permutation(m)) for _ in range(10000)])
            row.update({'adj_pct': 100 * fit.params['mutant'] / base, 'ci_low': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high': 100 * fit.conf_int().loc['mutant', 1] / base, 'p_perm': (np.sum(np.abs(null) >= abs(obs)) + 1) / 10001})
        rows.append(row)
R = pd.DataFrame(rows); R.to_csv('replication_LAML_GBM.csv', index=False); print(R.round(3).to_string(index=False))
