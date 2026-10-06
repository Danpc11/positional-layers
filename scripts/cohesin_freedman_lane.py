"""Freedman-Lane permutation for the cohesin/CTCF effect on cis excess.
Reduced model y ~ covariates gives fitted values and residuals; each permutation shuffles the RESIDUALS, adds them back to the
reduced fit, refits the FULL model y* ~ mutant + covariates and records the HC3 t statistic of 'mutant'. The p-value compares the
observed HC3 t with that null, so the permutation tests the same estimand that is reported, and the covariate design is kept.
Specification pre-declared as primary: any coding mutation, tumours above the 90th percentile of mutation burden excluded.
The other three specifications are sensitivity analyses."""
import os, sys
from poslayers.config import OUTDIR
from lib_tumour import coding, prepare, trunc
from cohesin_freedman_lane_lib import lagscore
import numpy as np, pandas as pd, statsmodels.api as sm
B = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
rows = []
for c, gm in [x for x in [('BLCA', ['STAG2']), ('UCEC', ['CTCF'])] if x[0] in sys.argv[1].split(',')]:
    D, chrs, samp, cov = prepare(c, cn='continuous'); e = lagscore(D, chrs, (1, 2, 3)) - lagscore(D, chrs, (20, 30))
    np.save(OUTDIR + f'cisexcess_{c}.npy', e); cov.assign(cis_excess=e).to_csv(OUTDIR + f'tumour_scores_{c}.csv')
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, hq); y = e[keep]
            m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(gm)].Sample_ID)).astype(float)[keep]
            Xs = cov[keep]; Xs = (Xs - Xs.mean()) / Xs.std()
            Xfull = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xs.index).join(Xs)).values; Xred = sm.add_constant(Xs.values)
            full = sm.OLS(y, Xfull).fit(cov_type='HC3'); t_obs = full.tvalues[1]; base = y[m == 0].mean()
            red = sm.OLS(y, Xred).fit(); fit_r, res_r = red.fittedvalues, red.resid
            rng = np.random.default_rng(1); tnull = np.empty(B)
            for b in range(B): tnull[b] = sm.OLS(fit_r + rng.permutation(res_r), Xfull).fit(cov_type='HC3').tvalues[1]
            rows.append({'cohort': c, 'gene': gm[0], 'class': mclass, 'tmb_q': hq, 'primary': mclass == 'all coding' and hq == 0.9,
                         'n_mut': int(m.sum()), 'n': len(y), 'adj_pct': 100 * full.params[1] / base,
                         'ci_low': 100 * full.conf_int()[1, 0] / base, 'ci_high': 100 * full.conf_int()[1, 1] / base,
                         'p_HC3': full.pvalues[1], 'p_freedman_lane': (np.sum(np.abs(tnull) >= abs(t_obs)) + 1) / (B + 1),
                         'corr_mutant_with_covariates_max': float(np.max(np.abs([np.corrcoef(m, Xs[k])[0, 1] for k in Xs.columns])))})
            print(rows[-1], flush=True)
# Idempotent output: rows for the cohorts run now replace any earlier rows for the same cohort; one row per specification.
out = OUTDIR + 'cohesin_freedman_lane.csv'; new = pd.DataFrame(rows)
if os.path.exists(out):
    prev = pd.read_csv(out); new = pd.concat([prev[~prev.cohort.isin(new.cohort.unique())], new], ignore_index=True)
key = ['cohort', 'gene', 'class', 'tmb_q']
assert not new.duplicated(key).any(), 'duplicate specifications in cohesin_freedman_lane.csv'
new.to_csv(out, index=False)
