"""Sensitivity of the tumour cis-excess analysis to the baseline lags.

The cis excess of each tumour is the mean standardised product of genes 1-3 positions apart minus that of genes at the
baseline lags. The paper uses 20-30 genes, which lies within the domain scale of coupling (7-24 genes in GTEx); here the
baseline is moved beyond it (60, 70 and 80 genes). Primary specification of cohesin_freedman_lane.py (any coding
mutation, tumours above the 90th percentile of mutation burden excluded), same covariates and Freedman-Lane permutation.
Usage: python tumour_baseline_sensitivity.py BLCA,UCEC [permutations]   Output: OUTDIR/tumour_baseline_sensitivity.csv
"""
import sys
import numpy as np, pandas as pd, statsmodels.api as sm
from poslayers.config import OUTDIR, replace_groups_csv
from lib_tumour import coding, prepare
from cohesin_freedman_lane_lib import lagscore

B = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
for c, gene in [x for x in [('BLCA', 'STAG2'), ('UCEC', 'CTCF')] if x[0] in sys.argv[1].split(',')]:
    D, chrs, samp, cov = prepare(c, cn='continuous'); short = lagscore(D, chrs, (1, 2, 3)); rows = []
    keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, 0.9)
    m = pd.Index(samp).isin(set(coding[coding.gene == gene].Sample_ID)).astype(float)[keep]
    Xs = cov[keep]; Xs = (Xs - Xs.mean()) / Xs.std(); Xfull = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xs.index).join(Xs)).values; Xred = sm.add_constant(Xs.values)
    for label, base_lags in (('20-30 genes (paper)', (20, 30)), ('60-80 genes', (60, 70, 80))):
        y = (short - lagscore(D, chrs, base_lags))[keep]
        full = sm.OLS(y, Xfull).fit(cov_type='HC3'); t_obs = full.tvalues[1]; base = y[m == 0].mean()
        red = sm.OLS(y, Xred).fit(); rng = np.random.default_rng(1)
        tnull = np.array([sm.OLS(red.fittedvalues + rng.permutation(red.resid), Xfull).fit(cov_type='HC3').tvalues[1] for _ in range(B)])
        rows.append({'cohort': c, 'gene': gene, 'baseline': label, 'n_mut': int(m.sum()), 'n': len(y), 'wild_type_mean': base,
                     'adj_pct': 100 * full.params[1] / base, 'ci_low': 100 * full.conf_int()[1, 0] / base, 'ci_high': 100 * full.conf_int()[1, 1] / base,
                     'p_freedman_lane': (np.sum(np.abs(tnull) >= abs(t_obs)) + 1) / (B + 1)})
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    replace_groups_csv(pd.DataFrame(rows), OUTDIR + 'tumour_baseline_sensitivity.csv', ['cohort'], [c], key=['cohort', 'baseline'])
