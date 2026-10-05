from poslayers.config import OUTDIR
import numpy as np
import pandas as pd
from lib_tumour import CHR, prepare

def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
rows = []
for c in ['BLCA', 'UCEC', 'STAD', 'COAD', 'GBM']:
    D, Draw, chrs, samp, cov = prepare(c, codes=('01', '03'), return_raw=True); cna = cov.cna_burden.values
    for name, M in [('raw', Draw), ('copy-number corrected', D)]:
        far = lagscore(M, chrs, (20, 30)); near = lagscore(M, chrs, (1, 2, 3))
        rows.append({'cohort': c, 'expression': name, 'tumours': len(far), 'mean_far': far.mean(), 'mean_near': near.mean(),
                     'rho_far_cna': pd.Series(far).corr(pd.Series(cna), method='spearman'), 'rho_cis_excess_cna': pd.Series(near - far).corr(pd.Series(cna), method='spearman'), 'cohort_mean_cna': cna.mean()})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'aneuploidy_scaling.csv', index=False); pd.set_option('display.width', 220); print(R.round(3).to_string(index=False))
raw = R[R.expression == 'raw']; print('\nacross cohorts (raw): Spearman(mean far floor, mean CNA burden) = %.2f' % raw.mean_far.corr(raw.cohort_mean_cna, method='spearman'))
