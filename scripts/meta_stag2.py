"""Review point 4: inverse-variance combination of the STAG2 effect in bladder cancer (TCGA) and AML (BeatAML2), computed from the
exported estimates rather than typed into the figure. Primary specification in both cohorts: STAG2, any coding mutation.
Percent effects are combined on the percent scale; SE is recovered from the HC3 confidence interval (BLCA) or reported directly (AML)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd
from scipy import stats
b = pd.read_csv(OUTDIR + 'cohesin_freedman_lane.csv'); b = b[(b.cohort == 'BLCA') & b.primary].iloc[0]
a = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv'); a = a[a.group == 'STAG2, any coding'].iloc[0]
est = np.array([b.adj_pct, a.adj_pct]); se = np.array([(b.ci_high - b.ci_low) / (2 * 1.959964), a.se_pct]); w = 1 / se ** 2
m = np.sum(w * est) / w.sum(); s = 1 / np.sqrt(w.sum()); z = m / s; Q = np.sum(w * (est - m) ** 2)
row = {'combination': 'STAG2 BLCA + AML, any coding', 'blca_pct': b.adj_pct, 'blca_se': se[0], 'aml_pct': a.adj_pct, 'aml_se': se[1],
       'meta_pct': m, 'meta_se': s, 'ci_low': m - 1.959964 * s, 'ci_high': m + 1.959964 * s, 'p': 2 * stats.norm.sf(abs(z)),
       'cochran_Q': Q, 'p_heterogeneity': stats.chi2.sf(Q, 1)}
pd.DataFrame([row]).to_csv(OUTDIR + 'stag2_meta.csv', index=False); print({k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()})
