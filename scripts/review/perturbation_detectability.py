"""How large a coupling-dependent neighbour response could the perturbation analyses have detected?

For CRISPRi (cis and trans) the coupling x knockdown coefficient and its perturbation-clustered SE are refitted from the pair
table; for each drug comparison the slope SE comes from the genomic-block bootstrap in drug_test.py. For each test:
  mde80  : minimum detectable effect at 80% power, two-sided alpha 0.05, = 2.80 x SE;
  tost_* : two one-sided tests of equivalence at a margin equal to the MDE, and (CRISPRi cis) at the trans estimate,
           i.e. could an effect as large as the one seen in trans have gone undetected in cis?
Output: OUTDIR/perturbation_detectability.csv
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # scripts/ for lib_*.py
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
from poslayers.config import OUTDIR

D = pd.read_csv(OUTDIR + 'crispri_pairs.csv.gz'); D = D[np.isfinite(D.resp) & np.isfinite(D.r)]; D['resp'] = D.resp.clip(-20, 20); D['rk'] = D.r * D.kd
CI = D[D.type == 'cis'].copy(); TR = D[D.type == 'trans'].copy()
CI['dbin'] = pd.cut(CI.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['a', 'b', 'c', 'd', 'e']).astype(str)
gc = lambda d: {'groups': d.pert.astype('category').cat.codes}
mc = smf.ols('resp ~ rk + r + kd + C(dbin) + C(dbin):kd', data=CI).fit(cov_type='cluster', cov_kwds=gc(CI))
mt = smf.ols('resp ~ rk + r + kd', data=TR).fit(cov_type='cluster', cov_kwds=gc(TR))
def tost(est, se, margin):
    return max(stats.norm.sf((est + margin) / se), stats.norm.cdf((est - margin) / se))
rows = []
for name, m in (('CRISPRi cis', mc), ('CRISPRi trans', mt)):
    est, se = m.params['rk'], m.bse['rk']; mde = 2.80 * se
    rows.append({'test': name, 'estimate': est, 'se': se, 'ci_low': est - 1.96 * se, 'ci_high': est + 1.96 * se, 'mde80': mde, 'p_tost_at_mde': tost(est, se, mde)})
rows[0]['trans_estimate_as_margin'] = mt.params['rk']; rows[0]['p_tost_at_trans_estimate'] = tost(mc.params['rk'], mc.bse['rk'], abs(mt.params['rk']))
Dr = pd.read_csv(OUTDIR + 'drug_cis_results_thr3.csv')
for _, r in Dr.iterrows():
    se = (r.slope_ci_high - r.slope_ci_low) / (2 * 1.96)
    rows.append({'test': f'{r.cell} · {r.perturbation}', 'estimate': r.slope_on_baseline_coupling, 'se': se, 'ci_low': r.slope_ci_low, 'ci_high': r.slope_ci_high,
                 'mde80': 2.80 * se, 'p_tost_at_mde': tost(r.slope_on_baseline_coupling, se, 2.80 * se)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'perturbation_detectability.csv', index=False); print(R.round(4).to_string(index=False))
