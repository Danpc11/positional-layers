import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # scripts/ for lib_*.py
import os, sys
from poslayers.config import OUTDIR
import numpy as np, pandas as pd, statsmodels.formula.api as smf
D = pd.read_csv(OUTDIR + 'crispri_pairs.csv.gz'); D = D[np.isfinite(D.resp) & np.isfinite(D.r)]; D['resp'] = D.resp.clip(-20, 20); D['rk'] = D.r * D.kd
CI = D[D.type == 'cis'].copy(); TR = D[D.type == 'trans'].copy()
CI['dbin'] = pd.cut(CI.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['a<10kb', 'b10-50kb', 'c50-200kb', 'd200-500kb', 'e0.5-1Mb']).astype(str)
print(f'cis pairs {len(CI):,} | trans pairs {len(TR):,} | perturbations {D.pert.nunique():,}')
gc = lambda d: {'groups': d.pert.astype('category').cat.codes}
mc = smf.ols('resp ~ rk + r + kd + C(dbin) + C(dbin):kd', data=CI).fit(cov_type='cluster', cov_kwds=gc(CI))
mt = smf.ols('resp ~ rk + r + kd', data=TR).fit(cov_type='cluster', cov_kwds=gc(TR))
for nm, m in (('cis', mc), ('trans', mt)):
    ci = m.conf_int().loc['rk']; print(f"{nm:5s}: coupling x knockdown {m.params['rk']:.3f} (95% CI {ci[0]:.3f} to {ci[1]:.3f}), P = {m.pvalues['rk']:.1e}")
# by distance: does coupling matter beyond the KRAB spreading distance?
for b, d in CI.groupby('dbin'):
    m = smf.ols('resp ~ rk + r + kd', data=d).fit(cov_type='cluster', cov_kwds=gc(d)); print(f"  {b[1:]:>10}: n={len(d):,} slope {m.params['rk']:.3f} P={m.pvalues['rk']:.1e} | mean response {d.resp.mean():.3f}")
S = CI[CI.kd > 2].copy(); qs = S.r.quantile([1 / 3, 2 / 3]).values
S['rt'] = pd.cut(S.r, [-1, qs[0], qs[1], 1], labels=['low r', 'mid r', 'high r']); TS = TR[TR.kd > 2].copy(); TS['rt'] = pd.cut(TS.r, [-1, qs[0], qs[1], 1], labels=['low r', 'mid r', 'high r'])
print('\nstrong knockdowns (>4-fold): mean response by coupling tertile')
tab = S.pivot_table(index='dbin', columns='rt', values='resp', aggfunc='mean', observed=True); tab.loc['trans (other chromosomes)'] = TS.groupby('rt', observed=True).resp.mean(); print(tab.round(3).to_string())

# The trans control is random, not matched on coupling: report both distributions so the comparison can be judged.
dist = pd.DataFrame({'cis': CI.r.describe(percentiles=[.1, .25, .5, .75, .9]), 'trans': TR.r.describe(percentiles=[.1, .25, .5, .75, .9])})
dist.to_csv(OUTDIR + 'crispri_coupling_distribution_cis_vs_trans.csv'); print('coupling distribution, cis vs random trans:'); print(dist.round(3).to_string())
