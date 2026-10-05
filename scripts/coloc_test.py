import os
from poslayers.config import OUTDIR
from lib_boot import block_bootstrap
import glob, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
norm = lambda s: s.lower().replace('-', '_'); BINS = [-np.inf, 0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
Cc = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'coloc_*.csv.gz')}
rows, dose = [], []
for t in sorted(set(P) & set(S) & set(Cc)):
    d = P[t].merge(S[t], on=['g1', 'g2']).merge(Cc[t], on=['g1', 'g2']); d = d[(d.dist > 0) & d.fm1 & d.fm2].copy()
    d['bin'] = pd.cut(d.dist, BINS).astype(str)
    d['coloc_same'] = ((d.p_coloc >= 0.5) & (d.direction == 1)).astype(int); d['coloc_opp'] = ((d.p_coloc >= 0.5) & (d.direction == -1)).astype(int)
    d['sig_same'] = d.share_same.astype(int); d['sig_opp'] = (d.share_any & ~d.share_same).astype(int)
    m1 = smf.ols('r ~ sig_same + sig_opp + C(bin) + C(orientation)', data=d).fit()
    m2 = smf.ols('r ~ coloc_same + coloc_opp + sig_same + sig_opp + C(bin) + C(orientation)', data=d).fit()
    bb = block_bootstrap('r ~ coloc_same + coloc_opp + sig_same + sig_opp + C(bin) + C(orientation)', d, ['coloc_same'], seed=2)
    rows.append({'tissue': t, 'fine_mapped_pairs': len(d), 'n_coloc_same': int(d.coloc_same.sum()), 'n_coloc_opp': int(d.coloc_opp.sum()),
                 'r_coloc_same': d.r[d.coloc_same == 1].mean(), 'r_sig_shared_not_coloc': d.r[(d.sig_same == 1) & (d.p_coloc < 0.1)].mean(), 'r_no_sharing': d.r[(d.sig_same == 0) & (d.sig_opp == 0)].mean(),
                 'b_sig_same_alone': m1.params['sig_same'], 'b_coloc_same': m2.params['coloc_same'], 'coloc_ci_low': bb['coloc_same']['ci_low'], 'coloc_ci_high': bb['coloc_same']['ci_high'], 'p_coloc_block': bb['coloc_same']['p_block'], 'p_coloc_ols_pairs_independent': m2.pvalues['coloc_same'],
                 'b_coloc_opp': m2.params['coloc_opp'], 'b_sig_same_given_coloc': m2.params['sig_same']})
    sm = d[(d.direction == 1) | (d.p_coloc == 0)]
    for lo, hi, lab in [(-1, 1e-9, '0'), (1e-9, 0.1, '0-0.1'), (0.1, 0.5, '0.1-0.5'), (0.5, 0.8, '0.5-0.8'), (0.8, 1.01, '>0.8')]:
        x = sm[(sm.p_coloc > lo) & (sm.p_coloc <= hi)]; dose.append({'tissue': t, 'p_coloc_bin': lab, 'n': len(x), 'mean_r': x.r.mean()})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'coloc_cis_test.csv', index=False)
pd.set_option('display.width', 250)
print(R[['tissue', 'fine_mapped_pairs', 'n_coloc_same', 'n_coloc_opp', 'r_coloc_same', 'r_sig_shared_not_coloc', 'r_no_sharing', 'b_sig_same_alone', 'b_coloc_same', 'b_coloc_opp', 'b_sig_same_given_coloc']].round(3).to_string(index=False))
print('\nmedians:'); print(R.drop(columns='tissue').median().round(4).to_string())
print(f"\ncoloc_same > 0 in {(R.b_coloc_same > 0).sum()}/{len(R)} tissues (P<0.05 in {(R.p_coloc_block < 0.05).sum()}); coloc_opp < 0 in {(R.b_coloc_opp < 0).sum()}/{len(R)}")
print(f"significant-variant sharing coefficient: alone {R.b_sig_same_alone.median():.3f} -> with colocalisation in the model {R.b_sig_same_given_coloc.median():.3f}")
D = pd.DataFrame(dose); D.to_csv(OUTDIR + 'coloc_dose_response.csv', index=False)
print('\ndose-response (same-direction or none), median over tissues of mean r:'); print(D.groupby('p_coloc_bin', sort=False).agg(n=('n', 'sum'), median_r=('mean_r', 'median')).round(3).to_string())
