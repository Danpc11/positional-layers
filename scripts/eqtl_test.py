import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import glob, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
norm = lambda s: s.lower().replace('-', '_')
BINS = [-np.inf, 0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
tis = sorted(set(P) & set(S)); print('tissues with coupling and eQTL:', len(tis))
rows = []
for t in tis:
    d = P[t].merge(S[t], on=['g1', 'g2']); d = d[d.dist > 0]                       # non-overlapping pairs only (avoids read sharing)
    d['bin'] = pd.cut(d.dist, BINS).astype(str); be = d[d.egene1 & d.egene2].copy()
    be['same'] = be.share_same.astype(int); be['opp_only'] = (be.share_any & ~be.share_same).astype(int)
    m = smf.ols('r ~ same + opp_only + C(bin) + C(orientation)', data=be).fit()
    grp = lambda q: d.loc[q, 'r'].mean()
    r_all = d.r.mean(); delta = m.params['same']; f_same = be.same.sum() / len(d)
    rows.append({'tissue': t, 'pairs': len(d), 'both_eGenes': len(be), 'frac_share_same_of_all_pairs': f_same,
                 'r_share_same': grp(d.egene1 & d.egene2 & d.share_same), 'r_eGenes_no_share': grp(d.egene1 & d.egene2 & ~d.share_any), 'r_not_both_eGenes': grp(~(d.egene1 & d.egene2)),
                 'delta_same_adj': delta, 'p_same': m.pvalues['same'], 'delta_opposite_only_adj': m.params['opp_only'],
                 'share_of_cis_from_shared_eQTL': f_same * delta / r_all})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'eqtl_cis_test.csv', index=False)
pd.set_option('display.width', 250); print(R.round(3).to_string(index=False))
print('\nmedian over tissues:'); print(R.drop(columns='tissue').median().round(4).to_string())
# tissue specificity: own-tissue sharing vs other-tissue sharing, on pairs that are eGene pairs in both tissues
spec = []
for t in tis:
    d0 = P[t][P[t].dist > 0][['g1', 'g2', 'r', 'dist', 'orientation']]; d0['bin'] = pd.cut(d0.dist, BINS).astype(str)
    own = S[t][['g1', 'g2', 'egene1', 'egene2', 'share_same']].rename(columns={'share_same': 'own', 'egene1': 'oe1', 'egene2': 'oe2'})
    for u in tis:
        if u == t: continue
        oth = S[u][['g1', 'g2', 'egene1', 'egene2', 'share_same']].rename(columns={'share_same': 'other'})
        x = d0.merge(own, on=['g1', 'g2']).merge(oth, on=['g1', 'g2']); x = x[x.oe1 & x.oe2 & x.egene1 & x.egene2]
        if len(x) < 500: continue
        x['own'] = x.own.astype(int); x['other'] = x.other.astype(int)
        m = smf.ols('r ~ own + other + C(bin) + C(orientation)', data=x).fit()
        spec.append({'coupling_tissue': t, 'eqtl_tissue': u, 'n': len(x), 'b_own': m.params['own'], 'b_other': m.params['other']})
SP = pd.DataFrame(spec); SP.to_csv(OUTDIR + 'eqtl_tissue_specificity.csv', index=False)
print(f'\ntissue specificity (joint models, {len(SP)} tissue pairs): own-tissue sharing coefficient median {SP.b_own.median():.3f}; other-tissue sharing coefficient median {SP.b_other.median():.3f}; own > other in {100*(SP.b_own > SP.b_other).mean():.0f}% of tissue pairs')
