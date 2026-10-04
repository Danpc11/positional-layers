"""Does the cis layer respect TAD boundaries? Adjacent gene pairs within one TAD vs across a boundary, at matched distance,
using the tissue's own TAD partition (McArthur & Capra 20-bin landscape, hg19; genes placed with GRCh37 coordinates)."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import glob, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
from lib_tad import CHR, MATCH, tads, assign
TADA = {n: assign(tads(n)) for n in sorted(set(MATCH.values()))}
norm = lambda s: s.lower().replace('-', '_'); BINS = [0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
rows, spec = [], []
for t, tn in MATCH.items():
    d = P[t].merge(S[t][['g1', 'g2', 'egene1', 'egene2', 'share_any', 'share_same']], on=['g1', 'g2']); d = d[d.dist > 0].copy(); d['bin'] = pd.cut(d.dist, BINS).astype(str)
    def fit(tadname, sub=None):
        a = TADA[tadname]; x = d.copy(); x['t1'] = a.reindex(x.g1).values; x['t2'] = a.reindex(x.g2).values; x = x[(x.t1 >= 0) & (x.t2 >= 0)]
        x['same_tad'] = (x.t1 == x.t2).astype(int)
        if sub is not None: x = x[sub(x)]
        m = smf.ols('r ~ same_tad + C(bin) + C(orientation)', data=x).fit(); return m.params['same_tad'], m.pvalues['same_tad'], len(x), x.same_tad.mean()
    b_all, p_all, n_all, f_same = fit(tn)
    b_ng, p_ng, n_ng, _ = fit(tn, sub=lambda x: ~x.share_any)                      # no shared eQTL: non-genetic component
    others = [fit(o)[0] for o in TADA if o != tn]
    rows.append({'tissue': t, 'tad_map': tn, 'pairs': n_all, 'frac_same_tad': f_same, 'b_same_tad': b_all, 'p': p_all, 'b_same_tad_no_shared_eQTL': b_ng, 'p_no_shared_eQTL': p_ng,
                 'b_same_tad_other_maps_median': np.median(others), 'own_map_rank_among_maps': 1 + sum(o > b_all for o in others), 'n_maps': 1 + len(others)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_cis_test.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
print(f"\nsame-TAD effect > 0 in {(R.b_same_tad > 0).sum()}/{len(R)} (P<0.05 in {(R.p < 0.05).sum()}); median {R.b_same_tad.median():.3f}; without shared eQTL median {R.b_same_tad_no_shared_eQTL.median():.3f}; own map > median of other maps in {(R.b_same_tad > R.b_same_tad_other_maps_median).sum()}/{len(R)}")
