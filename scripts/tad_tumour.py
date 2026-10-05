from poslayers.config import DATA, OUTDIR
import glob
import numpy as np
import pandas as pd
import pyannotables as pa
import statsmodels.api as sm
from lib_tumour import CHR, coding, prepare, trunc

B = pd.read_csv(glob.glob(DATA + 'TAD-full/*/data/boundariesByStability/100kbBookendBoundaries_mainText/100kbBookendBoundaries_byStability.bed')[0], sep='\t')
B['chr'] = B.chr.str.replace('chr', ''); B['mid'] = (B['loc'] + B['loc2']) / 2
G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()]; mid37 = ((G37.Start + G37.End) / 2); chr37 = G37.Chromosome.astype(str)
def pair_classes(gids, chrs):
    i = np.where(chrs[:-1] == chrs[1:])[0]; a = gids[i]; b = gids[i + 1]
    m1 = mid37.reindex(a).values; m2 = mid37.reindex(b).values; c1 = chr37.reindex(a).values
    cls = np.full(len(i), 'na', dtype=object)
    for c, bb in B.groupby('chr'):
        sel = np.where((c1 == c) & np.isfinite(m1) & np.isfinite(m2))[0]
        if not len(sel): continue
        lo = np.minimum(m1[sel], m2[sel]); hi = np.maximum(m1[sel], m2[sel]); bm = bb['mid'].values; st = bb.stability_percentile.values
        for k, (l, h) in zip(sel, zip(lo, hi)):
            inside = (bm > l) & (bm < h)
            cls[k] = 'within' if not inside.any() else ('stable' if st[inside].max() >= 0.75 else 'weak')
    return i, cls
def score(D, chrs, i, mask):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9); near = (Z[i[mask]] * Z[i[mask] + 1]).mean(0)
    far = np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in (20, 30)], 0)
    return near - far
rows = []
for c, gm in [('BLCA', ['STAG2']), ('UCEC', ['CTCF'])]:
    D, chrs, samp, cov, gids = prepare(c, return_gids=True); i, cls = pair_classes(gids, chrs)
    print(c, pd.Series(cls).value_counts().to_dict(), flush=True)
    S = {k: score(D, chrs, i, cls == k) for k in ('within', 'stable')}
    for mclass, srcm in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            keep = cov.log_tmb.values <= np.quantile(cov.log_tmb.values, hq); m = pd.Index(samp).isin(set(srcm[srcm.gene.isin(gm)].Sample_ID)).astype(float)[keep]
            Xc = cov[keep]; Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join((Xc - Xc.mean()) / Xc.std()))
            res = {'cohort': c, 'gene': gm[0], 'mutation_class': mclass, 'tmb_q': hq, 'n_mut': int(m.sum())}
            for k in ('within', 'stable'):
                y = S[k][keep]; f = sm.OLS(y, Xd).fit(cov_type='HC3'); base = y[m == 0].mean()
                res[f'{k}_wt'] = base; res[f'{k}_effect_pct'] = 100 * f.params['mutant'] / abs(base); res[f'{k}_p'] = f.pvalues['mutant']
            dd = S['within'][keep] - S['stable'][keep]; f = sm.OLS(dd, Xd).fit(cov_type='HC3'); res['within_minus_stable_effect'] = f.params['mutant']; res['p_difference'] = f.pvalues['mutant']
            rows.append(res)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_tumour_results.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
