"""Friction 1: is the TAD effect explained by 3D contact? All gene pairs within 2 Mb (GM12878 / IMR-90 Hi-C, GTEx EBV / fibroblasts)."""
from poslayers.config import OUTDIR
import pandas as pd
import statsmodels.formula.api as smf
from lib_tad import tads, assign
rows = []
for cell, tadname in [('GM12878', 'GM12878_lymphoblastoid_Lieberman'), ('IMR90', 'IMR90_fetalLungFibroblast_Lieberman')]:
    H = pd.read_csv(OUTDIR + f'hic_coupling_{cell}.csv.gz'); a = assign(tads(tadname))
    H['t1'] = a.reindex(H.g1).values; H['t2'] = a.reindex(H.g2).values; H = H[(H.t1 >= 0) & (H.t2 >= 0)].copy(); H['same_tad'] = (H.t1 == H.t2).astype(int)
    H['db'] = pd.cut(H.tss_distance, [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6]).astype(str)
    for b, d in H.groupby('db'):
        if d.same_tad.sum() < 50 or (1 - d.same_tad).sum() < 50: continue
        m0 = smf.ols('r ~ same_tad + log_d', data=d).fit(); m1 = smf.ols('r ~ same_tad + log_d + log_oe', data=d).fit()
        rows.append({'cell': cell, 'distance': b, 'pairs': len(d), 'frac_same_tad': d.same_tad.mean(), 'r_same': d.r[d.same_tad == 1].mean(), 'r_diff': d.r[d.same_tad == 0].mean(),
                     'TAD_effect': m0.params['same_tad'], 'p_TAD': m0.pvalues['same_tad'], 'TAD_effect_given_contact': m1.params['same_tad'], 'p_TAD_given_contact': m1.pvalues['same_tad'],
                     'contact_effect_given_TAD': m1.params['log_oe'], 'p_contact_given_TAD': m1.pvalues['log_oe']})
    m0 = smf.ols('r ~ same_tad + bs(log_d, df=5)', data=H).fit(); m1 = smf.ols('r ~ same_tad + bs(log_d, df=5) + log_oe', data=H).fit()
    print(f"{cell} all distances: TAD effect {m0.params['same_tad']:.4f} (t {m0.tvalues['same_tad']:.1f}) -> given contact {m1.params['same_tad']:.4f} (t {m1.tvalues['same_tad']:.1f}); contact given TAD t {m1.tvalues['log_oe']:.1f}")
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'friction1_tad_vs_contact.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
# short distances (< 20 kb, unresolved by 25-kb Hi-C): adjacent pairs in GTEx, same vs different TAD
P = []
for tissue, tadname in [('cells_ebv-transformed_lymphocytes', 'GM12878_lymphoblastoid_Lieberman'), ('cells_cultured_fibroblasts', 'IMR90_fetalLungFibroblast_Lieberman')]:
    d = pd.read_csv(OUTDIR + f'pairs_{tissue}.csv.gz'); a = assign(tads(tadname)); d['t1'] = a.reindex(d.g1).values; d['t2'] = a.reindex(d.g2).values
    d = d[(d.t1 >= 0) & (d.t2 >= 0) & (d.dist > 0)]; d['same'] = d.t1 == d.t2; d['db'] = pd.cut(d.dist, [0, 5e3, 2e4, 5e4, 2e5]).astype(str)
    for b, x in d.groupby('db'): P.append({'tissue': tissue, 'distance': b, 'n_same': int(x.same.sum()), 'n_diff': int((~x.same).sum()), 'r_same': x.r[x.same].mean(), 'r_diff': x.r[~x.same].mean()})
print(pd.DataFrame(P).round(3).to_string(index=False)); pd.DataFrame(P).to_csv(OUTDIR + 'friction1_short_adjacent.csv', index=False)
