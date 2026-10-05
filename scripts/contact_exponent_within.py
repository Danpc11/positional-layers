"""Test of the source theory, cis layer: one contact exponent for two independent sources of contact variation.

Theory: if a gene's loading on a regulatory element scales as contact^k, coupling scales as contact^k whether contact
varies because of genomic distance or, at a FIXED distance, because of 3D folding. k_within is estimated only from
variation of observed/expected contact inside distance bins (quintiles of O/E within each bin, bin fixed effects);
k_across only from the change of mean contact between distance bins. The theory predicts k_within = k_across, so the
decay of coupling with distance follows from the Hi-C decay of contact with no fitted constant beyond the scale.
95% intervals from 300 resamples of 10-Mb genomic blocks. Pairs with zero contact are excluded, as in boot_hic.py.
Output: OUTDIR/contact_exponent_within.csv (Supplementary Note 1, section 6)
"""
import numpy as np, pandas as pd
from poslayers.config import OUTDIR
from lib_boot import blocks_for

BINS = [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6]
def exponents(d):
    d = d.assign(db=pd.cut(d.tss_distance, BINS, labels=False))
    g = d.groupby(['db', pd.qcut(d.groupby('db').oe.rank(pct=True), 5, labels=False)], observed=True).agg(r=('r', 'mean'), c=('contact_KR', 'mean')).reset_index()
    g = g[g.r > 0]; X = np.column_stack([np.log(g.c), pd.get_dummies(g.db).values.astype(float)])
    k_within = np.linalg.lstsq(X, np.log(g.r), rcond=None)[0][0]
    a = d.groupby('db').agg(r=('r', 'mean'), c=('contact_KR', 'mean')); a = a[a.r > 0]
    k_across = np.polyfit(np.log(a.c), np.log(a.r), 1)[0]
    return k_within, k_across
if __name__ == '__main__':
    rows = []
    for cell in ('GM12878', 'IMR90'):
        H = pd.read_csv(OUTDIR + f'hic_coupling_{cell}.csv.gz', dtype={'chr': str}); H = H[(H.contact_KR > 0) & H.tss_distance.between(BINS[0], BINS[-1])].reset_index(drop=True)
        kw, ka = exponents(H); blk = blocks_for(H.g1.values); grp = [np.where(blk == b)[0] for b in np.unique(blk)]; rng = np.random.default_rng(9); bs = []
        for _ in range(300):
            bs.append(exponents(H.iloc[np.concatenate([grp[i] for i in rng.integers(0, len(grp), len(grp))])]))
        bs = np.array(bs); diff = bs[:, 0] - bs[:, 1]
        rows.append({'cell': cell, 'pairs': len(H), 'k_within_distance': kw, 'kw_ci_low': np.percentile(bs[:, 0], 2.5), 'kw_ci_high': np.percentile(bs[:, 0], 97.5),
                     'k_across_distance': ka, 'ka_ci_low': np.percentile(bs[:, 1], 2.5), 'ka_ci_high': np.percentile(bs[:, 1], 97.5),
                     'difference': kw - ka, 'diff_ci_low': np.percentile(diff, 2.5), 'diff_ci_high': np.percentile(diff, 97.5)})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUTDIR + 'contact_exponent_within.csv', index=False)
