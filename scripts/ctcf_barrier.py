"""Coupling across convergent CTCF loop anchors at equal distance and contact (Fig. 4e; Extended Data Fig. 7c,d; Supplementary Table 9B,C).
Rao et al. 2014 HiCCUPS loops with CTCF motifs; GTEx coupling of the matching cell type. Usage: python ctcf_barrier.py GM12878|IMR90

Convergent CTCF loops (forward motif at the left anchor, reverse at the right) mark where extrusion stops. For each gene
pair (GRCh37 TSSs):  same_loop = both TSSs inside the span [x1, y2] of at least one convergent loop;
                     n_cross   = number of convergent loops with exactly one of the two TSSs inside their span.
Prediction fixed in advance: at equal distance AND equal Hi-C contact, coupling is higher for same_loop and lower with
each crossed loop. Model: r ~ log contact + distance-bin fixed effects + same_loop + min(n_cross, 3).
95% intervals from 200 resamples of 10-Mb genomic blocks.
Outputs: OUTDIR/ctcf_barrier_<cell>.csv (coefficients), OUTDIR/ctcf_barrier_means_<cell>.csv (adjusted means by crossed loops)
"""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA, OUTDIR
import sys
from lib_boot import blocks_for
CELL = sys.argv[1] if len(sys.argv) > 1 else 'GM12878'
LOOPS = {'GM12878': 'GSE63525_GM12878_primary_replicate_HiCCUPS_looplist_with_motifs_txt.gz', 'IMR90': 'GSE63525_IMR90_HiCCUPS_looplist_with_motifs_txt.gz'}[CELL]

G = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G = G[~G.index.duplicated()]
TSS = pd.Series(np.where(G.Strand.astype(str).isin(['1', '+']), G.Start, G.End), index=G.index)
L = pd.read_csv(DATA + 'hic/' + LOOPS, sep='\t', dtype={'chr1': str}); L.columns = L.columns.str.strip()
L = L[(L.orientation_1 == 'p') & (L.orientation_2 == 'n') & (L.chr1 == L['chr2'].astype(str))]
H = pd.read_csv(OUTDIR + f'hic_coupling_{CELL}.csv.gz', dtype={'chr': str})
H = H[(H.contact_KR > 0) & H.tss_distance.between(25e3, 2e6) & H.chr.isin([str(i) for i in range(1, 23)])].reset_index(drop=True)
t1, t2 = TSS.reindex(H.g1).values, TSS.reindex(H.g2).values; H = H[np.isfinite(t1) & np.isfinite(t2)].reset_index(drop=True)
a, b = np.minimum(TSS.reindex(H.g1).values, TSS.reindex(H.g2).values), np.maximum(TSS.reindex(H.g1).values, TSS.reindex(H.g2).values)
same, ncross = np.zeros(len(H), bool), np.zeros(len(H), int)
for ch, l in L.groupby('chr1'):
    idx = np.where(H.chr.values == ch)[0]
    if not len(idx): continue
    lo, hi = l.x1.values, l.y2.values
    for s in range(0, len(idx), 20000):
        j = idx[s:s + 20000]; ia = (a[j, None] >= lo) & (a[j, None] <= hi); ib = (b[j, None] >= lo) & (b[j, None] <= hi)
        same[j] = (ia & ib).any(1); ncross[j] = (ia ^ ib).sum(1)
H['same_loop'] = same.astype(float); H['ncross'] = np.minimum(ncross, 3); H['logc'] = np.log(H.contact_KR)
H['db'] = pd.cut(H.tss_distance, [25e3, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6], labels=False)
print(f'{len(H):,} pairs; in a shared convergent loop {H.same_loop.mean():.1%}; crossing >= 1 loop {np.mean(H.ncross > 0):.1%}', flush=True)
def fit(d, with_contact=True):
    X = np.column_stack(([d.logc.values] if with_contact else []) + [d.same_loop.values, d.ncross.values, pd.get_dummies(d.db).values.astype(float)])
    beta = np.linalg.lstsq(X, d.r.values, rcond=None)[0]; return beta[:3] if with_contact else np.r_[np.nan, beta[:2]]
blk = blocks_for(H.g1.values); grp = [np.where(blk == x)[0] for x in np.unique(blk)]; rng = np.random.default_rng(8); rows = []
for wc in (False, True):
    est = fit(H, wc); bs = np.array([fit(H.iloc[np.concatenate([grp[i] for i in rng.integers(0, len(grp), len(grp))])], wc) for _ in range(200)])
    for k, nm in ((1, 'same_loop'), (2, 'per crossed loop')):
        rows.append({'adjusted_for_contact': wc, 'term': nm, 'coef': est[k], 'ci_low': np.percentile(bs[:, k], 2.5), 'ci_high': np.percentile(bs[:, k], 97.5)})
    if wc: rows.append({'adjusted_for_contact': wc, 'term': 'log contact', 'coef': est[0], 'ci_low': np.percentile(bs[:, 0], 2.5), 'ci_high': np.percentile(bs[:, 0], 97.5)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + f'ctcf_barrier_{CELL}.csv', index=False); print(R.round(4).to_string(index=False))
print('mean coupling by crossings (all distances):', H.groupby('ncross').r.mean().round(4).to_dict(), '| same loop vs not:', H.groupby('same_loop').r.mean().round(4).to_dict())

# ---- adjusted means by number of crossed loops (Fig. 4e): residual within distance-bin x contact-quintile strata
H['cq'] = H.groupby('db').logc.rank(pct=True).mul(5).clip(upper=4.999).astype(int); H['stratum'] = H.db.astype(str) + '_' + H.cq.astype(str)
def adjusted(d): return (d.r - d.groupby('stratum').r.transform('mean') + d.r.mean()).groupby(d.ncross).mean()
rng_m = np.random.default_rng(23); est = adjusted(H); bs = pd.DataFrame([adjusted(H.iloc[np.concatenate([grp[i] for i in rng_m.integers(0, len(grp), len(grp))])]) for _ in range(200)])
rel = bs.div(bs[0], axis=0)                     # ratio to 'no loop crossed' within each replicate: includes denominator uncertainty
pd.DataFrame({'cell': CELL, 'crossed_loops': est.index, 'pairs': H.ncross.value_counts().reindex(est.index).values, 'coupling': est.values,
              'ci_low': bs.quantile(0.025).values, 'ci_high': bs.quantile(0.975).values, 'relative': (est / est[0]).values,
              'relative_ci_low': rel.quantile(0.025).values, 'relative_ci_high': rel.quantile(0.975).values}).to_csv(OUTDIR + f'ctcf_barrier_means_{CELL}.csv', index=False)
