"""Review point 9: genomic-block bootstrap (10-Mb blocks) of the same-TAD effect in each tissue, with the model of Fig. 4b:
r ~ same_tad + distance bin + orientation, on adjacent pairs with both genes in a TAD."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import os, glob, numpy as np, pandas as pd, pyannotables as pa
from lib_tad import MATCH, tads, assign
TADA = {n: assign(tads(n)) for n in sorted(set(MATCH.values()))}
BINS = [0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]; G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
rng = np.random.default_rng(5); rows = []
for t, tn in MATCH.items():
    f = [x for x in glob.glob(OUTDIR + 'pairs_*.csv.gz') if os.path.basename(x)[6:-7].lower().replace('-', '_') == t][0]
    d = pd.read_csv(f); d = d[d.dist > 0].copy(); a = TADA[tn]; d['t1'] = a.reindex(d.g1).values; d['t2'] = a.reindex(d.g2).values
    d = d[(d.t1 >= 0) & (d.t2 >= 0)].reset_index(drop=True); d['same'] = (d.t1 == d.t2).astype(float)
    X = np.column_stack([d.same.values, pd.get_dummies(pd.cut(d.dist, BINS).astype(str)).values.astype(float), pd.get_dummies(d.orientation, drop_first=True).values.astype(float)])
    y = d.r.values; est = np.linalg.lstsq(X, y, rcond=None)[0][0]
    blk = (G.Chromosome.reindex(d.g1).astype(str).values + ':' + (G.Start.reindex(d.g1).fillna(0).values // 1e7).astype(int).astype(str))
    ub, inv = np.unique(blk, return_inverse=True); grp = [np.where(inv == k)[0] for k in range(len(ub))]
    bs = np.array([np.linalg.lstsq(X[ix], y[ix], rcond=None)[0][0] for ix in (np.concatenate([grp[k] for k in rng.integers(0, len(grp), len(grp))]) for _ in range(500))])
    rows.append({'tissue': t, 'pairs': len(d), 'blocks': len(ub), 'b_same_tad': est, 'ci_low': np.percentile(bs, 2.5), 'ci_high': np.percentile(bs, 97.5), 'boot_se': bs.std(),
                 'p_boot_normal': 2 * __import__('scipy.stats', fromlist=['norm']).norm.sf(abs(est / bs.std()))})
    print(rows[-1]['tissue'], round(est, 4), round(rows[-1]['ci_low'], 4), round(rows[-1]['ci_high'], 4), flush=True)
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_block_bootstrap.csv', index=False)
print(f"CI excludes 0 in {(R.ci_low > 0).sum()}/{len(R)} tissues; median effect {R.b_same_tad.median():.3f}")
