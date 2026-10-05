"""Approximate colocalisation of adjacent genes with SuSiE credible sets: P(same causal variant) = max over credible-set
pairs of sum_v pip1(v)*pip2(v); direction from the shared variant with the highest pip product (afc sign)."""
import os, sys
from poslayers.config import OUTDIR
import sys, os, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
norm = lambda s: s.lower().replace('-', '_')
PAIRS = pd.concat([pd.read_csv(f, usecols=['g1', 'g2']) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')]).drop_duplicates(); genes = set(PAIRS.g1) | set(PAIRS.g2)
atlas = {norm(os.path.basename(f)[6:-7]) for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
for f in sys.argv[1:]:
    t = norm(os.path.basename(f).replace('_v11_eQTLs_SuSiE_summary.parquet', '')); out = OUTDIR + f'coloc_{t}.csv.gz'
    if t not in atlas or os.path.exists(out): continue
    E = pq.read_table(f, columns=['phenotype_id', 'variant_id', 'pip', 'cs_id', 'afc']).to_pandas(); E['g'] = E.phenotype_id.str.split('.').str[0]; E = E[E.g.isin(genes)]
    CS = {g: [(dict(zip(c.variant_id, c.pip)), dict(zip(c.variant_id, np.sign(c.afc)))) for _, c in d.groupby('cs_id')] for g, d in E.groupby('g')}
    rows = []
    for g1, g2 in zip(PAIRS.g1.values, PAIRS.g2.values):
        a, b = CS.get(g1), CS.get(g2)
        if a is None or b is None: rows.append((g1, g2, a is not None, b is not None, 0.0, np.nan)); continue
        best, sgn = 0.0, np.nan
        for p1, s1 in a:
            for p2, s2 in b:
                sh = p1.keys() & p2.keys()
                if not sh: continue
                pr = sum(p1[v] * p2[v] for v in sh)
                if pr > best:
                    best = pr; vs = sorted(sh, key=lambda v: -p1[v] * p2[v]); sgn = np.nan
                    for v in vs:
                        if np.isfinite(s1.get(v, np.nan)) and np.isfinite(s2.get(v, np.nan)): sgn = s1[v] * s2[v]; break
        rows.append((g1, g2, True, True, best, sgn))
    S = pd.DataFrame(rows, columns=['g1', 'g2', 'fm1', 'fm2', 'p_coloc', 'direction']); S.to_csv(out, index=False)
    print(t, f'fine-mapped pairs {int((S.fm1 & S.fm2).sum())} | coloc>=0.5 {int((S.p_coloc >= 0.5).sum())}', flush=True)
