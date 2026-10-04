import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
"""Per tissue: for every adjacent gene pair of the atlas, are both genes eGenes, and do they share a significant eQTL variant
(any direction; same direction of effect)? Reads only the needed parquet columns."""
import sys, os, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
PAIRS = pd.concat([pd.read_csv(f, usecols=['g1', 'g2']) for f in glob.glob(DATA + 'pairs_*.csv.gz')]).drop_duplicates()
genes = set(PAIRS.g1) | set(PAIRS.g2)
norm = lambda s: s.lower().replace('-', '_')
for f in sys.argv[1:]:
    t = norm(os.path.basename(f).replace('_v11_eQTLs_signif_pairs.parquet', '')); out = fDATA + 'shares_{t}.csv.gz'
    if os.path.exists(out): continue
    pf = pq.ParquetFile(f); sign = {}; nE = 0
    for rg in range(pf.num_row_groups):          # stream row groups to bound memory
        E = pf.read_row_group(rg, columns=['phenotype_id', 'variant_id', 'slope']).to_pandas()
        E['g'] = E.phenotype_id.str.split('.').str[0]; E = E[E.g.isin(genes)]; nE += len(E)
        for g, d in E.groupby('g'):
            sign.setdefault(g, {}).update(zip(d.variant_id.values, np.sign(d.slope.values).astype(np.int8)))
        del E
    E = range(nE)
    rows = []
    for g1, g2 in zip(PAIRS.g1.values, PAIRS.g2.values):
        a, b = sign.get(g1), sign.get(g2)
        if a is None or b is None: rows.append((g1, g2, a is not None, b is not None, False, False, 0)); continue
        sh = a.keys() & b.keys(); same = sum(1 for v in sh if a[v] == b[v])
        rows.append((g1, g2, True, True, len(sh) > 0, same > 0, len(sh)))
    S = pd.DataFrame(rows, columns=['g1', 'g2', 'egene1', 'egene2', 'share_any', 'share_same', 'n_shared'])
    S.to_csv(out, index=False); print(t, f'{len(E):,} eQTL pairs | both eGenes {int((S.egene1 & S.egene2).sum())} | share {int(S.share_any.sum())}', flush=True)
