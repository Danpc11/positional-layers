"""eQTL law across tissues (ST4A): Pearson r and the slope of observed on predicted coupling with tissue intercepts, on both observed
scales and for all shared variants versus the lead variant only, with 95% intervals from 1,000 resamples of 10-Mb genomic blocks."""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR
d = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
d['block'] = d.tissue + ':' + G.Chromosome.reindex(d.g1).astype(str).values + ':' + (G.Start.reindex(d.g1).fillna(0).values // 1e7).astype(int).astype(str)
def slope(df, y, x): X = np.column_stack([df[x].values, pd.get_dummies(df.tissue).values.astype(float)]); return np.linalg.lstsq(X, df[y].values, rcond=None)[0][0]
blocks = d.block.unique(); grp = {b: np.where(d.block.values == b)[0] for b in blocks}; rng = np.random.default_rng(7); res = []
for y in ('obs_r_int_pc15', 'obs_r_int_raw'):
    for x in ('pred_r', 'pred_r_lead_only'):
        bs, rs = [], []
        for _ in range(1000):
            s = d.iloc[np.concatenate([grp[b] for b in rng.choice(blocks, len(blocks))])]; bs.append(slope(s, y, x)); rs.append(s[[x, y]].corr().iloc[0, 1])
        res.append({'observed_scale': y, 'prediction': x, 'pairs': len(d), 'pearson_r': d[[x, y]].corr().iloc[0, 1], 'r_ci': '%.2f-%.2f' % tuple(np.percentile(rs, [2.5, 97.5])),
                    'slope': slope(d, y, x), 'slope_ci': '%.2f-%.2f' % tuple(np.percentile(bs, [2.5, 97.5]))})
R = pd.DataFrame(res); R.to_csv(OUTDIR + 'eqtl_law_v2_summary.csv', index=False); print(R.round(3).to_string(index=False))
