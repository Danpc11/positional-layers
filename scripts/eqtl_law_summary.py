"""eQTL law across tissues (Supplementary Table 4A), for two inclusion thresholds on the colocalisation score.

eqtl_law.py includes every adjacent pair with a colocalisation score >= 0.1 (the inclusion threshold of this analysis);
the 'colocalised' class of Fig. 3d uses >= 0.5. Both are reported here. For each threshold and observed scale:
  - Pearson r and the slope of observed on predicted coupling with tissue intercepts;
  - 95% intervals from 1,000 bootstrap replicates that resample 10-Mb genomic blocks WITHIN each tissue (stratified);
  - sign concordance: among pairs with predicted genetic correlation < 0, the number and fraction with observed coupling < 0.
"""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR

d0 = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
d0['block'] = G.Chromosome.reindex(d0.g1).astype(str).values + ':' + (G.Start.reindex(d0.g1).fillna(0).values // 1e7).astype(int).astype(str)


def slope(df, y, x):
    X = np.column_stack([df[x].values, pd.get_dummies(df.tissue).values.astype(float)])
    return np.linalg.lstsq(X, df[y].values, rcond=None)[0][0]


res = []
for thr in (0.1, 0.5):
    d = d0[d0.coloc_score >= thr].reset_index(drop=True)
    strata = {t: [g.index.values for _, g in sub.groupby('block')] for t, sub in d.groupby('tissue')}
    rng = np.random.default_rng(7)
    reps = [np.concatenate([blk[k] for t, blk in strata.items() for k in rng.integers(0, len(blk), len(blk))]) for _ in range(1000)]
    for y in ('obs_r_int_pc15', 'obs_r_int_raw'):
        for x in ('pred_r', 'pred_r_lead_only'):
            bs = [slope(d.iloc[ix], y, x) for ix in reps]; rs = [np.corrcoef(d[x].values[ix], d[y].values[ix])[0, 1] for ix in reps]
            neg = d[d[x] < 0]; k = int((neg[y] < 0).sum())
            res.append({'score_threshold': thr, 'observed_scale': y, 'prediction': x, 'pairs': len(d), 'tissues': d.tissue.nunique(),
                        'pearson_r': np.corrcoef(d[x], d[y])[0, 1], 'r_ci_low': np.percentile(rs, 2.5), 'r_ci_high': np.percentile(rs, 97.5),
                        'slope_tissue_intercepts': slope(d, y, x), 'slope_ci_low': np.percentile(bs, 2.5), 'slope_ci_high': np.percentile(bs, 97.5),
                        'slope_pooled': np.polyfit(d[x], d[y], 1)[0],
                        'pred_negative_pairs': len(neg), 'pred_negative_obs_negative': k, 'sign_concordance': k / len(neg) if len(neg) else np.nan})
R = pd.DataFrame(res); R.to_csv(OUTDIR + 'eqtl_law_v2_summary.csv', index=False)
print(R[R.prediction == 'pred_r'].round(3).to_string(index=False))
