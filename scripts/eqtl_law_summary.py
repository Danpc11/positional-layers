"""eQTL law across tissues (Supplementary Table 4A), for two inclusion thresholds on the colocalisation score.

eqtl_law.py includes every adjacent pair with a colocalisation score >= 0.1 (the inclusion threshold of this analysis);
the 'colocalised' class of Fig. 3d uses >= 0.5. Both are reported here. For each threshold and observed scale:
  - Pearson r and the slope of observed on predicted coupling with tissue intercepts;
  - 95% intervals from 1,000 bootstrap replicates that resample 10-Mb genomic blocks WITHIN each tissue (stratified);
  - sign concordance: among pairs with predicted genetic correlation < 0, the number and fraction with observed coupling < 0.
Also written:
  - eqtl_law_strata.csv: mean observed coupling in eight strata of the prediction (Fig. 3f), with 95% intervals from the
    same tissue-stratified block bootstrap (strata boundaries fixed on the full data);
  - eqtl_law_block_sensitivity.csv: r and slope with intervals from (i) tissue-stratified blocks of 5, 10 and 20 Mb and
    (ii) blocks coordinated across tissues, in which the same genomic blocks are drawn for every tissue so that shared
    regions stay together. Donors shared between tissues are not resampled (genotypes are not available), so all
    intervals are conditional on the observed donors and describe genomic variation only.
"""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR

d0 = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
def blocks(size):
    st = G.Start.reindex(d0.g1).values; ch = G.Chromosome.reindex(d0.g1).astype(str).values
    return np.array([f'{c}:{int(x // size)}' if np.isfinite(x) else f'unplaced:{g}' for c, x, g in zip(ch, st, d0.g1)], dtype=object)
d0['block'] = blocks(1e7)


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

# ---- Fig. 3f: stratum means with tissue-stratified block-bootstrap intervals
d = d0[d0.coloc_score >= 0.1].reset_index(drop=True); d['stratum'] = pd.qcut(d.pred_r, 8, labels=False, duplicates='drop')
strata = {t: [g.index.values for _, g in sub.groupby('block')] for t, sub in d.groupby('tissue')}; rng = np.random.default_rng(11)
reps = [np.concatenate([blk[k] for t, blk in strata.items() for k in rng.integers(0, len(blk), len(blk))]) for _ in range(1000)]
est = d.groupby('stratum').agg(pred=('pred_r', 'mean'), obs=('obs_r_int_pc15', 'mean'), pairs=('pred_r', 'size'))
bs = np.array([d.iloc[ix].groupby('stratum').obs_r_int_pc15.mean().reindex(est.index).values for ix in reps])
est['ci_low'], est['ci_high'] = np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)
est.reset_index().to_csv(OUTDIR + 'eqtl_law_strata.csv', index=False); print(est.round(3).to_string())
# ---- sensitivity: block size, and blocks coordinated across tissues
sens = []
for size in (5e6, 1e7, 2e7):
    dd = d.assign(block=blocks(size)[d0.coloc_score.values >= 0.1])
    st_ = {t: [g.index.values for _, g in sub.groupby('block')] for t, sub in dd.groupby('tissue')}; r_ = np.random.default_rng(13)
    rp = [np.concatenate([blk[k] for t, blk in st_.items() for k in r_.integers(0, len(blk), len(blk))]) for _ in range(1000)]
    for scheme, R_ in (('within tissue', rp), ('coordinated across tissues', None)):
        if R_ is None:
            ub = {b: g.index.values for b, g in dd.groupby('block')}; keys = list(ub); r2 = np.random.default_rng(17)
            R_ = [np.concatenate([ub[keys[k]] for k in r2.integers(0, len(keys), len(keys))]) for _ in range(1000)]
        rs = [np.corrcoef(dd.pred_r.values[ix], dd.obs_r_int_pc15.values[ix])[0, 1] for ix in R_]; sl = [slope(dd.iloc[ix], 'obs_r_int_pc15', 'pred_r') for ix in R_]
        sens.append({'block_mb': size / 1e6, 'scheme': scheme, 'r': np.corrcoef(dd.pred_r, dd.obs_r_int_pc15)[0, 1], 'r_ci_low': np.percentile(rs, 2.5), 'r_ci_high': np.percentile(rs, 97.5),
                     'slope': slope(dd, 'obs_r_int_pc15', 'pred_r'), 'slope_ci_low': np.percentile(sl, 2.5), 'slope_ci_high': np.percentile(sl, 97.5)})
S = pd.DataFrame(sens); S.to_csv(OUTDIR + 'eqtl_law_block_sensitivity.csv', index=False); print(S.round(3).to_string(index=False))
