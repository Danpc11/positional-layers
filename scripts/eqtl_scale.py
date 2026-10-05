"""Scale of the eQTL relation (Supplementary Table 4C): does the slope of observed on predicted coupling depend on how the
observed quantity is normalised and on how many expression PCs are removed?

GTEx slopes are estimated on inverse-normal (INT) expression with unit variance, so the predicted genetic contribution
sum_v w_v 2p(1-p) b1 b2 is a covariance on that scale. Removing K PCs lowers each gene's variance; re-standardising (a
correlation) then inflates the observed value relative to the prediction. For K = 0, 5, 10, 15, 25 and 40 we report:
  corr   : correlation after removing K PCs (re-standardised; what eqtl_law.py reports for K = 15)
  cov_sd0: covariance after removing K PCs divided by the ORIGINAL standard deviations (common scale with the prediction)
with Pearson r and the tissue-intercept slope. The LD structure of marginal slopes is not modelled; the comparison
separates normalisation from the remaining gap, it does not calibrate the prediction.

Input : OUTDIR/eqtl_law_v2_pairs.csv (eqtl_law.py), DATA/gtex_eqtl/<Tissue>_v11_normalized_expression_bed.gz
Output: OUTDIR/eqtl_scale_pairs.csv (per pair and K), OUTDIR/eqtl_scale_sensitivity.csv (summary)
"""
import os, numpy as np, pandas as pd
from poslayers.config import DATA, OUTDIR, upsert_csv, RESUME

NAME = {'nerve_tibial': 'Nerve_Tibial', 'thyroid': 'Thyroid', 'cells_cultured_fibroblasts': 'Cells_Cultured_fibroblasts', 'artery_tibial': 'Artery_Tibial',
        'whole_blood': 'Whole_Blood', 'testis': 'Testis', 'skin_sun_exposed_lower_leg': 'Skin_Sun_Exposed_Lower_leg', 'esophagus_mucosa': 'Esophagus_Mucosa',
        'adipose_subcutaneous': 'Adipose_Subcutaneous', 'lung': 'Lung'}
KS = [0, 5, 10, 15, 25, 40]
P = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv'); OUT = OUTDIR + 'eqtl_scale_pairs.csv'
done = set(pd.read_csv(OUT).tissue) if (RESUME and os.path.exists(OUT)) else set()
for t, N in NAME.items():
    if t in done or t not in set(P.tissue): continue
    bed = pd.read_csv(DATA + f'gtex_eqtl/{N}_v11_normalized_expression_bed.gz', sep='\t'); bed['gid'] = bed.iloc[:, 3].str.split('.').str[0]
    X = bed.drop_duplicates('gid').set_index('gid').iloc[:, 4:].astype(float); idx = {g: i for i, g in enumerate(X.index)}
    Xc = X.values - X.values.mean(1, keepdims=True); sd0 = Xc.std(1); U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    pr = P[P.tissue == t]; pr = pr[pr.g1.isin(idx) & pr.g2.isin(idx)]; i1 = pr.g1.map(idx).values; i2 = pr.g2.map(idx).values
    rows = []
    for K in KS:
        R = Xc - (U[:, :K] * s[:K]) @ Vt[:K] if K else Xc
        cov = np.mean((R[i1] - R[i1].mean(1, keepdims=True)) * (R[i2] - R[i2].mean(1, keepdims=True)), 1); sd = R.std(1)
        rows.append(pd.DataFrame({'tissue': t, 'g1': pr.g1.values, 'g2': pr.g2.values, 'K': K, 'pred_r': pr.pred_r.values, 'coloc_score': pr.coloc_score.values,
                                  'corr': cov / (sd[i1] * sd[i2]), 'cov_sd0': cov / (sd0[i1] * sd0[i2])}))
    upsert_csv(pd.concat(rows), OUT, ['tissue', 'g1', 'g2', 'K']); print(t, len(pr), 'pairs', flush=True)
    del bed, X, Xc, U, Vt

D = pd.read_csv(OUT)
def slope(df, y):
    X = np.column_stack([df.pred_r.values, pd.get_dummies(df.tissue).values.astype(float)]); return np.linalg.lstsq(X, df[y].values, rcond=None)[0][0]
S = []
for (K, thr), _ in pd.DataFrame([(k, t) for k in KS for t in (0.1, 0.5)]).groupby([0, 1]):
    d = D[(D.K == K) & (D.coloc_score >= thr)]
    for y in ('corr', 'cov_sd0'):
        S.append({'K_PCs_removed': K, 'score_threshold': thr, 'measure': y, 'pairs': len(d), 'pearson_r': np.corrcoef(d.pred_r, d[y])[0, 1],
                  'slope_tissue_intercepts': slope(d, y), 'median_observed': d[y].median()})
S = pd.DataFrame(S); S.to_csv(OUTDIR + 'eqtl_scale_sensitivity.csv', index=False)
print(S[S.score_threshold == 0.1].round(3).to_string(index=False))
