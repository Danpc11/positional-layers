"""eQTL law on a common scale, with every shared variant's own effect.
Observed coupling: correlation of the two genes in GTEx's normalized (inverse-normal) expression, the scale on which GTEx slopes are
estimated, after removing 15 expression PCs as a stand-in for the hidden factors in the eQTL model (raw INT also kept).
Predicted genetic correlation: sum over pairs of credible sets, and over each shared variant v, of
    PIP_1(v) * PIP_2(v) * 2 p_v (1 - p_v) * beta_1(v) * beta_2(v)
using each variant's own GTEx slope for each gene. The overlap max over credible-set pairs of sum_v PIP_1 PIP_2 is reported as a
colocalisation SCORE, not a posterior probability of a shared causal variant."""
import os
from poslayers.config import DATA, OUTDIR, upsert_csv
import os
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
U = DATA + 'gtex_eqtl/'; OUT = OUTDIR + 'eqtl_law_v2_pairs.csv'
done = set(pd.read_csv(OUT).tissue) if os.path.exists(OUT) else set()
NAME = {'nerve_tibial': 'Nerve_Tibial', 'thyroid': 'Thyroid', 'cells_cultured_fibroblasts': 'Cells_Cultured_fibroblasts', 'artery_tibial': 'Artery_Tibial',
        'whole_blood': 'Whole_Blood', 'testis': 'Testis', 'skin_sun_exposed_lower_leg': 'Skin_Sun_Exposed_Lower_leg', 'esophagus_mucosa': 'Esophagus_Mucosa',
        'adipose_subcutaneous': 'Adipose_Subcutaneous', 'lung': 'Lung'}
for t, N in NAME.items():
    if t in done: continue
    P = pd.read_csv(OUTDIR + f'pred_pairs_{t}.csv.gz'); P = P[(P.dist > 0) & (P.p_coloc >= 0.1)]; genes = set(P.g1) | set(P.g2)
    # observed coupling on GTEx's normalized expression
    bed = pd.read_csv(U + f'{N}_v11_normalized_expression_bed.gz', sep='\t'); bed['gid'] = bed.iloc[:, 3].str.split('.').str[0]
    X = bed.drop_duplicates('gid').set_index('gid').iloc[:, 4:].astype(float)
    Xc = X.values - X.values.mean(1, keepdims=True)
    U_, s_, Vt = np.linalg.svd(Xc, full_matrices=False); Xp = Xc - (U_[:, :15] * s_[:15]) @ Vt[:15]
    idx = {g: i for i, g in enumerate(X.index)}
    def corr(M, a, b):
        x, y = M[a] - M[a].mean(), M[b] - M[b].mean(); return float(x @ y / np.sqrt((x @ x) * (y @ y)))
    # credible sets and slopes
    E = pq.read_table(U + f'{N}_v11_eQTLs_SuSiE_summary.parquet', columns=['phenotype_id', 'variant_id', 'pip', 'cs_id']).to_pandas()
    E['gid'] = E.phenotype_id.str.split('.').str[0]; E = E[E.gid.isin(genes)]
    CS = {g: [dict(zip(c.variant_id, c.pip)) for _, c in d.groupby('cs_id')] for g, d in E.groupby('gid')}
    need = set(E.variant_id); pf = pq.ParquetFile(U + f'{N}_v11_eQTLs_signif_pairs.parquet'); slope, af = {}, {}
    for rg in range(pf.num_row_groups):
        S = pf.read_row_group(rg, columns=['phenotype_id', 'variant_id', 'slope', 'af']).to_pandas(); S = S[S.variant_id.isin(need)]
        S['gid'] = S.phenotype_id.str.split('.').str[0]; S = S[S.gid.isin(genes)]
        slope.update({(a, b): s for a, b, s in zip(S.gid, S.variant_id, S.slope)}); af.update(dict(zip(S.variant_id, S.af)))
    rows = []
    for g1, g2 in zip(P.g1, P.g2):
        if g1 not in idx or g2 not in idx: continue
        score, pred, lead_pred, cov_w, tot_w = 0.0, 0.0, 0.0, 0.0, 0.0
        for p1 in CS.get(g1, []):
            for p2 in CS.get(g2, []):
                sh = p1.keys() & p2.keys()
                if not sh: continue
                score = max(score, sum(p1[v] * p2[v] for v in sh))
                lv = max(sh, key=lambda v: p1[v] * p2[v]); ws = sum(p1[v] * p2[v] for v in sh)
                for v in sh:
                    w = p1[v] * p2[v]; tot_w += w
                    if (g1, v) in slope and (g2, v) in slope:
                        p = af[v]; pred += w * 2 * p * (1 - p) * slope[(g1, v)] * slope[(g2, v)]; cov_w += w
                if (g1, lv) in slope and (g2, lv) in slope:
                    p = af[lv]; lead_pred += ws * 2 * p * (1 - p) * slope[(g1, lv)] * slope[(g2, lv)]
        if tot_w == 0: continue
        rows.append({'tissue': t, 'g1': g1, 'g2': g2, 'coloc_score': score, 'pred_r': pred, 'pred_r_lead_only': lead_pred,
                     'weight_with_slopes': cov_w / tot_w, 'obs_r_int_pc15': corr(Xp, idx[g1], idx[g2]), 'obs_r_int_raw': corr(Xc, idx[g1], idx[g2])})
    R = pd.DataFrame(rows); upsert_csv(R, OUT, ['tissue', 'g1', 'g2'])
    print(t, len(R), 'pairs | slope coverage %.2f | r(pred, obs_pc15) %.2f' % (R.weight_with_slopes.mean(), R[['pred_r', 'obs_r_int_pc15']].corr().iloc[0, 1]), flush=True)
    del bed, X, Xc, Xp, U_, Vt
