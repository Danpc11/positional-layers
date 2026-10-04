import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
"""Calibrated eQTL law: lead shared fine-mapped variant -> its GTEx slopes (inverse-normal expression, variance ~1) for both genes
-> predicted genetic correlation 2p(1-p)*b1*b2 (times P(same variant)). Compared with the observed coupling."""
import sys, os, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
norm = lambda s: s.lower().replace('-', '_')
for t in sys.argv[1:]:
    out = fDATA + 'predcal_{t}.csv.gz'
    if os.path.exists(out): continue
    P = pd.read_csv(fDATA + 'pred_pairs_{t}.csv.gz'); P = P[P.dist > 0]; genes = set(P.g1) | set(P.g2)
    sf = [f for f in glob.glob(DATA + 'susie/*.parquet') if norm(os.path.basename(f).replace('_v11_eQTLs_SuSiE_summary.parquet', '')) == t][0]
    E = pq.read_table(sf, columns=['phenotype_id', 'variant_id', 'pip', 'cs_id']).to_pandas(); E['gid'] = E.phenotype_id.str.split('.').str[0]; E = E[E.gid.isin(genes)]
    CS = {g: [dict(zip(c.variant_id, c.pip)) for _, c in d.groupby('cs_id')] for g, d in E.groupby('gid')}
    lead = {}
    for g1, g2 in zip(P.g1, P.g2):
        best = (0, None)
        for p1 in CS.get(g1, []):
            for p2 in CS.get(g2, []):
                sh = p1.keys() & p2.keys()
                if sh:
                    ps = sum(p1[v] * p2[v] for v in sh); v = max(sh, key=lambda v: p1[v] * p2[v])
                    if ps > best[0]: best = (ps, v)
        lead[(g1, g2)] = best
    need = {v for _, v in lead.values() if v}
    qf = [f for f in glob.glob(DATA + 'eqtl/*.parquet') if norm(os.path.basename(f).replace('_v11_eQTLs_signif_pairs.parquet', '')) == t][0]
    pf = pq.ParquetFile(qf); slope = {}; af = {}
    for rg in range(pf.num_row_groups):
        S = pf.read_row_group(rg, columns=['phenotype_id', 'variant_id', 'slope', 'af']).to_pandas(); S = S[S.variant_id.isin(need)]
        S['gid'] = S.phenotype_id.str.split('.').str[0]; S = S[S.gid.isin(genes)]
        slope.update({(a, b): s for a, b, s in zip(S.gid, S.variant_id, S.slope)}); af.update(dict(zip(S.variant_id, S.af)))
    rows = []
    for (g1, g2), (ps, v) in lead.items():
        if v is None or (g1, v) not in slope or (g2, v) not in slope: continue
        p = af[v]; rows.append({'g1': g1, 'g2': g2, 'p_same': ps, 'pred_genetic_r': ps * 2 * p * (1 - p) * slope[(g1, v)] * slope[(g2, v)]})
    R = pd.DataFrame(rows).merge(P[['g1', 'g2', 'obs_r']], on=['g1', 'g2']); R['tissue'] = t; R.to_csv(out, index=False)
    print(t, len(R), 'pairs | Pearson %.2f' % R[['pred_genetic_r', 'obs_r']].corr().iloc[0, 1], flush=True)
