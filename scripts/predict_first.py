"""STEP 1 — predictions fixed from GTEx whole blood only (baseline cis coupling), blind to the edited-cell data."""
import pandas as pd, numpy as np, pyannotables as pa
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
P = pd.read_csv('/home/claude/atlas/pairs_whole_blood.csv.gz')
P['s1'] = sym.reindex(P.g1).values; P['s2'] = sym.reindex(P.g2).values
def neighbours(target, k=6):
    """genes within k positions of the target on the GTEx-blood gene order, with the coupling of the intervening pairs"""
    i1 = P.index[(P.s1 == target) | (P.s2 == target)]
    if not len(i1): return pd.DataFrame()
    rows = []
    for i in i1:
        r = P.loc[i]; nb = r.s2 if r.s1 == target else r.s1
        rows.append({'target': target, 'neighbour': nb, 'steps': 1, 'coupling_r': r.r, 'dist_kb': r.dist / 1e3})
    # two steps away: product of consecutive couplings (expected transmitted coupling)
    for i in i1:
        r = P.loc[i]; mid = r.s2 if r.s1 == target else r.s1
        j = P.index[((P.s1 == mid) | (P.s2 == mid))]
        for jj in j:
            r2 = P.loc[jj]; nb2 = r2.s2 if r2.s1 == mid else r2.s1
            if nb2 in (target, mid): continue
            rows.append({'target': target, 'neighbour': nb2, 'steps': 2, 'coupling_r': r.r * r2.r, 'dist_kb': np.nan})
    return pd.DataFrame(rows)
pred = pd.concat([neighbours(t) for t in ['BCL11A', 'HBG1', 'HBG2', 'HBB', 'HBD']]).drop_duplicates(['target', 'neighbour'])
pred['predicted_direction'] = np.where(pred.coupling_r > 0, 'same as target', 'opposite to target')
pred = pred.sort_values(['target', 'coupling_r'], ascending=[True, False])
pred.to_csv('predictions_from_GTEx_blood.csv', index=False)
pd.set_option('display.width', 200); print(pred.round(3).to_string(index=False))
print('\nblood-wide reference: median |r| %.3f, 90th percentile %.3f' % (P.r.abs().median(), P.r.quantile(0.9)))
