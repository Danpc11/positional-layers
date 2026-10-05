"""Do therapeutic CRISPR edits perturb the cis neighbours of the edited locus, in proportion to their baseline coupling?
Baseline coupling: GTEx whole blood (GC- and technically corrected), computed earlier. Edits: BCL11A erythroid enhancer
(chr2) and HBG1/2 promoters (chr11, beta-globin locus). Predictions are made from genomic position alone, before looking
at the edited-cell data."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
C = pd.read_csv(DATA + 'edit/GSE264491_merged_counts.csv.gz', index_col=0)
grp = ['unedited'] * 3 + ['BCL11A_enhancer'] * 3 + ['HBG1_2_promoter'] * 3
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
ann = G.loc[G.index.intersection(C.index), ['Chromosome', 'Start', 'End', 'gene_name']]; ann.columns = ['chr', 'start', 'end', 'sym']; ann['chr'] = ann.chr.astype(str)
C = C.loc[ann.index]; C = C[C.median(axis=1) >= 10]; ann = ann.loc[C.index]
cpm = np.log2(C / C.sum(0) * 1e6 + 1)
def lfc(g): return cpm[[c for c, x in zip(C.columns, grp) if x == g]].mean(1) - cpm[[c for c, x in zip(C.columns, grp) if x == 'unedited']].mean(1)
L = pd.DataFrame({'BCL11A_enhancer': lfc('BCL11A_enhancer'), 'HBG1_2_promoter': lfc('HBG1_2_promoter')})
sd = cpm[[c for c, x in zip(C.columns, grp) if x == 'unedited']].std(1)
# baseline coupling from GTEx whole blood
P = pd.read_csv(OUTDIR + 'pairs_whole_blood.csv.gz'); nb = {}
for a, b, r in zip(P.g1, P.g2, P.r): nb.setdefault(a, []).append((b, r)); nb.setdefault(b, []).append((a, r))
def neighbours(target_sym, k=12):
    tid = ann.index[ann.sym == target_sym]
    if not len(tid): return pd.DataFrame()
    t = tid[0]; c = ann.loc[t, 'chr']; sub = ann[ann.chr == c].sort_values('start'); pos = list(sub.index).index(t)
    out = []
    for j in range(max(0, pos - k), min(len(sub), pos + k + 1)):
        if j == pos: continue
        gid = sub.index[j]; r = dict(nb.get(t, [])).get(gid, np.nan)
        out.append({'gene': sub.loc[gid, 'sym'], 'gid': gid, 'rank_from_target': j - pos, 'dist_kb': (sub.loc[gid, 'start'] - ann.loc[t, 'end']) / 1e3,
                    'baseline_coupling_r': r, 'lfc_BCL11A_edit': L.loc[gid, 'BCL11A_enhancer'], 'lfc_HBG_edit': L.loc[gid, 'HBG1_2_promoter']})
    return pd.DataFrame(out)
print('=== Neighbourhood of BCL11A (chr2), edited target of the approved therapy ===')
NB = neighbours('BCL11A'); print(NB.round(3).to_string(index=False))
print('\n=== Neighbourhood of HBG1 (chr11, beta-globin locus), target of reni-cel ===')
NH = neighbours('HBG1'); print(NH.round(3).to_string(index=False))
# local-window enrichment vs random genome-wide background (not matched on expression or local structure): |lfc| of genes within +/-k positions of the edited gene
def window_test(sym, col, k=10, nperm=20000):
    tid = ann.index[ann.sym == sym][0]; c = ann.loc[tid, 'chr']; sub = ann[ann.chr == c].sort_values('start'); pos = list(sub.index).index(tid)
    idx = [sub.index[j] for j in range(max(0, pos - k), min(len(sub), pos + k + 1)) if j != pos]
    obs = L.loc[idx, col].abs().mean(); rng = np.random.default_rng(0); allg = L.index
    null = np.array([L.loc[rng.choice(allg, len(idx), replace=False), col].abs().mean() for _ in range(nperm)])
    return obs, null.mean(), (np.sum(null >= obs) + 1) / (nperm + 1), len(idx)
rows = []
for sym in ['BCL11A', 'HBG1']:
    for col in ['BCL11A_enhancer', 'HBG1_2_promoter']:
        o, n, p, k = window_test(sym, col); rows.append({'edited_locus_window': sym, 'edit': col, 'n_neighbours': k, 'mean_|lfc|_window': o, 'mean_|lfc|_random': n, 'P': p})
W = pd.DataFrame(rows); print('\n=== Local perturbation around each edited locus ===\n' + W.round(4).to_string(index=False))
# genome-wide: does baseline coupling to a strongly changed gene predict a neighbour's change?
res = []
for col in ['BCL11A_enhancer', 'HBG1_2_promoter']:
    d = P.assign(l1=L[col].reindex(P.g1).values, l2=L[col].reindex(P.g2).values, s1=sd.reindex(P.g1).values, s2=sd.reindex(P.g2).values).dropna()
    d = d[d.dist > 0]; big = d[(d.l1.abs() > 0.5) | (d.l2.abs() > 0.5)].copy()
    drv = np.where(big.l1.abs() >= big.l2.abs(), big.l1, big.l2); nbv = np.where(big.l1.abs() >= big.l2.abs(), big.l2, big.l1)
    res.append({'edit': col, 'n_pairs_with_driver': len(big), 'pearson_pred_vs_obs': stats.pearsonr(big.r * drv, nbv)[0], 'slope': np.polyfit(big.r * drv, nbv, 1)[0],
                'corr_lfc_neighbours_all_pairs': stats.pearsonr(d.l1, d.l2)[0], 'corr_shuffled': stats.pearsonr(d.l1, np.random.default_rng(1).permutation(d.l2.values))[0]})
R = pd.DataFrame(res); print('\n=== Bystander model genome-wide ===\n' + R.round(3).to_string(index=False))
NB.to_csv(OUTDIR + 'edit_BCL11A_neighbours.csv', index=False); NH.to_csv(OUTDIR + 'edit_HBG_neighbours.csv', index=False)
W.to_csv(OUTDIR + 'edit_window_tests.csv', index=False); R.to_csv(OUTDIR + 'edit_bystander_model.csv', index=False)
