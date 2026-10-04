"""STEP 2 — test the fixed predictions in the edited erythroblasts (GSE264491)."""
import pandas as pd, numpy as np, pyannotables as pa
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
C = pd.read_csv('GSE264491_merged_counts_csv.gz', index_col=0)
grp = {'unedited': C.columns[:3], 'BCL11A': C.columns[3:6], 'HBG': C.columns[6:9]}
C = C[C.median(axis=1) >= 10]; cpm = np.log2(C / C.sum(0) * 1e6 + 1)
def lfc(g): return cpm[grp[g]].mean(1) - cpm[grp['unedited']].mean(1)
def tstat(g):
    a = cpm[grp[g]].values; b = cpm[grp['unedited']].values
    return (a.mean(1) - b.mean(1)) / np.sqrt(a.var(1, ddof=1) / 3 + b.var(1, ddof=1) / 3 + 1e-9)
R = pd.DataFrame({'sym': sym.reindex(C.index).values, 'lfc_BCL11A': lfc('BCL11A'), 'lfc_HBG': lfc('HBG'), 't_BCL11A': tstat('BCL11A'), 't_HBG': tstat('HBG')}, index=C.index)
print('== on-target effects'); print(R[R.sym.isin(['BCL11A', 'HBG1', 'HBG2', 'HBB', 'HBD', 'HBE1', 'HBBP1'])].round(2).to_string(index=False))
pred = pd.read_csv('predictions_from_GTEx_blood.csv'); s2i = {v: k for k, v in sym.items() if v in set(pred.neighbour)}
pred['gid'] = pred.neighbour.map(s2i); pred = pred[pred.gid.isin(R.index)]
pred['lfc_BCL11A'] = R.lfc_BCL11A.reindex(pred.gid).values; pred['lfc_HBG'] = R.lfc_HBG.reindex(pred.gid).values
print('\n== predicted neighbours, observed log2 fold change')
print(pred[['target', 'neighbour', 'steps', 'coupling_r', 'lfc_BCL11A', 'lfc_HBG']].round(3).to_string(index=False))
# genome-wide: do neighbours of the edited locus move with their baseline coupling?
P = pd.read_csv('/home/claude/atlas/pairs_whole_blood.csv.gz'); P = P[P.g1.isin(R.index) & P.g2.isin(R.index)]
for edit, tcol in [('BCL11A', 't_BCL11A'), ('HBG', 't_HBG')]:
    d = P.assign(t1=R[tcol].reindex(P.g1).values, t2=R[tcol].reindex(P.g2).values).dropna(subset=['t1', 't2'])
    d['bin'] = pd.qcut(d.r, 5)
    g = d.groupby('bin', observed=True).apply(lambda x: pd.Series({'n': len(x), 'baseline_r': x.r.mean(), 'corr_of_edit_effects': stats.pearsonr(x.t1, x.t2)[0]}), include_groups=False)
    print(f'\n== {edit} edit: concordance of neighbour responses by baseline coupling'); print(g.round(3).to_string())
    print('  slope of concordance on baseline coupling: %.2f' % np.polyfit(d.r, (d.t1 - d.t1.mean()) / d.t1.std() * (d.t2 - d.t2.mean()) / d.t2.std(), 1)[0])
R.to_csv('edit_effects.csv'); pred.to_csv('predictions_tested.csv', index=False)
# cis vs trans: distance of responding genes from the edited locus
loc = {'BCL11A': ('2', 60450000), 'HBG': ('11', 5250000)}
pos = pd.DataFrame({'chr': G.Chromosome.astype(str).reindex(R.index), 'mid': ((G.Start + G.End) / 2).reindex(R.index)})
for edit, tcol in [('BCL11A', 't_BCL11A'), ('HBG', 't_HBG')]:
    c, p0 = loc[edit]; near = (pos.chr == c) & ((pos.mid - p0).abs() < 2e6)
    big = R[tcol].abs() > 4
    print(f'\n{edit} edit: responding genes (|t|>4) = {int(big.sum())}; of these within 2 Mb of the edited site = {int((big & near).sum())}; genes within 2 Mb = {int(near.sum())}')
    if (big & near).sum(): print('  ', R.sym[big & near].tolist())
