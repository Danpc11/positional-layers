"""Gene-editing neighbour predictions, evaluated with an explicit rule (Fig. 5c,d; Supplementary Table 7B).

Groups are read from the GEO series matrix (genotype of each GSM, linked to the count columns through the supplementary
file names), not from column order. For each edit versus unedited cells (n = 3 v 3): log2 fold change, Welch t, Welch
degrees of freedom, two-sided P and Benjamini-Hochberg q over all expressed genes.

predictions_from_GTEx_blood.csv (20 neighbour predictions; target BCL11A -> BCL11A-enhancer edit, target HBB -> HBG1/2
promoter edit). Its SHA-256 and modification time are recorded; the file was written before the edit effects were
computed but after the count matrix had been downloaded. A prediction is
  'direction met'  when sign(neighbour lfc) x sign(target lfc) matches the predicted relation (same / opposite);
  'supported'      when the direction is met AND the neighbour responds at q < 0.10.
Outputs: OUTDIR/edit_stats_<edit>.csv, OUTDIR/edit_predictions_evaluated.csv
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # scripts/ for lib_*.py
import gzip, hashlib, os, re, datetime, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from poslayers.config import DATA, OUTDIR

lines = [l.rstrip('\n').split('\t') for l in gzip.open(DATA + 'edit/GSE264491_series_matrix_txt.gz', 'rt') if l.startswith('!Sample_')]
meta = {}
for l in lines:
    key = l[0]; vals = [v.strip('"') for v in l[1:]]
    if key == '!Sample_supplementary_file_1': meta['sna'] = [re.search(r'(SNA\d+)', v).group(1) for v in vals]
    if key == '!Sample_characteristics_ch1' and vals[0].startswith('genotype'): meta['genotype'] = [v.split(': ', 1)[1] for v in vals]
M = pd.DataFrame(meta).set_index('sna')
C = pd.read_csv(DATA + 'edit/GSE264491_merged_counts.csv.gz', index_col=0); C.index = C.index.str.split('.').str[0]
grp = {'unedited': M.index[M.genotype == 'Wild Type'], 'BCL11A': M.index[M.genotype.str.contains('BCL11A')], 'HBG': M.index[M.genotype.str.contains('HBG')]}
assert all(len(v) == 3 for v in grp.values()) and all(set(v) <= set(C.columns) for v in grp.values()), grp
cpm = np.log2(C / C.sum(0) * 1e6 + 1); keep = (C[grp['unedited']].mean(1) >= 10)
import pyannotables as pa
sym = pa.tables()['homo_sapiens-GRCh38-ensembl100']; sym = sym[~sym.index.duplicated()].gene_name.astype(str)
st = {}
for e in ('BCL11A', 'HBG'):
    a, b = cpm.loc[keep, grp[e]].values, cpm.loc[keep, grp['unedited']].values
    va, vb = a.var(1, ddof=1) / 3, b.var(1, ddof=1) / 3; t = (a.mean(1) - b.mean(1)) / np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / 2 + vb ** 2 / 2); p = 2 * stats.t.sf(np.abs(t), df)
    ok = np.isfinite(p); q = np.full(len(p), np.nan); q[ok] = multipletests(p[ok], method='fdr_bh')[1]
    S = pd.DataFrame({'gene_id': cpm.index[keep], 'symbol': sym.reindex(cpm.index[keep]).values, 'lfc': a.mean(1) - b.mean(1), 'welch_t': t, 'welch_df': df, 'p': p, 'q_BH': q})
    S.to_csv(OUTDIR + f'edit_stats_{e}.csv', index=False); st[e] = S.set_index('symbol')
pf = DATA + 'edit/predictions_from_GTEx_blood.csv'
sha = hashlib.sha256(open(pf, 'rb').read()).hexdigest(); mtime = datetime.datetime.fromtimestamp(os.path.getmtime(pf)).isoformat(timespec='seconds')
P = pd.read_csv(pf); rows = []
for _, r in P.iterrows():
    e = 'BCL11A' if r.target == 'BCL11A' else 'HBG'; S = st[e]
    tgt = S.lfc.get(r.target, np.nan); nb = S.loc[r.neighbour] if r.neighbour in S.index else None
    if nb is None or isinstance(nb, pd.DataFrame): rows.append({**r.to_dict(), 'edit': e, 'tested': False}); continue
    same = np.sign(nb.lfc) == np.sign(tgt); want_same = r.predicted_direction.startswith('same')
    rows.append({**r.to_dict(), 'edit': e, 'tested': True, 'target_lfc': tgt, 'neighbour_lfc': nb.lfc, 'welch_df': nb.welch_df, 'p': nb.p, 'q_BH': nb.q_BH,
                 'direction_met': bool(same == want_same), 'supported': bool(same == want_same and nb.q_BH < 0.10)})
R = pd.DataFrame(rows); R['predictions_sha256'] = sha; R['predictions_mtime'] = mtime
R.to_csv(OUTDIR + 'edit_predictions_evaluated.csv', index=False)
T = R[R.tested]
print(f'{len(P)} predictions, {len(T)} testable; direction met {int(T.direction_met.sum())}/{len(T)}; supported (direction and q < 0.10) {int(T.supported.sum())}/{len(T)}')
print(T[['target', 'neighbour', 'predicted_direction', 'neighbour_lfc', 'q_BH', 'direction_met', 'supported']].round(3).to_string(index=False))
