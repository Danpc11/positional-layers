import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
import sys
THR = float(sys.argv[1]) if len(sys.argv) > 1 else 10
import glob, os, re, numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
G['chr'] = G.Chromosome.astype(str); G['len'] = G.End - G.Start
# ---- SLAM-seq: newly synthesised reads per gene (Entrez), mapped to Ensembl by 3' UTR overlap on the same strand
def load(f):
    d = pd.read_csv(f, sep='\t'); return d.groupby('Name').agg(chr=('Chromosome', 'first'), s=('Start', 'min'), e=('End', 'max'), strand=('Strand', 'first'), tc=('TcReadCount', 'sum'), tot=('ReadCount', 'sum'))
files = {re.search(r'Sample(\d+)', f).group(1): f for f in glob.glob('GSM*_tcount.tsv.gz')}
ref = load(files['26']); ref['chr'] = ref.chr.str.replace('chr', '')
emap = {}
for c, d in ref.groupby('chr'):
    g = G[G.chr == c]
    for ez, r in d.iterrows():
        hit = g[(g.Start <= r.e) & (g.End >= r.s) & (g.Strand == r.strand)]
        if len(hit) == 1: emap[ez] = hit.index[0]
TC = pd.DataFrame({k: load(f).tc for k, f in files.items()}).fillna(0); TC = TC[TC.index.isin(emap)]; TC.index = [emap[i] for i in TC.index]; TC = TC.groupby(level=0).sum()
print('genes mapped:', len(TC))
# ---- baseline coupling in myeloid leukaemia (BeatAML2), adjacent genes in GRCh38 order, GC-corrected
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
B = pd.read_csv(DATA + 'beataml_waves1to4_counts_dbgap.txt', sep='\t'); B['gid'] = B.stable_id.str.split('.').str[0]
B = B.drop_duplicates('gid').set_index('gid')[[x for x in B.columns if x.startswith('BA')]]; B = B.loc[B.index.intersection(G.index).intersection(BM.index)]; B = B[B.median(axis=1) >= 10]
g = G.loc[B.index].sort_values(['chr', 'Start']); B = B.loc[g.index]
Y = np.log2(B.values / B.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[g.index, 'gc'].values; gcz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]; Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
same = g.chr.values[:-1] == g.chr.values[1:]; i = np.where(same)[0]
P = pd.DataFrame({'g1': g.index[i], 'g2': g.index[i + 1], 'r': (Z[i] * Z[i + 1]).mean(1), 'dist': g.Start.values[i + 1] - g.End.values[i]}); P = P[P.dist > 0]
print('baseline pairs:', len(P), '| median r %.3f' % P.r.median())
# ---- drug responses
COMP = [('MOLM-13', 'JQ1 (exp 1)', 'chromatin/transcription', ['26', '27', '28'], ['29', '30', '31']), ('MOLM-13', 'JQ1 (exp 2)', 'chromatin/transcription', ['35', '36', '37'], ['38', '39', '40']),
        ('MV4-11', 'JQ1', 'chromatin/transcription', ['20', '21', '22'], ['23', '24', '25']), ('K562', 'JQ1', 'chromatin/transcription', ['14', '15', '16'], ['17', '18', '19']),
        ('K562', 'BRD4 degradation', 'chromatin/transcription', ['8', '9', '10'], ['11', '12', '13']), ('MOLM-13', 'NVP-2 6 nM (CDK9i)', 'chromatin/transcription', ['35', '36', '37'], ['41', '42', '43']),
        ('MOLM-13', 'NVP-2 60 nM (CDK9i)', 'chromatin/transcription', ['35', '36', '37'], ['44', '45']), ('K562', 'flavopiridol (CDK9i)', 'chromatin/transcription', ['1'], ['2']),
        ('K562', 'nilotinib (BCR-ABL)', 'signalling', ['3'], ['7']), ('K562', 'trametinib (MEK)', 'signalling', ['3'], ['5']), ('K562', 'MK-2206 (AKT)', 'signalling', ['3'], ['4']), ('K562', 'MK-2206 + trametinib', 'signalling', ['3'], ['6'])]
lc = np.log10(G.reindex(TC.index)['len'].values.astype(float)); rows = []
for cell, drug, cls, ctrl, trt in COMP:
    X = np.log2(TC[ctrl + trt] / TC[ctrl + trt].sum() * 1e6 + 1); keep = TC[ctrl].mean(1) >= THR
    a, b = X[trt].values, X[ctrl].values; lfc = a.mean(1) - b.mean(1)
    stat = lfc / np.sqrt(a.var(1, ddof=1) / a.shape[1] + b.var(1, ddof=1) / b.shape[1] + 0.01) if len(trt) > 1 and len(ctrl) > 1 else lfc
    s = pd.Series(stat, index=TC.index)[keep]; base = b.mean(1)[keep.values]
    Xc = np.column_stack([np.ones(keep.sum()), lc[keep.values], base, base ** 2]); ok = np.isfinite(Xc).all(1)
    res = pd.Series(np.nan, index=s.index); res[ok] = s[ok].values - Xc[ok] @ np.linalg.lstsq(Xc[ok], s[ok].values, rcond=None)[0]   # remove length/expression-level trends
    d = P.assign(t1=res.reindex(P.g1).values, t2=res.reindex(P.g2).values).dropna()
    z1 = (d.t1 - d.t1.mean()) / d.t1.std(); z2 = (d.t2 - d.t2.mean()) / d.t2.std(); prod = z1 * z2
    slope, _, _, p, _ = stats.linregress(d.r, prod); q = pd.qcut(d.r, 5, labels=False)
    rows.append({'cell': cell, 'perturbation': drug, 'class': cls, 'replicates': f'{len(ctrl)}v{len(trt)}', 'pairs': len(d), 'slope_on_baseline_coupling': slope, 'p': p,
                 'concordance_low_coupling': np.corrcoef(d.t1[q == 0], d.t2[q == 0])[0, 1], 'concordance_high_coupling': np.corrcoef(d.t1[q == 4], d.t2[q == 4])[0, 1]})
R = pd.DataFrame(rows); R.to_csv(f'drug_cis_results_thr{int(THR)}.csv', index=False); pd.set_option('display.width', 220); print(R.round(3).to_string(index=False))
print('\nmedian slope by class:'); print(R.groupby('class').slope_on_baseline_coupling.median().round(3).to_string())
print('Mann-Whitney chromatin vs signalling P = %.3f' % stats.mannwhitneyu(R[R['class'] == 'chromatin/transcription'].slope_on_baseline_coupling, R[R['class'] == 'signalling'].slope_on_baseline_coupling).pvalue)
