"""Reduce the Schmitt et al. 2016 Fit-Hi-C files (GSE87112, 40-kb bins, hg19) to one gene-pair table per tissue.

Usage (on your computer, in the folder that contains FitHiC_primary_cohort/):
    pip install pyannotables pandas numpy
    python3 reduce_fithic.py FitHiC_primary_cohort schmitt_gene_pairs

Output: schmitt_gene_pairs/<CODE>_gene_pairs.csv.gz (about 10 MB each), one row per pair of protein-coding genes on the
same chromosome whose TSSs are 1 to 50 bins apart (40 kb to 2 Mb):
    g1, g2, chr, bin_distance, observed, expected, oe
Bin pairs absent from the sparse file are given observed = 0 and the median expected value at that distance.
The column layout is detected from the data: either RowID, ColumnID (1-based 40-kb bin numbers, as in GSE87112) or
base-pair positions, followed by observed, expected, O/E, p and q; header optional.
"""
import sys, os, re, glob, gzip, numpy as np, pandas as pd, pyannotables as pa

src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True); RES, MAXB = 40000, 50
CODES = sys.argv[3].split(',') if len(sys.argv) > 3 else ['AD', 'AO', 'BL', 'CO', 'HC', 'LG', 'LI', 'LV', 'OV', 'PA', 'PO', 'RV', 'SB', 'SX', 'GM12878', 'imr90']
G = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G = G[~G.index.duplicated()]
G = G[(G.gene_biotype.astype(str) == 'protein_coding') & G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)] + ['X'])]
TSS = pd.DataFrame({'gid': G.index, 'chr': G.Chromosome.astype(str).values, 'tss': np.where(G.Strand.astype(str).isin(['1', '+']), G.Start, G.End).astype(int)})
TSS['bin'] = TSS.tss // RES

def read_fithic(path):
    with gzip.open(path, 'rt') as fh: first = fh.readline().split()
    header = not all(re.fullmatch(r'-?[0-9.eE+]+', x) for x in first[:2])
    d = pd.read_csv(path, sep=r'\s+', header=0 if header else None, engine='python' if header else 'c')
    num = d.select_dtypes('number')
    if num.shape[1] < 5: raise SystemExit(f'unexpected format in {path}: {first}')
    num = num.iloc[:, :5]; num.columns = ['p1', 'p2', 'observed', 'expected', 'oe']
    cols = [str(c).lower() for c in d.columns]
    as_index = ('rowid' in cols and 'columnid' in cols) or num.p2.max() < 1e5     # bin numbers (1-based), not base pairs
    if as_index: b1, b2 = num.p1.astype(int) - 1, num.p2.astype(int) - 1          # bin i covers [(i-1)*40 kb, i*40 kb)
    else: b1, b2 = (num.p1 // RES).astype(int), (num.p2 // RES).astype(int)
    return pd.DataFrame({'lo': np.minimum(b1, b2), 'k': np.abs(b2 - b1), 'observed': num.observed, 'expected': num.expected})

for code in CODES:
    files = sorted(glob.glob(os.path.join(src, f'FitHiC_output.{code}_chr*.sparse.matrix.gz')))
    if not files: print('no files for', code); continue
    parts = []
    for f in files:
        ch = re.search(r'_chr(\w+)\.sparse', f).group(1); ch = 'X' if ch in ('23', 'X') else ch
        t = TSS[TSS.chr == ch].sort_values('tss')
        if len(t) < 2: continue
        M = read_fithic(f); M = M[(M.k >= 1) & (M.k <= MAXB)]
        exp_k = M.groupby('k').expected.median()
        key = M.lo.values.astype(np.int64) * 100 + M.k.values; obs = pd.Series(M.observed.values, index=key); ex = pd.Series(M.expected.values, index=key)
        obs = obs[~obs.index.duplicated()]; ex = ex[~ex.index.duplicated()]
        g, b = t.gid.values, t.bin.values; I, J = [], []
        for i in range(len(t)):
            j = np.arange(i + 1, min(len(t), i + 400)); kk = b[j] - b[i]; sel = j[(kk >= 1) & (kk <= MAXB)]
            I.append(np.full(len(sel), i)); J.append(sel)
        I, J = np.concatenate(I), np.concatenate(J); lo, kk = np.minimum(b[I], b[J]), np.abs(b[J] - b[I]); kq = lo.astype(np.int64) * 100 + kk
        o = obs.reindex(kq).values; e = ex.reindex(kq).values; e = np.where(np.isnan(e), exp_k.reindex(kk).values, e); o = np.nan_to_num(o, nan=0.0)
        parts.append(pd.DataFrame({'g1': g[I], 'g2': g[J], 'chr': ch, 'bin_distance': kk, 'observed': o, 'expected': e, 'oe': o / e}))
    R = pd.concat(parts); R.to_csv(os.path.join(out, f'{code}_gene_pairs.csv.gz'), index=False)
    print(f'{code}: {len(files)} chromosome files, {len(R):,} gene pairs, observed > 0 in {np.mean(R.observed > 0):.1%}', flush=True)
