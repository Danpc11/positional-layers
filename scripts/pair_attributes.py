"""Attributes of adjacent gene pairs and their coupling, for the UpSet plot in Fig. 2d.
 
For each of the 14 tissues with a matched TAD map, every adjacent pair with both genes inside a TAD is flagged for three
attributes: gene family, classic cluster or read-through transcript (as in robust_families.py); a shared significant eQTL
(shares_<tissue>.csv.gz from eqtl_share.py); and membership of the same TAD (lib_tad.py). Coupling is adjusted for distance
by subtracting the mean coupling of the pair's distance bin and adding back the tissue mean, because same-TAD pairs are closer.
 
Inputs : OUTDIR/pairs_<tissue>.csv.gz (orient.py), OUTDIR/shares_<tissue>.csv.gz (eqtl_share.py), DATA/TAD-full/
Output : OUTDIR/upset_pair_attributes.csv, one row per tissue and attribute combination
"""
import glob, os, re, numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR
from lib_tad import MATCH, tads, assign
 
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
stem = lambda s: (re.match(r'^([A-Z]+)\d', s).group(1) if re.match(r'^([A-Z]+)\d', s) else None)
CLUST = ('HIST', 'H2A', 'H2B', 'H3C', 'H4C', 'H1-', 'OR', 'KRT', 'KRTAP', 'PCDH', 'ZNF', 'UGT', 'CYP', 'HLA', 'IGH', 'IGK', 'IGL', 'TRB', 'TRA', 'KIR', 'LCE', 'SPRR', 'MT1', 'MT2',
         'DEF', 'SERPIN', 'APO', 'GST', 'CLEC', 'LILR', 'SIGLEC', 'TAS2R', 'PRAME', 'GAGE', 'MAGE', 'SPANX', 'USP17', 'NBPF', 'TBC1D3', 'GOLGA6', 'FAM90')
def readthrough(s): return ('-' in s) and not re.search(r'-(AS|DT|IT|OT|OS)\d*$', s) and not s.startswith(('HLA-', 'H1-', 'H2A', 'H2B', 'H3-', 'H4-'))
BINS = [0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
norm = lambda s: s.lower().replace('-', '_')
PAIRS = {norm(os.path.basename(f)[6:-7]): f for f in glob.glob(OUTDIR + 'pairs_*.csv.gz')}
SHARES = {norm(os.path.basename(f)[7:-7]): f for f in glob.glob(OUTDIR + 'shares_*.csv.gz')}
 
rows = []
for t, tad_map in MATCH.items():
    if t not in PAIRS or t not in SHARES:
        print('missing inputs for', t); continue
    d = pd.read_csv(PAIRS[t]).merge(pd.read_csv(SHARES[t])[['g1', 'g2', 'share_any']], on=['g1', 'g2']); d = d[d.dist > 0].copy()
    a = assign(tads(tad_map)); d['t1'] = a.reindex(d.g1).values; d['t2'] = a.reindex(d.g2).values; d = d[(d.t1 >= 0) & (d.t2 >= 0)]
    s1 = sym.reindex(d.g1).fillna('').values; s2 = sym.reindex(d.g2).fillna('').values
    d['family'] = [(stem(x) is not None and stem(x) == stem(y)) or x.startswith(CLUST) or y.startswith(CLUST) or readthrough(x) or readthrough(y) for x, y in zip(s1, s2)]
    d['shared_eqtl'] = d.share_any.astype(bool); d['same_tad'] = d.t1 == d.t2
    d['bin'] = pd.cut(d.dist, BINS); d['r_adj'] = d.r - d.groupby('bin', observed=True).r.transform('mean') + d.r.mean()
    for key, g in d.groupby(['family', 'shared_eqtl', 'same_tad']):
        rows.append({'tissue': t, 'family': key[0], 'shared_eqtl': key[1], 'same_tad': key[2], 'pairs': len(g), 'frac': len(g) / len(d), 'r_adj': g.r_adj.mean(), 'r': g.r.mean()})
U = pd.DataFrame(rows); U.to_csv(OUTDIR + 'upset_pair_attributes.csv', index=False)
print(U.groupby(['family', 'shared_eqtl', 'same_tad']).agg(pairs=('pairs', 'median'), r_adj=('r_adj', 'median'), tissues=('tissue', 'size')).round(3).to_string())
 
