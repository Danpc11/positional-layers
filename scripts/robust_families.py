"""Gene families, classic clusters and read-through transcripts (Fig. 2d). Input: pairs_<tissue>.csv.gz from atlas.py."""
import glob, os, re, numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import OUTDIR
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; sym = G.gene_name.astype(str)
stem = lambda s: (re.match(r'^([A-Z]+)\d', s).group(1) if re.match(r'^([A-Z]+)\d', s) else None)
CLUST = ('HIST', 'H2A', 'H2B', 'H3C', 'H4C', 'H1-', 'OR', 'KRT', 'KRTAP', 'PCDH', 'ZNF', 'UGT', 'CYP', 'HLA', 'IGH', 'IGK', 'IGL', 'TRB', 'TRA', 'KIR', 'LCE', 'SPRR', 'MT1', 'MT2',
         'DEF', 'SERPIN', 'APO', 'GST', 'CLEC', 'LILR', 'SIGLEC', 'TAS2R', 'PRAME', 'GAGE', 'MAGE', 'SPANX', 'USP17', 'NBPF', 'TBC1D3', 'GOLGA6', 'FAM90')
def readthrough(s): return ('-' in s) and not re.search(r'-(AS|DT|IT|OT|OS)\d*$', s) and not s.startswith(('HLA-', 'H1-', 'H2A', 'H2B', 'H3-', 'H4-'))
rows = []
for f in sorted(glob.glob(OUTDIR + 'pairs_*.csv.gz')):
    t = os.path.basename(f)[6:-7]; P = pd.read_csv(f); P = P[P.dist > 0].copy()
    s1 = sym.reindex(P.g1).fillna('').values; s2 = sym.reindex(P.g2).fillna('').values
    st1 = np.array([stem(x) for x in s1], dtype=object); st2 = np.array([stem(x) for x in s2], dtype=object)
    same_stem = (st1 == st2) & pd.notna(st1)
    clus = np.array([a.startswith(CLUST) or b.startswith(CLUST) for a, b in zip(s1, s2)])
    rt = np.array([readthrough(a) or readthrough(b) for a, b in zip(s1, s2)])
    keep = ~(same_stem | clus | rt); base = P.r.mean()
    rows.append({'tissue': t, 'pairs': len(P), 'excluded_pct': 100 * (1 - keep.mean()), 'r_all': base, 'r_same_stem': P.r[same_stem].mean(), 'r_clusters': P.r[clus].mean(),
                 'r_readthrough': P.r[rt].mean(), 'r_clean': P.r[keep].mean(), 'clean_vs_all_pct': 100 * (P.r[keep].mean() / base - 1),
                 'r_clean_le20kb': P.r[keep & (P.dist <= 2e4)].mean(), 'r_clean_gt100kb': P.r[keep & (P.dist > 1e5)].mean()})
pd.DataFrame(rows).to_csv(OUTDIR + 'robust_families.csv', index=False); print(len(rows), 'tissues')
