"""Autocorrelation of gene GC content along the genome (Extended Data Fig. 2d).

Protein-coding genes on chromosomes 1-22 (Ensembl 100, GRCh38) in genomic order; GC content standardised within each
chromosome; lag-L autocorrelation averaged over chromosomes with more than L + 5 genes, for L = 1..60.

Input : DATA/biomart_GRCh38_gene_gc.txt (with the 'Gene type' column)
Output: OUTDIR/gc_autocorrelation.csv (lag, autocorrelation, genes)
"""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA, OUTDIR

BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(
    columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])].join(BM[['gc', 'type']], how='inner')
G = G[G.type == 'protein_coding'].sort_values(['Chromosome', 'Start'])
rows = []
for L in range(1, 61):
    v = []
    for _, g in G.groupby('Chromosome'):
        x = (g.gc.values - g.gc.mean()) / g.gc.std()
        if len(x) > L + 5: v.append(np.mean(x[:-L] * x[L:]))
    rows.append({'lag': L, 'autocorrelation': float(np.mean(v)), 'genes': len(G)})
pd.DataFrame(rows).to_csv(OUTDIR + 'gc_autocorrelation.csv', index=False)
print(f'{len(G):,} genes; lag 1 {rows[0]["autocorrelation"]:.2f}, lag 10 {rows[9]["autocorrelation"]:.2f}, lag 30 {rows[29]["autocorrelation"]:.2f}')
