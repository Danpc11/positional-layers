"""Gene TSS table (hg19) used by hic_extract.py: one row per Ensembl gene on chromosomes 1-22 and X, strand-aware TSS from Ensembl GRCh37."""
import pandas as pd, pyannotables as pa
import os
from poslayers.config import DATA
os.makedirs(DATA + 'hic', exist_ok=True)
G = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G = G[~G.index.duplicated()]
G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)] + ['X'])]
tss = G.Start.where(G.Strand.astype(str).isin(['1', '+']), G.End)
pd.DataFrame({'gene_id': G.index, 'chr': G.Chromosome.astype(str).values, 'tss': tss.astype(int).values}).to_csv(DATA + 'hic/genes_hg19_tss.tsv', sep='\t', index=False)
