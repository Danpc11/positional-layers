"""Continuous gene-level copy number from GDC segment log2 ratios, for the cohorts already extracted."""
from poslayers.config import DATA, OUTDIR
import numpy as np, pandas as pd, pyannotables as pa
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
mid = ((G.Start + G.End) / 2).values; gch = G.Chromosome.astype(str).values; gid = G.index.values
for c in ['BLCA', 'UCEC', 'STAD', 'COAD', 'GBM']:
    z = np.load(OUTDIR + f'expr_{c}.npz', allow_pickle=True); samp = set(z['samples']); genes = set(z['genes'])
    sel = np.isin(gid, list(genes)); gi, gm, gc = gid[sel], mid[sel], gch[sel]
    parts = [ch[ch['sample'].isin(samp)] for ch in pd.read_csv(DATA + 'tcga/GDC-PANCAN_cnv.tsv', sep='\t', chunksize=1_000_000, dtype={'Chrom': str})]
    S = pd.concat(parts); S['Chrom'] = S.Chrom.str.replace('chr', ''); smp = sorted(S['sample'].unique())
    M = np.zeros((len(gi), len(smp)), dtype=np.float32)
    for k, (s, d) in enumerate(S.groupby('sample')):
        for ch, dd in d.groupby('Chrom'):
            idx = np.where(gc == ch)[0]
            if not len(idx): continue
            dd = dd.sort_values('Start'); j = np.searchsorted(dd.Start.values, gm[idx], side='right') - 1
            ok = (j >= 0) & (gm[idx] <= dd.End.values[np.clip(j, 0, None)]); M[idx[ok], smp.index(s)] = dd.value.values[j[ok]]
    np.savez_compressed(OUTDIR + f'cn_cont_{c}.npz', CN=M, genes=gi, samples=np.array(smp)); print(c, M.shape, 'segments', len(S), flush=True)
