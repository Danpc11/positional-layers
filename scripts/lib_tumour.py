"""Shared tumour preprocessing for the TCGA analyses (cohesin_v2, cont_tests, aneuploidy, replicate, tad_tumour, cohesin_freedman_lane).
Previously each script read cohesin_v2.py as text, edited it and ran it with exec; the variants are now options of prepare().

prepare(c, cn='gistic' | 'continuous', codes=('01',), return_raw=False, return_gids=False)
    -> (D, [D_raw,] chrs, samp, cov[, gids])
    D       GC-corrected deviations with each gene's own copy-number effect removed
    D_raw   the same before the copy-number correction
    cn      'gistic': remove the mean of each GISTIC level per gene; 'continuous': regress each gene on its own segment log2 ratio (linear + quadratic)
"""
import numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA, OUTDIR
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene name': 'sym'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc', 'sym']], how='inner')
mu = pd.read_csv(DATA + 'tcga/GDC-PANCAN.mutect2_snv.tsv', sep='\t', usecols=['Sample_ID', 'gene', 'effect', 'filter']); mu = mu[mu['filter'] == 'PASS']
NONCOD = ['synonymous_variant', 'intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant', 'non_coding_transcript_exon_variant']
coding = mu[~mu.effect.isin(NONCOD)]; trunc = coding[coding.effect.str.contains('stop_gained|frameshift|splice_acceptor|splice_donor', regex=True)]
burden = mu.groupby('Sample_ID').size()
IMM = ['PTPRC', 'CD2', 'CD3E', 'CD3D', 'CD247', 'LCK', 'CD48', 'CD53', 'CD52', 'CORO1A', 'LAPTM5', 'CD37', 'FCER1G', 'TYROBP', 'CD74', 'HLA-DRA', 'CCL5', 'IL2RG', 'CXCR4', 'SELL']
STR = ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A3', 'DCN', 'LUM', 'FBN1', 'FAP', 'PDGFRB', 'THY1', 'SPARC', 'VCAN', 'MMP2', 'CDH11', 'FN1', 'POSTN', 'COL11A1', 'SFRP2', 'ACTA2']
SUB = {'LAML': (['MPO', 'ELANE', 'AZU1'], ['CD14', 'LYZ', 'CSF1R']), 'GBM': (['OLIG2', 'SOX2', 'PDGFRA'], ['CHI3L1', 'CD44', 'MET']), 'BLCA': (['GATA3', 'KRT20', 'UPK1B', 'UPK2', 'FOXA1', 'PPARG', 'FGFR3'], ['KRT5', 'KRT6A', 'KRT14', 'CD44', 'CDH3']),
       'UCEC': (['PGR', 'ESR1', 'PAX8', 'MSX1'], ['CDKN2A', 'L1CAM', 'WT1']), 'STAD': (['CDX2', 'MUC2', 'TFF3'], ['VIM', 'ZEB1']), 'COAD': (['CDX2', 'VIL1', 'LGALS4'], ['VIM', 'ZEB1'])}
def prepare(c, cn='gistic', codes=('01',), return_raw=False, return_gids=False):
    z = np.load(OUTDIR + f'expr_{c}.npz', allow_pickle=True); X = z['X']; genes = z['genes']; samp = z['samples']
    keep = np.array([s[13:15] in codes for s in samp]) & pd.Index(samp).isin(burden.index); X = X[:, keep]; samp = samp[keep]
    C = np.clip(2 ** X - 1, 0, None); ok = pd.Index(genes).isin(G38.index); C = C[ok]; genes = genes[ok]
    good = np.median(C, 1) >= 10; C = C[good]; genes = genes[good]
    g = G38.loc[genes].reset_index(); o = np.lexsort((g.start.values, g.chr.values)); C = C[o]; g = g.iloc[o].reset_index(drop=True)
    Y = np.log2(C / C.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
    if cn == 'continuous':
        Zc = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CNd = pd.DataFrame(Zc['CN'], index=Zc['genes'], columns=Zc['samples'])
        have = pd.Index(samp).isin(CNd.columns); D = D[:, have]; Y = Y[:, have]; samp = samp[have]
        CNv = CNd.reindex(index=g.iloc[:, 0].values, columns=samp).fillna(0).values; D_raw = D.copy()
        for i in range(D.shape[0]):                          # remove each gene's own continuous dosage effect (linear + quadratic)
            x = CNv[i]
            if np.std(x) < 1e-6: continue
            X1 = np.column_stack([np.ones(len(x)), x, x ** 2]); D[i] -= X1 @ np.linalg.lstsq(X1, D[i], rcond=None)[0]
        cna_burden = np.mean(np.abs(CNv) > 0.2, 0)
    else:
        hdr = pd.read_csv(DATA + 'tcga/GDC-PANCAN.gistic.tsv', sep='\t', nrows=0).columns
        CN = pd.read_csv(DATA + 'tcga/GDC-PANCAN.gistic.tsv', sep='\t', usecols=[hdr[0]] + [x for x in samp if x in set(hdr)], index_col=0); CN.index = CN.index.str.split('.').str[0]
        CN = CN.reindex(index=g.iloc[:, 0].values, columns=samp); have = (CN.notna().mean(axis=0) > 0.5).values
        D = D[:, have]; Y = Y[:, have]; samp = samp[have]; CNv = np.nan_to_num(CN.values[:, have].astype(float), nan=0.0); D_raw = D.copy()
        for i in range(D.shape[0]):
            cnv = CNv[i]
            for lv in np.unique(cnv): D[i, cnv == lv] -= D[i, cnv == lv].mean()
        cna_burden = np.mean(CNv != 0, 0)
    sym = g.sym.astype(str).values; sc = lambda gs: np.nanmean([(Y[sym == s][0] - Y[sym == s][0].mean()) / (Y[sym == s][0].std() + 1e-9) for s in gs if (sym == s).any()], 0)
    up, dn = SUB[c]
    cov = pd.DataFrame({'immune': sc(IMM), 'stromal': sc(STR), 'subtype': sc(up) - sc(dn), 'cna_burden': cna_burden, 'log_tmb': np.log10(burden.reindex(samp).values + 1)}, index=samp)
    out = [D] + ([D_raw] if return_raw else []) + [g.chr.values, samp, cov] + ([g.iloc[:, 0].values] if return_gids else [])
    return tuple(out)
def scores(D, chrs, lags_cis=(1, 2, 3), lags_far=(20, 30)):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    def lagmean(L):
        parts = [Z[idx[:-L]] * Z[idx[L:]] for c in CHR for idx in [np.where(chrs == c)[0]] if len(idx) > L + 5]
        return np.concatenate(parts, 0).mean(0)
    return np.mean([lagmean(L) for L in lags_cis], 0) - np.mean([lagmean(L) for L in lags_far], 0)
