import numpy as np, pandas as pd, pyannotables as pa
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
mu = pd.read_csv('GDC-PANCAN.mutect2_snv.tsv', sep='\t', usecols=['Sample_ID', 'gene', 'effect', 'filter']); mu = mu[mu['filter'] == 'PASS']
coding = mu[~mu.effect.isin(['synonymous_variant', 'intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant', 'non_coding_transcript_exon_variant'])]
burden = mu.groupby('Sample_ID').size()
LAGS = np.array([1, 2, 3, 5, 10, 20, 30])
def prepare(c):
    z = np.load(f'expr_{c}.npz', allow_pickle=True); X = z['X']; genes = z['genes']; samp = z['samples']
    keep = np.array([s[13:15] == '01' for s in samp]) & pd.Index(samp).isin(burden.index); X = X[:, keep]; samp = samp[keep]
    C = np.clip(2 ** X - 1, 0, None); g = pd.DataFrame({'gid': genes}).set_index('gid').join(G38, how='inner'); ok = pd.Index(genes).isin(g.index)
    C = C[ok]; genes = genes[ok]; good = np.median(C, 1) >= 10; C = C[good]; genes = genes[good]
    g = G38.loc[genes].reset_index(); o = np.lexsort((g.start.values, g.chr.values)); C = C[o]; g = g.iloc[o].reset_index(drop=True)
    Y = np.log2(C / C.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
    if CNADJ:
        hdr = pd.read_csv('GDC-PANCAN.gistic.tsv', sep='\t', nrows=0).columns; use = [hdr[0]] + [x for x in samp if x in set(hdr)]
        CN = pd.read_csv('GDC-PANCAN.gistic.tsv', sep='\t', usecols=use, index_col=0); CN.index = CN.index.str.split('.').str[0]
        CN = CN.reindex(index=g.iloc[:, 0].values, columns=samp)
        have = (CN.notna().mean(axis=0) > 0.5).values; print('samples with copy number:', int(have.sum()), 'of', len(have), '| genes with copy number:', int(CN.notna().all(axis=1).sum()), flush=True); D = D[:, have]; samp = samp[have]; CNv = np.nan_to_num(CN.values[:, have].astype(float), nan=0.0)   # missing copy-number calls treated as neutral
        for i in range(D.shape[0]):                         # remove each gene's own dosage effect (categorical copy-number state)
            cn = CNv[i]
            for lv in np.unique(cn): D[i, cn == lv] -= D[i, cn == lv].mean()
    return D, g.chr.values, samp
def lagprof(D, chrs, cols):
    Z = D[:, cols]; Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9); out = []
    for L in LAGS:
        v = [np.mean(Z[idx[:-L]] * Z[idx[L:]], axis=1) for c in CHR for idx in [np.where(chrs == c)[0]] if len(idx) > L + 5]
        out.append(np.mean(np.concatenate(v)))
    return np.array(out)
def test(c, genes_mut, hyper_q=0.9, nperm=200, seed=0, truncating=False):
    D, chrs, samp = prepare(c); pat = np.array([s[:16] for s in samp])
    b = burden.reindex(samp).values; nonhyper = b <= np.quantile(b, hyper_q)
    src = coding[coding.effect.str.contains('stop_gained|frameshift|splice_acceptor|splice_donor', regex=True)] if truncating else coding
    mut = pd.Index(samp).isin(set(src[src.gene.isin(genes_mut)].Sample_ID))
    idx = np.where(nonhyper)[0]; m = mut[idx]; im, iw = idx[m], idx[~m]
    rng = np.random.default_rng(seed); k = len(im)
    pm = lagprof(D, chrs, im); pw = np.mean([lagprof(D, chrs, rng.choice(iw, k, replace=False)) for _ in range(20)], 0)
    floor = np.mean([lagprof(D, chrs, rng.choice(idx, k, replace=False)) for _ in range(1)], 0)
    cisx = lambda p: p[[0, 1, 2]].mean() - p[[5, 6]].mean()          # cis excess over the long-range (copy-number/purity) floor
    null = []
    for _ in range(nperm):
        sh = rng.permutation(idx); null.append(cisx(lagprof(D, chrs, sh[:k])) - cisx(lagprof(D, chrs, rng.choice(sh[k:], k, replace=False))))
    obs = cisx(pm) - cisx(pw); null = np.array(null)
    return {'cohort': c, 'cn_adjusted': CNADJ, 'truncating_only': truncating, 'hypermutation_quantile': hyper_q, 'genes': '/'.join(genes_mut), 'n_mutant': k, 'n_wildtype': len(iw), 'excluded_hypermutated': int((~nonhyper).sum()),
            **{f'mut_L{L}': v for L, v in zip(LAGS, pm)}, **{f'wt_L{L}': v for L, v in zip(LAGS, pw)},
            'cis_excess_mut': cisx(pm), 'cis_excess_wt': cisx(pw), 'diff_cis_excess': obs, 'perm_p_two_sided': (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)}
import sys
CNADJ = len(sys.argv) > 4 and sys.argv[4] == 'cn'
c, gm = sys.argv[1], sys.argv[2].split(','); r = test(c, gm, nperm=int(sys.argv[3]) if len(sys.argv) > 3 else 100, truncating=len(sys.argv) > 5 and sys.argv[5] == 'trunc', hyper_q=float(sys.argv[6]) if len(sys.argv) > 6 else 0.9)
pd.DataFrame([r]).to_csv('cohesin_results.csv', mode='a', header=not __import__('os').path.exists('cohesin_results.csv'), index=False)
print(pd.Series(r).to_string())
