"""Per-tumour cis-excess score and covariate-adjusted test of cohesin/CTCF loss.
Score: for each tumour, mean product of pooled-standardised residual expression of genes 1-3 positions apart minus that of
genes 20-30 apart (GC- and own-copy-number-corrected expression). Covariates: expression-based immune and stromal scores
(purity proxy), copy-number burden, log mutation burden, cohort-specific subtype score."""
import sys, os, numpy as np, pandas as pd, pyannotables as pa, statsmodels.api as sm
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene name': 'sym'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc', 'sym']], how='inner')
mu = pd.read_csv('GDC-PANCAN.mutect2_snv.tsv', sep='\t', usecols=['Sample_ID', 'gene', 'effect', 'filter']); mu = mu[mu['filter'] == 'PASS']
NONCOD = ['synonymous_variant', 'intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant', 'non_coding_transcript_exon_variant']
coding = mu[~mu.effect.isin(NONCOD)]; trunc = coding[coding.effect.str.contains('stop_gained|frameshift|splice_acceptor|splice_donor', regex=True)]
burden = mu.groupby('Sample_ID').size()
IMM = ['PTPRC', 'CD2', 'CD3E', 'CD3D', 'CD247', 'LCK', 'CD48', 'CD53', 'CD52', 'CORO1A', 'LAPTM5', 'CD37', 'FCER1G', 'TYROBP', 'CD74', 'HLA-DRA', 'CCL5', 'IL2RG', 'CXCR4', 'SELL']
STR = ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A3', 'DCN', 'LUM', 'FBN1', 'FAP', 'PDGFRB', 'THY1', 'SPARC', 'VCAN', 'MMP2', 'CDH11', 'FN1', 'POSTN', 'COL11A1', 'SFRP2', 'ACTA2']
SUB = {'BLCA': (['GATA3', 'KRT20', 'UPK1B', 'UPK2', 'FOXA1', 'PPARG', 'FGFR3'], ['KRT5', 'KRT6A', 'KRT14', 'CD44', 'CDH3']),
       'UCEC': (['PGR', 'ESR1', 'PAX8', 'MSX1'], ['CDKN2A', 'L1CAM', 'WT1']), 'STAD': (['CDX2', 'MUC2', 'TFF3'], ['VIM', 'ZEB1']), 'COAD': (['CDX2', 'VIL1', 'LGALS4'], ['VIM', 'ZEB1'])}
def prepare(c):
    z = np.load(f'expr_{c}.npz', allow_pickle=True); X = z['X']; genes = z['genes']; samp = z['samples']
    keep = np.array([s[13:15] == '01' for s in samp]) & pd.Index(samp).isin(burden.index); X = X[:, keep]; samp = samp[keep]
    C = np.clip(2 ** X - 1, 0, None); ok = pd.Index(genes).isin(G38.index); C = C[ok]; genes = genes[ok]
    good = np.median(C, 1) >= 10; C = C[good]; genes = genes[good]
    g = G38.loc[genes].reset_index(); o = np.lexsort((g.start.values, g.chr.values)); C = C[o]; g = g.iloc[o].reset_index(drop=True)
    Y = np.log2(C / C.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True)
    gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
    hdr = pd.read_csv('GDC-PANCAN.gistic.tsv', sep='\t', nrows=0).columns
    CN = pd.read_csv('GDC-PANCAN.gistic.tsv', sep='\t', usecols=[hdr[0]] + [x for x in samp if x in set(hdr)], index_col=0); CN.index = CN.index.str.split('.').str[0]
    CN = CN.reindex(index=g.iloc[:, 0].values, columns=samp); have = (CN.notna().mean(axis=0) > 0.5).values
    D = D[:, have]; Y = Y[:, have]; samp = samp[have]; CNv = np.nan_to_num(CN.values[:, have].astype(float), nan=0.0)
    for i in range(D.shape[0]):
        cn = CNv[i]
        for lv in np.unique(cn): D[i, cn == lv] -= D[i, cn == lv].mean()
    sym = g.sym.astype(str).values; sc = lambda gs: np.nanmean([(Y[sym == s][0] - Y[sym == s][0].mean()) / (Y[sym == s][0].std() + 1e-9) for s in gs if (sym == s).any()], 0)
    up, dn = SUB[c]
    cov = pd.DataFrame({'immune': sc(IMM), 'stromal': sc(STR), 'subtype': sc(up) - sc(dn), 'cna_burden': np.mean(CNv != 0, 0), 'log_tmb': np.log10(burden.reindex(samp).values + 1)}, index=samp)
    return D, g.chr.values, samp, cov
def scores(D, chrs, lags_cis=(1, 2, 3), lags_far=(20, 30)):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    def lagmean(L):
        parts = [Z[idx[:-L]] * Z[idx[L:]] for c in CHR for idx in [np.where(chrs == c)[0]] if len(idx) > L + 5]
        return np.concatenate(parts, 0).mean(0)
    return np.mean([lagmean(L) for L in lags_cis], 0) - np.mean([lagmean(L) for L in lags_far], 0)
rows = []
for c, genes_mut in [('BLCA', ['STAG2']), ('UCEC', ['CTCF']), ('UCEC', ['STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']), ('STAD', ['CTCF', 'STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']), ('COAD', ['CTCF', 'STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL'])]:
    D, chrs, samp, cov = prepare(c); e = scores(D, chrs)
    for mclass, src in [('all coding', coding), ('truncating', trunc)]:
        for hq in (0.9, 0.7):
            b = cov.log_tmb.values; keep = b <= np.quantile(b, hq)
            mut = pd.Index(samp).isin(set(src[src.gene.isin(genes_mut)].Sample_ID)).astype(float)
            y = e[keep]; m = mut[keep]; Xc = cov[keep]
            if m.sum() < 10: continue
            Xd = sm.add_constant(pd.DataFrame({'mutant': m}, index=Xc.index).join((Xc - Xc.mean()) / Xc.std()))
            fit = sm.OLS(y, Xd).fit(cov_type='HC3'); base = y[m == 0].mean()
            # permutation of the mutant label on covariate-residualised scores
            r = y - sm.OLS(y, sm.add_constant(((Xc - Xc.mean()) / Xc.std()).values)).fit().fittedvalues
            obs = r[m == 1].mean() - r[m == 0].mean(); rng = np.random.default_rng(0)
            null = np.array([(lambda p: r[p == 1].mean() - r[p == 0].mean())(rng.permutation(m)) for _ in range(10000)])
            rows.append({'cohort': c, 'genes': '/'.join(genes_mut) if len(genes_mut) < 3 else 'cohesin+CTCF' if 'CTCF' in genes_mut else 'cohesin', 'mutation_class': mclass, 'tmb_quantile_kept': hq,
                         'n_mutant': int(m.sum()), 'n_wildtype': int((m == 0).sum()), 'mean_cis_excess_wt': base,
                         'raw_diff_pct': 100 * (y[m == 1].mean() - base) / base, 'adj_effect': fit.params['mutant'], 'adj_effect_pct': 100 * fit.params['mutant'] / base,
                         'ci_low_pct': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high_pct': 100 * fit.conf_int().loc['mutant', 1] / base,
                         'p_adj_HC3': fit.pvalues['mutant'], 'p_perm_10000': (np.sum(np.abs(null) >= abs(obs)) + 1) / 10001})
    print(c, '/'.join(genes_mut)[:20], 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv('cohesin_v2_results.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
