"""Replication of the cohesin effect on the cis layer in BeatAML2 (open data). Per-sample cis-excess score as in TCGA:
neighbour products at 1-3 genes minus 20-30 genes, on GC-corrected expression. Covariates: blasts, monocytic score, log TMB, sex,
specimen type. 10,000 label permutations on covariate-residualised scores."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.api as sm
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
G38 = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G38 = G38[~G38.index.duplicated()][['Chromosome', 'Start']]; G38.columns = ['chr', 'start']; G38['chr'] = G38.chr.astype(str); G38 = G38[G38.chr.isin(CHR)].join(BM[['gc']], how='inner')
C = pd.read_csv(DATA + 'beataml/beataml_waves1to4_counts_dbgap.txt', sep='\t'); C['gid'] = C.stable_id.str.split('.').str[0]; sym = dict(zip(C.gid, C.display_label.astype(str)))
C = C.drop_duplicates('gid').set_index('gid')[[x for x in C.columns if x.startswith('BA')]]
mp = pd.read_excel(DATA + 'beataml/beataml_waves1to4_sample_mapping.xlsx'); cl = pd.read_excel(DATA + 'beataml/beataml_wv1to4_clinical.xlsx')
mp = mp[(mp.rna_control != 'yes') & mp.dbgap_rnaseq_sample.notna() & mp.dbgap_dnaseq_sample.notna()] if 'rna_control' in mp else mp
pairs = mp[['dbgap_rnaseq_sample', 'dbgap_dnaseq_sample']].dropna().drop_duplicates('dbgap_rnaseq_sample')
m = pd.read_csv(DATA + 'beataml/beataml_wes_wv1to4_mutations_dbgap.txt', sep='\t', low_memory=False); tmb = m.groupby('dbgap_sample_id').size()
pairs = pairs[pairs.dbgap_rnaseq_sample.isin(C.columns) & pairs.dbgap_dnaseq_sample.isin(tmb.index)]
C = C[pairs.dbgap_rnaseq_sample.values]; C = C.loc[C.index.intersection(G38.index)]; C = C[C.median(axis=1) >= 10]
g = G38.loc[C.index].sort_values(['chr', 'start']); C = C.loc[g.index]; chrs = g.chr.values
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True)
gcz = (g.gc.values - g.gc.mean()) / g.gc.std(); G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
lag = lambda L: np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0)
e = np.mean([lag(L) for L in (1, 2, 3)], 0) - np.mean([lag(L) for L in (20, 30)], 0)
syms = np.array([sym.get(x, '') for x in g.index]); zsc = lambda gs: np.nanmean([(Y[syms == s][0] - Y[syms == s][0].mean()) / (Y[syms == s][0].std() + 1e-9) for s in gs if (syms == s).any()], 0)
dna = pairs.dbgap_dnaseq_sample.values; rna = pairs.dbgap_rnaseq_sample.values
clin = cl.drop_duplicates('dbgap_rnaseq_sample').set_index('dbgap_rnaseq_sample').reindex(rna)
cov = pd.DataFrame({'monocytic': zsc(['CD14', 'LYZ', 'CSF1R', 'FCGR1A', 'CD68']) - zsc(['CD34', 'KIT', 'PROM1']), 'log_tmb': np.log10(tmb.reindex(dna).values + 1),
                    'blasts': pd.to_numeric(clin['%.Blasts.in.BM'], errors='coerce').fillna(pd.to_numeric(clin['%.Blasts.in.PB'], errors='coerce')).values,
                    'male': (clin.consensus_sex.astype(str).str.lower() == 'male').astype(float).values, 'pb_specimen': clin.specimenType.astype(str).str.contains('Peripheral', case=False).astype(float).values}, index=rna)
cov['blasts'] = cov.blasts.fillna(cov.blasts.median())
COH = ['STAG2', 'RAD21', 'SMC1A', 'SMC3', 'STAG1', 'NIPBL']; TR = ['frameshift_variant', 'stop_gained', 'splice_acceptor_variant', 'splice_donor_variant']
rows = []
for name, genes, cls in [('STAG2, any coding', ['STAG2'], None), ('STAG2, truncating', ['STAG2'], TR), ('cohesin, any coding', COH, None), ('cohesin, truncating', COH, TR)]:
    mm = m[m.symbol.isin(genes)]; mm = mm[mm.variant_classification.isin(cls)] if cls else mm
    mut = pd.Index(dna).isin(set(mm.dbgap_sample_id)).astype(float)
    Xs = (cov - cov.mean()) / cov.std(); Xd = sm.add_constant(pd.DataFrame({'mutant': mut}, index=cov.index).join(Xs))
    fit = sm.OLS(e, Xd).fit(cov_type='HC3'); base = e[mut == 0].mean(); r = e - sm.OLS(e, sm.add_constant(Xs.values)).fit().fittedvalues
    red = sm.OLS(e, sm.add_constant(Xs.values)).fit(); Xf = Xd.values; t_obs = fit.tvalues['mutant']; rng = np.random.default_rng(1)
    tnull = np.array([sm.OLS(red.fittedvalues + rng.permutation(red.resid), Xf).fit(cov_type='HC3').tvalues[1] for _ in range(2000)])
    rows.append({'group': name, 'n_samples': len(e), 'n_mutant': int(mut.sum()), 'cis_excess_wt': base, 'raw_pct': 100 * (e[mut == 1].mean() - base) / base,
                 'adj_pct': 100 * fit.params['mutant'] / base, 'ci_low': 100 * fit.conf_int().loc['mutant', 0] / base, 'ci_high': 100 * fit.conf_int().loc['mutant', 1] / base, 'se_pct': 100 * fit.bse['mutant'] / base, 'p_HC3': fit.pvalues['mutant'], 'p_freedman_lane': (np.sum(np.abs(tnull) >= abs(t_obs)) + 1) / 2001})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'beataml_freedman_lane.csv', index=False); pd.set_option('display.width', 200); print(R.round(4).to_string(index=False))
