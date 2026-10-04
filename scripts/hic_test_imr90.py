import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
"""Distance vs 3D contact (GM12878 in situ Hi-C, KR, 25 kb) as predictors of cis coupling in GTEx EBV lymphocytes."""
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.formula.api as smf
from scipy import stats
CHR = [str(i) for i in range(1, 23)] + ['X']
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
C = pd.read_csv(DATA + 'gtex/gene_reads_adult_gtex_v11_cells_cultured_fibroblasts_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description')
C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]; samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]
C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True); gc = BM.loc[C.index, 'gc'].values; gcz = (gc - gc.mean()) / gc.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T
Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9); pos = {g: i for i, g in enumerate(C.index)}
H = pd.read_csv('/mnt/user-data/uploads/hic_contacts_IMR90_csv.gz', dtype={'chr': str}); E = pd.read_csv('/mnt/user-data/uploads/hic_expected_IMR90_csv.gz', dtype={'chr': str})
H = H[H.g1.isin(pos) & H.g2.isin(pos)].copy()
i1 = H.g1.map(pos).values; i2 = H.g2.map(pos).values; r = np.empty(len(H))
for s in range(0, len(H), 200000): r[s:s + 200000] = (Z[i1[s:s + 200000]] * Z[i2[s:s + 200000]]).mean(1)
H['r'] = r; H = H.merge(E, on=['chr', 'bin_distance'], how='left'); H = H[(H.bin_distance > 0) & (H.expected_KR > 0) & (H.tss_distance > 0)]
H['oe'] = H.contact_KR / H.expected_KR; H['log_oe'] = np.log2(H.oe + 0.05); H['log_d'] = np.log10(H.tss_distance); H['log_c'] = np.log10(H.contact_KR + 0.1)
print(f'gene pairs within 2 Mb with Hi-C and coupling: {len(H):,}')
# decay laws
H['dbin'] = pd.cut(H.tss_distance, [25e3, 50e3, 100e3, 200e3, 500e3, 1e6, 2e6])
dec = H.groupby('dbin', observed=True).agg(n=('r', 'size'), mean_r=('r', 'mean'), mean_contact=('contact_KR', 'mean'), mid=('tss_distance', 'median'))
print(dec.round(4).to_string())
print('log-log slope: coupling %.2f | contact %.2f' % (np.polyfit(np.log10(dec.mid), np.log10(dec.mean_r.clip(lower=1e-4)), 1)[0], np.polyfit(np.log10(dec.mid), np.log10(dec.mean_contact), 1)[0]))
# dissociation: at equal distance, does observed/expected contact predict coupling?
rows = []
for b, d in H.groupby('dbin', observed=True):
    q = pd.qcut(d.oe.rank(method='first'), 5, labels=False)
    rows.append({'distance': str(b), 'pairs': len(d), 'spearman_r_vs_OE': stats.spearmanr(d.oe, d.r)[0], 'r_lowest_OE_quintile': d.r[q == 0].mean(), 'r_highest_OE_quintile': d.r[q == 4].mean()})
print(pd.DataFrame(rows).round(4).to_string(index=False))
m1 = smf.ols('r ~ bs(log_d, df=5)', data=H).fit(); m2 = smf.ols('r ~ bs(log_d, df=5) + log_oe', data=H).fit(); m3 = smf.ols('r ~ bs(log_c, df=5)', data=H).fit(); m4 = smf.ols('r ~ bs(log_c, df=5) + bs(log_d, df=5)', data=H).fit()
print(f"\nR2 distance only {m1.rsquared:.4f} | + O/E contact {m2.rsquared:.4f} (O/E coef {m2.params['log_oe']:.4f}, t {m2.tvalues['log_oe']:.1f})")
print(f"R2 contact only {m3.rsquared:.4f} | contact + distance {m4.rsquared:.4f} | AIC distance {m1.aic:.0f} vs contact {m3.aic:.0f}")
H.drop(columns='dbin').to_csv('hic_coupling_IMR90.csv.gz', index=False)
