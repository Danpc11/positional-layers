"""Contact exponents in 40-kb Hi-C maps of human tissues (Schmitt et al. 2016; Supplementary Note 1, section 6).
Gene-pair tables are made from the Fit-Hi-C files of GSE87112 with reduce_fithic.py (DATA/schmitt/schmitt_gene_pairs/).

Predictions fixed in advance (as established in GM12878 and IMR-90 at 25 kb): k_across = 1/2, k_within = 1,
k_within / k_across = 2, with the same values in every tissue.
k_across: log mean coupling on log mean observed contact across distance bins (80 kb - 2 Mb).
k_within: inside each distance bin, quintiles of Fit-Hi-C O/E; log mean coupling on log mean observed contact, with
distance-bin fixed effects.
Coupling: GTEx GC- and technically-corrected correlations of protein-coding gene pairs, same pipeline as hic_test.py.
Bin biases: Fit-Hi-C observed counts are not normalised, so per-bin biases (coverage, mappability) are removed by
matrix balancing on the available bin pairs (as KR/ICE): b_i = sum_j obs_ij / sum_j exp_ij b_j, iterated per chromosome;
contact = obs / (b_i b_j), O/E = contact / expected. Calibration rule: GM12878 and IMR-90 must reproduce the Rao et al.
25-kb KR estimates (k_within 0.95 and 0.78; k_across 0.51 and 0.52) before tissues are interpreted.
Inclusion rule fixed in advance: median observed contact >= 5 at 80 kb (low-depth maps reported but flagged).
95% intervals: 100 resamples of 10-Mb genomic blocks. Usage: python contact_exponent_schmitt.py AD,AO,...
Output: OUTDIR/contact_exponent_schmitt.csv (one row per tissue; re-runs replace rows)
"""
import sys, numpy as np, pandas as pd
from poslayers.config import DATA, OUTDIR, upsert_csv
from lib_boot import blocks_for
import pyannotables as pa

G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()]
TSSBIN = pd.Series(np.where(G37.Strand.astype(str).isin(['1', '+']), G37.Start, G37.End) // 40000, index=G37.index)
def balance(H, iters=30):
    H = H.assign(b1=TSSBIN.reindex(H.g1).values, b2=TSSBIN.reindex(H.g2).values); out = np.full(len(H), np.nan)
    for ch, d in H.groupby('chr'):
        u = d.drop_duplicates(['b1', 'b2']); bins = np.unique(np.concatenate([u.b1, u.b2])); ix = {b: i for i, b in enumerate(bins)}
        i, j = u.b1.map(ix).values, u.b2.map(ix).values; o, e = u.observed.values.astype(float), u.expected.values.astype(float); bvec = np.ones(len(bins))
        for _ in range(iters):
            num = np.bincount(i, o, len(bins)) + np.bincount(j, o, len(bins)); den = np.bincount(i, e * bvec[j], len(bins)) + np.bincount(j, e * bvec[i], len(bins))
            bvec = np.where(den > 0, num / np.maximum(den, 1e-12), 1.0); bvec = np.where(bvec > 0, bvec, np.nan); bvec /= np.nanmean(bvec)
        out[d.index.values] = bvec[d.b1.map(ix).values] * bvec[d.b2.map(ix).values]
    return out

GTEX = {'AD': 'adrenal_gland', 'AO': 'artery_aorta', 'BL': 'bladder', 'CO': 'brain_frontal_cortex_ba9', 'HC': 'brain_hippocampus', 'LG': 'lung',
        'LI': 'liver', 'LV': 'heart_left_ventricle', 'PA': 'pancreas', 'PO': 'muscle_skeletal', 'RV': 'heart_left_ventricle', 'SB': 'small_intestine_terminal_ileum',
        'SX': 'spleen', 'GM12878': 'cells_ebv-transformed_lymphocytes', 'imr90': 'cells_cultured_fibroblasts'}
DB = [1.5, 3.5, 5.5, 9.5, 15.5, 25.5, 35.5, 50.5]
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')

def zmat(tis):
    C = pd.read_csv(DATA + f'gtex/gene_reads_adult_gtex_v11_{tis}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
    samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(BM.index)]; C = C[C.median(axis=1) >= 10]
    Y = np.log2(C.values / C.values.sum(0) * 1e6 + 1); D = Y - Y.mean(1, keepdims=True); gz = BM.loc[C.index, 'gc'].values; gz = (gz - gz.mean()) / gz.std()
    G1 = np.column_stack([np.ones(len(gz)), gz, gz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
    a = SA.loc[samp]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
    T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
    D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T; return (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9), C.index, len(samp)
def k_pair(d):
    a = d.groupby('db').agg(r=('r', 'mean'), c=('observed', 'mean')); a = a[(a.r > 0) & (a.c > 0)]; ka = np.polyfit(np.log(a.c), np.log(a.r), 1)[0]
    g = d.groupby(['db', 'q']).agg(r=('r', 'mean'), c=('observed', 'mean')).reset_index(); g = g[(g.r > 0) & (g.c > 0)]
    kw = np.linalg.lstsq(np.column_stack([np.log(g.c), pd.get_dummies(g.db).values.astype(float)]), np.log(g.r), rcond=None)[0][0]; return kw, ka
if __name__ == '__main__':
    for code in sys.argv[1].split(','):
        H = pd.read_csv(DATA + f'schmitt/schmitt_gene_pairs/{code}_gene_pairs.csv.gz', dtype={'chr': str})
        depth = float(H[H.bin_distance == 2].observed.median())
        Z, genes, ns = zmat(GTEX[code]); pos = {g: i for i, g in enumerate(genes)}
        H = H[(H.bin_distance >= 2) & H.g1.isin(pos) & H.g2.isin(pos) & H.chr.isin([str(i) for i in range(1, 23)])].reset_index(drop=True)
        H['r'] = (Z[H.g1.map(pos).values] * Z[H.g2.map(pos).values]).mean(1)
        bb = balance(H); H = H[np.isfinite(bb) & (bb > 0)].copy(); H['observed'] = H.observed / bb[np.isfinite(bb) & (bb > 0)]; H['oe'] = H.observed / H.expected; H = H.reset_index(drop=True)
        H['db'] = pd.cut(H.bin_distance, DB, labels=False); H['q'] = H.groupby('db').oe.rank(pct=True, method='first').mul(5).clip(upper=4.999).astype(int)
        kw, ka = k_pair(H); blk = blocks_for(H.g1.values); grp = [np.where(blk == b)[0] for b in np.unique(blk)]; rng = np.random.default_rng(6)
        bs = np.array([k_pair(H.iloc[np.concatenate([grp[i] for i in rng.integers(0, len(grp), len(grp))])]) for _ in range(100)]); ratio = bs[:, 0] / bs[:, 1]
        row = {'code': code, 'gtex_tissue': GTEX[code], 'samples': ns, 'pairs': len(H), 'median_contact_80kb': depth, 'included': depth >= 5,
               'k_within': kw, 'kw_ci_low': np.percentile(bs[:, 0], 2.5), 'kw_ci_high': np.percentile(bs[:, 0], 97.5),
               'k_across': ka, 'ka_ci_low': np.percentile(bs[:, 1], 2.5), 'ka_ci_high': np.percentile(bs[:, 1], 97.5),
               'ratio': kw / ka, 'ratio_ci_low': np.percentile(ratio, 2.5), 'ratio_ci_high': np.percentile(ratio, 97.5)}
        upsert_csv(pd.DataFrame([row]), OUTDIR + 'contact_exponent_schmitt.csv', ['code'])
        print(f"{code:8s} n={ns:4d} depth {depth:5.0f}  k_within {kw:.2f} [{row['kw_ci_low']:.2f},{row['kw_ci_high']:.2f}]  k_across {ka:.2f} [{row['ka_ci_low']:.2f},{row['ka_ci_high']:.2f}]  ratio {kw/ka:.2f} [{row['ratio_ci_low']:.2f},{row['ratio_ci_high']:.2f}]", flush=True)
