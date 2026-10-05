"""Does CRISPRi silencing of a gene propagate to its cis-coupled neighbours? (Replogle 2022 K562 genome-wide Perturb-seq)
Response: z-normalised pseudobulk (relative to non-targeting controls), robustly rescaled per perturbation.
Baseline coupling: correlation across 615 BeatAML2 leukaemias (GC-corrected), for target-gene pairs.
Cis: measured genes within 1 Mb of the target. Trans control: up to 60 genes drawn at random from other chromosomes for
each perturbation; they are NOT matched on coupling (crispri_analyse.py reports both coupling distributions).
This script only builds the pair table; every model is fitted in crispri_analyse.py."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd, anndata as ad, pyannotables as pa, statsmodels.formula.api as smf
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
G['chr'] = G.Chromosome.astype(str); G['tss'] = np.where(G.Strand.astype(str) == '+', G.Start, G.End)
A = ad.read_h5ad(DATA + 'perturb/K562_gwps_normalized_bulk_01.h5ad'); X = np.asarray(A.X, dtype=np.float32)
obs = A.obs.copy(); obs['target'] = obs.index.str.split('_').str[-1]; genes = np.array(A.var.index)
keep = (obs.fold_expr < 0.5) & obs.target.isin(G.index) & ~obs.index.str.contains('non-targeting') & (obs.num_cells_filtered >= 25)
obs = obs[keep]; X = X[keep.values]; print('effective perturbations (>=50% knockdown):', len(obs))
med = np.median(X, 1, keepdims=True); mad = np.median(np.abs(X - med), 1, keepdims=True) * 1.4826 + 1e-6; R = (X - med) / mad    # robust per-perturbation scaling
# baseline coupling from BeatAML2
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
B = pd.read_csv(DATA + 'beataml/beataml_waves1to4_counts_dbgap.txt', sep='\t', low_memory=False); B['gid'] = B.stable_id.str.split('.').str[0]
B = B.drop_duplicates('gid').set_index('gid')[[x for x in B.columns if x.startswith('BA')]]; B = B.loc[B.index.intersection(BM.index)]; B = B[B.median(axis=1) >= 5]
Yb = np.log2(B.values / B.values.sum(0) * 1e6 + 1); Db = Yb - Yb.mean(1, keepdims=True); gcz = BM.loc[B.index, 'gc'].values; gcz = (gcz - gcz.mean()) / gcz.std()
G1 = np.column_stack([np.ones(len(gcz)), gcz, gcz ** 2]); Db = Db - G1 @ np.linalg.lstsq(G1, Db, rcond=None)[0]; Zb = (Db - Db.mean(1, keepdims=True)) / (Db.std(1, keepdims=True) + 1e-9)
bpos = {g: i for i, g in enumerate(B.index)}; vpos = {g: i for i, g in enumerate(genes)}
meas = [g for g in genes if g in G.index and g in bpos]; mg = G.loc[meas]
rng = np.random.default_rng(0); rows = []
for k, (pid, o) in enumerate(obs.iterrows()):
    t = o.target
    if t not in bpos: continue
    tc, tt = G.at[t, 'chr'], G.at[t, 'tss']; kd = -np.log2(max(o.fold_expr, 0.01))
    same = mg[(mg.chr == tc) & ((mg.tss - tt).abs() <= 1e6) & (mg.index != t)]
    other = mg[mg.chr != tc]; other = other.iloc[rng.choice(len(other), min(60, len(other)), replace=False)]
    for typ, sub in (('cis', same), ('trans', other)):
        if not len(sub): continue
        r = (Zb[[bpos[g] for g in sub.index]] * Zb[bpos[t]]).mean(1); resp = R[k, [vpos[g] for g in sub.index]]
        d = (sub.tss - tt).abs().values if typ == 'cis' else np.full(len(sub), np.nan)
        rows.append(pd.DataFrame({'pert': pid, 'target': t, 'type': typ, 'kd': kd, 'r': r, 'resp': resp, 'dist': d}))
D = pd.concat(rows, ignore_index=True); D.to_csv(OUTDIR + 'crispri_pairs.csv.gz', index=False)
print(f"pairs written: {(D.type == 'cis').sum():,} cis, {(D.type == 'trans').sum():,} trans")
