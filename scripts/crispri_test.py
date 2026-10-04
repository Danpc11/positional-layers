"""Does CRISPRi silencing of a gene propagate to its cis-coupled neighbours? (Replogle 2022 K562 genome-wide Perturb-seq)
Response: z-normalised pseudobulk (relative to non-targeting controls), robustly rescaled per perturbation.
Baseline coupling: correlation across 615 BeatAML2 leukaemias (GC-corrected), for target-gene pairs.
Cis: measured genes within 1 Mb of the target. Trans control: genes on other chromosomes with the same coupling range."""
import numpy as np, pandas as pd, anndata as ad, pyannotables as pa, statsmodels.formula.api as smf
CHR = [str(i) for i in range(1, 23)] + ['X']
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin(CHR)]
G['chr'] = G.Chromosome.astype(str); G['tss'] = np.where(G.Strand.astype(str) == '+', G.Start, G.End)
A = ad.read_h5ad('/mnt/user-data/uploads/K562_gwps_normalized_bulk_01.h5ad'); X = np.asarray(A.X, dtype=np.float32)
obs = A.obs.copy(); obs['target'] = obs.index.str.split('_').str[-1]; genes = np.array(A.var.index)
keep = (obs.fold_expr < 0.5) & obs.target.isin(G.index) & ~obs.index.str.contains('non-targeting') & (obs.num_cells_filtered >= 25)
obs = obs[keep]; X = X[keep.values]; print('effective perturbations (>=50% knockdown):', len(obs))
med = np.median(X, 1, keepdims=True); mad = np.median(np.abs(X - med), 1, keepdims=True) * 1.4826 + 1e-6; R = (X - med) / mad    # robust per-perturbation scaling
# baseline coupling from BeatAML2
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
B = pd.read_csv('/home/claude/beataml/beataml_waves1to4_counts_dbgap.txt', sep='\t', low_memory=False); B['gid'] = B.stable_id.str.split('.').str[0]
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
D = pd.concat(rows, ignore_index=True); D.to_csv('crispri_pairs.csv.gz', index=False)
D['rk'] = D.r * D.kd; C = D[D.type == 'cis'].copy(); T = D[D.type == 'trans']
C['dbin'] = pd.cut(C.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['<10kb', '10-50kb', '50-200kb', '200-500kb', '0.5-1Mb'])
print(f'cis pairs {len(C):,} | trans pairs {len(T):,}')
mc = smf.ols('resp ~ rk + r + kd + C(dbin) + C(dbin):kd', data=C).fit(cov_type='cluster', cov_kwds={'groups': C.pert.astype('category').cat.codes})
mt = smf.ols('resp ~ rk + r + kd', data=T).fit(cov_type='cluster', cov_kwds={'groups': T.pert.astype('category').cat.codes})
print(f"cis:   coupling x knockdown {mc.params['rk']:.3f} (95% CI {mc.conf_int().loc['rk', 0]:.3f} to {mc.conf_int().loc['rk', 1]:.3f}), P = {mc.pvalues['rk']:.1e}")
print(f"trans: coupling x knockdown {mt.params['rk']:.3f} (95% CI {mt.conf_int().loc['rk', 0]:.3f} to {mt.conf_int().loc['rk', 1]:.3f}), P = {mt.pvalues['rk']:.1e}")
print('\nmean response of cis neighbours by distance (all) and by coupling tertile, strong knockdowns (kd>2):')
S = C[C.kd > 2].copy(); S['rt'] = pd.qcut(S.r, 3, labels=['low r', 'mid r', 'high r'])
print(S.pivot_table(index='dbin', columns='rt', values='resp', aggfunc='mean', observed=True).round(3).to_string())
St = T[T.kd > 2].copy(); St['rt'] = pd.cut(St.r, [-1, S.r.quantile(1 / 3), S.r.quantile(2 / 3), 1], labels=['low r', 'mid r', 'high r'])
print('trans genes, same coupling cut-offs:', St.groupby('rt', observed=True).resp.mean().round(3).to_dict())
