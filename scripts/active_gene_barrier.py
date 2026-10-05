"""Coupling across an active intervening gene (Fig. 4f; Extended Data Fig. 5a,b; Supplementary Table 9D,E).
Valton et al. 2022: termination sites of active genes stall and unload cohesin.

Prediction fixed in advance: for consecutive genes (i, j, k), the coupling of the flanking genes i and k decreases with the
expression of the intervening gene j, at equal i-k distance and given the expression of i and k.
Secondary: an expression-weighted count of active termination sites lying between the TSSs of i and k
(i on the + strand, j always, k on the - strand).
Model: r_ik ~ expr_j + expr_i + expr_k + distance-bin fixed effects (expressions standardised; expression from donor half A,
coupling from half B, to avoid shared sampling noise). 95% intervals: 200 resamples of 10-Mb genomic blocks.
Usage: python active_gene_barrier.py TISSUE
Outputs: OUTDIR/active_gene_barrier_<tissue>.csv (coefficients), OUTDIR/active_gene_barrier_means_<tissue>.csv (adjusted means by quartile)
"""
import sys, numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA, OUTDIR
from lib_boot import blocks_for

tis = sys.argv[1]
BM = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc'}).drop_duplicates('gid').set_index('gid')
SA = pd.read_csv(DATA + 'GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt', sep='\t', low_memory=False, usecols=['SAMPID', 'SMRIN', 'SMTSISCH', 'SMNABTCH', 'SMGEBTCH']).set_index('SAMPID')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])]
C = pd.read_csv(DATA + f'gtex/gene_reads_adult_gtex_v11_{tis}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description'); C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]
samp = [x for x in C.columns if x in SA.index and pd.notna(SA.loc[x, 'SMRIN'])]; C = C[samp]; C = C.loc[C.index.intersection(BM.index).intersection(G.index)]
lib = C.values.sum(0); cpm_all = np.log2(C.values / lib * 1e6 + 1)
rng = np.random.default_rng(12); donors = np.array(['-'.join(s.split('-')[:2]) for s in samp]); u = np.unique(donors); hA = np.isin(donors, rng.choice(u, len(u) // 2, replace=False))
level = pd.Series(cpm_all[:, hA].mean(1), index=C.index)                         # expression of every gene, half A
E = C[C.median(axis=1) >= 10]                                                     # genes with measurable coupling
g = G.loc[E.index].sort_values(['Chromosome', 'Start']); E = E.loc[g.index]
Y = np.log2(E.values / lib * 1e6 + 1)[:, ~hA]; D = Y - Y.mean(1, keepdims=True); gz = BM.loc[E.index, 'gc'].values; gz = (gz - gz.mean()) / gz.std()
G1 = np.column_stack([np.ones(len(gz)), gz, gz ** 2]); D = D - G1 @ np.linalg.lstsq(G1, D, rcond=None)[0]
a = SA.loc[np.array(samp)[~hA]]; rin = a.SMRIN.astype(float).values; isch = pd.to_numeric(a.SMTSISCH, errors='coerce'); isch = isch.fillna(isch.median()).fillna(0).values
T = np.column_stack([np.ones(len(rin)), rin, rin ** 2, isch, pd.get_dummies(a.SMNABTCH.astype(str), drop_first=True).values.astype(float), pd.get_dummies(a.SMGEBTCH.astype(str), drop_first=True).values.astype(float)])
D = (D.T - T @ np.linalg.lstsq(T, D.T, rcond=None)[0]).T; Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
# triples of consecutive genes among ALL annotated genes (the intervening gene j may itself be lowly expressed)
A = G.sort_values(['Chromosome', 'Start']); A = A[A.index.isin(level.index)]; pos = {x: i for i, x in enumerate(E.index)}
tss = pd.Series(np.where(A.Strand.astype(str).isin(['1', '+']), A.Start, A.End), index=A.index); plus = A.Strand.astype(str).isin(['1', '+'])
rows = []
for ch, d in A.groupby('Chromosome', observed=True):
    ids = d.index.values
    for n in range(len(ids) - 2):
        i, j, k = ids[n], ids[n + 1], ids[n + 2]
        if i in pos and k in pos:
            barrier = level[i] * plus[i] + level[j] + level[k] * (not plus[k])
            rows.append((i, k, j, abs(tss[k] - tss[i]), level[i], level[j], level[k], barrier, float(plus[i]), float(not plus[k]), float(A.loc[j, 'End'] - A.loc[j, 'Start'])))
P = pd.DataFrame(rows, columns=['gi', 'gk', 'gj', 'dist', 'ei', 'ej', 'ek', 'barrier', 'i_tts_between', 'k_tts_between', 'len_j'])
P['r'] = (Z[P.gi.map(pos).values] * Z[P.gk.map(pos).values]).mean(1)
P['db'] = pd.qcut(P.dist, 8, labels=False)
for c in ('ei', 'ej', 'ek', 'barrier'): P[c + '_z'] = (P[c] - P[c].mean()) / P[c].std()
def fit(d, x):
    X = np.column_stack([d[x].values, d.ei_z.values, d.ek_z.values, pd.get_dummies(d.db).values.astype(float)]); return np.linalg.lstsq(X, d.r.values, rcond=None)[0][0]
blk = blocks_for(P.gi.values); grp = [np.where(blk == b)[0] for b in np.unique(blk)]; out = []
for x, lab in (('ej_z', 'expression of the intervening gene (primary)'), ('barrier_z', 'active termination sites between i and k (secondary)')):
    est = fit(P, x); r_ = np.random.default_rng(1); bs = [fit(P.iloc[np.concatenate([grp[q] for q in r_.integers(0, len(grp), len(grp))])], x) for _ in range(200)]
    out.append({'tissue': tis, 'triples': len(P), 'term': lab, 'coef_per_SD': est, 'ci_low': np.percentile(bs, 2.5), 'ci_high': np.percentile(bs, 97.5), 'mean_coupling': P.r.mean()})
    print(f"{tis}: {lab}: {est:+.4f} per SD [{np.percentile(bs, 2.5):+.4f}, {np.percentile(bs, 97.5):+.4f}]  (mean r_ik {P.r.mean():.4f}, {len(P):,} triples)", flush=True)
# robustness: continuous log distance, log length of j, and orientation of the flanks (barrier vs promoter competition)
P['logd'] = np.log(P.dist); P['loglen_j'] = np.log(P.len_j.clip(lower=100)); P['i_tts_x_ei'] = P.i_tts_between * P.ei_z; P['k_tts_x_ek'] = P.k_tts_between * P.ek_z
def fit2(d):
    cols = ['ej_z', 'ei_z', 'ek_z', 'logd', 'loglen_j', 'i_tts_between', 'k_tts_between', 'i_tts_x_ei', 'k_tts_x_ek']
    X = np.column_stack([d[c].values for c in cols] + [pd.get_dummies(d.db).values.astype(float)]); return np.linalg.lstsq(X, d.r.values, rcond=None)[0][:len(cols)]
est = fit2(P); r_ = np.random.default_rng(2); bs = np.array([fit2(P.iloc[np.concatenate([grp[q] for q in r_.integers(0, len(grp), len(grp))])]) for _ in range(200)])
names = ['expression of j', 'expression of i', 'expression of k', 'log distance', 'log length of j', 'i termination site between', 'k termination site between', 'i TTS between x expression of i', 'k TTS between x expression of k']
for n_, e_, lo_, hi_ in zip(names, est, np.percentile(bs, 2.5, 0), np.percentile(bs, 97.5, 0)):
    out.append({'tissue': tis, 'triples': len(P), 'term': 'robust: ' + n_, 'coef_per_SD': e_, 'ci_low': lo_, 'ci_high': hi_, 'mean_coupling': P.r.mean()})
    print(f'   robust {n_:36s} {e_:+.4f} [{lo_:+.4f}, {hi_:+.4f}]', flush=True)
pd.DataFrame(out).to_csv(OUTDIR + f'active_gene_barrier_{tis}.csv', index=False)
q = pd.qcut(P.ej, 4, labels=['j lowest', 'j low', 'j high', 'j highest']); print(P.groupby(q, observed=True).r.mean().round(4).to_dict())

# ---- adjusted means by expression quartile of the intervening gene (Fig. 4f): residual within distance x flank-expression strata
P['jq'] = pd.qcut(P.ej, 4, labels=False)
for c in ('ei', 'ek'): P[c + 'q'] = pd.qcut(P[c].rank(method='first'), 3, labels=False)
P['stratum'] = P.db.astype(str) + '_' + P.eiq.astype(str) + '_' + P.ekq.astype(str)
def adjusted(d): return (d.r - d.groupby('stratum').r.transform('mean') + d.r.mean()).groupby(d.jq).mean()
est = adjusted(P); r_ = np.random.default_rng(23); bs = pd.DataFrame([adjusted(P.iloc[np.concatenate([grp[q] for q in r_.integers(0, len(grp), len(grp))])]) for _ in range(200)])
pd.DataFrame({'tissue': tis, 'quartile': est.index + 1, 'median_expression_log2cpm': P.groupby('jq').ej.median().values, 'triples': P['jq'].value_counts().sort_index().values,
              'coupling': est.values, 'ci_low': bs.quantile(0.025).values, 'ci_high': bs.quantile(0.975).values, 'relative': (est / est[0]).values,
              'relative_ci_low': bs.div(bs[0], axis=0).quantile(0.025).values, 'relative_ci_high': bs.div(bs[0], axis=0).quantile(0.975).values}).to_csv(OUTDIR + f'active_gene_barrier_means_{tis}.csv', index=False)
