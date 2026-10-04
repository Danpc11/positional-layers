"""Disease bystanders (Fig. 5b): concordance of fibrosis stage effects between adjacent genes versus their baseline coupling in GTEx liver.
Inputs: pairs_liver.csv.gz (atlas.py) and the per-gene stage effects of the liver biopsy cohort (DATA/liver_tables/Table_S6c_...csv).
The 95% interval resamples 10-Mb genomic blocks (2,000 replicates); donor-level data are not available, so it does not reflect donor sampling."""
import numpy as np, pandas as pd, pyannotables as pa, statsmodels.formula.api as smf
from scipy import stats
from poslayers.config import DATA, OUTDIR
P = pd.read_csv(OUTDIR + 'pairs_liver.csv.gz'); P = P[P.dist > 0]
T = pd.read_csv(DATA + 'liver_tables/Table_S6c_stage_effect_per_gene_with_without_composition.csv').set_index('gene_id')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
rows = []
for col in ('t_stage', 't_stage_comp_adjusted'):
    d = P.assign(t1=T[col].reindex(P.g1).values, t2=T[col].reindex(P.g2).values).dropna(subset=['t1', 't2']).reset_index(drop=True)
    d['prod'] = (d.t1 - d.t1.mean()) / d.t1.std() * (d.t2 - d.t2.mean()) / d.t2.std()
    d.to_csv(OUTDIR + f'bystander_liver_{col}.csv', index=False)
    m = smf.ols('prod ~ r + C(orientation) + np.log10(dist)', data=d).fit()
    blk = (G.Chromosome.reindex(d.g1).astype(str).values + ':' + (G.Start.reindex(d.g1).fillna(0).values // 1e7).astype(int).astype(str))
    grp = [np.where(blk == b)[0] for b in np.unique(blk)]; rng = np.random.default_rng(11)
    sl = lambda ix: np.polyfit(d.r.values[ix], d['prod'].values[ix], 1)[0]
    est = sl(np.arange(len(d))); bs = np.array([sl(np.concatenate([grp[k] for k in rng.integers(0, len(grp), len(grp))])) for _ in range(2000)])
    rows.append({'stage_effect': col, 'pairs': len(d), 'blocks': len(grp), 'slope': est, 'slope_adjusted_orientation_distance': m.params['r'],
                 'ci_low': np.percentile(bs, 2.5), 'ci_high': np.percentile(bs, 97.5), 'boot_se': bs.std(), 'p_normal_approx': 2 * stats.norm.sf(abs(est / bs.std()))})
pd.DataFrame(rows).to_csv(OUTDIR + 'bystander_block_bootstrap.csv', index=False); print(pd.DataFrame(rows).round(4).to_string(index=False))
