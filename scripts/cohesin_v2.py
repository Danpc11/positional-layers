"""Per-tumour cis-excess score and covariate-adjusted test of cohesin/CTCF loss.
Score: for each tumour, mean product of pooled-standardised residual expression of genes 1-3 positions apart minus that of
genes 20-30 apart (GC- and own-copy-number-corrected expression). Covariates: expression-based immune and stromal scores
(purity proxy), copy-number burden, log mutation burden, cohort-specific subtype score."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm
from lib_tumour import coding, prepare, scores, trunc
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
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'cohesin_v2_results.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
