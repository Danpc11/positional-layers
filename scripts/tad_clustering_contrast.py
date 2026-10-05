"""Does excluding copy-number-altered genes reduce the tumour-normal difference in TAD clustering of active genes?
(Extended Data Fig. 2c; Supplementary Table 6G)

Tests the change of the contrast directly instead of comparing significant / non-significant labels:
  delta_before = median z(tumour) - median z(normal)
  delta_after  = median z(tumour, CNA genes removed) - median z(normal, matched gene exclusion)
  reduction    = delta_before - delta_after
95% intervals from 5,000 bootstrap replicates that resample tumours and normals (each sample keeps both of its values, so
the before/after pairing within a sample is preserved); two-sided bootstrap P for reduction = 0.
Input : OUTDIR/friction2_per_sample_z.csv (fric2.py)   Output: OUTDIR/tad_clustering_contrast.csv
"""
import numpy as np, pandas as pd
from poslayers.config import OUTDIR

Z = pd.read_csv(OUTDIR + 'friction2_per_sample_z.csv'); rows = []
for c, d in Z.groupby('cohort'):
    W = d.pivot_table(index='sample', columns='group', values='z')
    T = W[['tumour', 'tumour_CNA_removed']].dropna(); N = W[['normal', 'normal_matched_exclusion']].dropna()
    def stat(t, n):
        b = np.median(t[:, 0]) - np.median(n[:, 0]); a = np.median(t[:, 1]) - np.median(n[:, 1]); return b, a, b - a
    est = stat(T.values, N.values); rng = np.random.default_rng(4)
    bs = np.array([stat(T.values[rng.integers(0, len(T), len(T))], N.values[rng.integers(0, len(N), len(N))]) for _ in range(5000)])
    p = 2 * min((bs[:, 2] <= 0).mean(), (bs[:, 2] >= 0).mean())
    rows.append({'cohort': c, 'tumours': len(T), 'normals': len(N), 'delta_before': est[0], 'delta_before_ci_low': np.percentile(bs[:, 0], 2.5), 'delta_before_ci_high': np.percentile(bs[:, 0], 97.5),
                 'delta_after': est[1], 'delta_after_ci_low': np.percentile(bs[:, 1], 2.5), 'delta_after_ci_high': np.percentile(bs[:, 1], 97.5),
                 'reduction': est[2], 'reduction_ci_low': np.percentile(bs[:, 2], 2.5), 'reduction_ci_high': np.percentile(bs[:, 2], 97.5), 'p_reduction': max(p, 1 / 5000)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'tad_clustering_contrast.csv', index=False); print(R.round(3).to_string(index=False))
