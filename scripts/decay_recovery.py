"""Recovery of the simulated cis decay length, naive versus GC-corrected (Extended Data Fig. 1c).

Genomes of 4 x 1,000 genes and 200 samples (poslayers.simulate_genome; AR(1) GC track with phi = 0.97, per-sample GC bias
SD 0.25) for lambda = 2, 3, 5, 8 and 12 genes, three replicates each. The decay length is fitted to the linear lag
profile over the lags at which the profile is still above 5% of its lag-1 value (no calibration factor; see
Supplementary Note 1, section 7).

Output: OUTDIR/sim_decay_recovery.csv (true_lambda, rep, naive, corrected)
"""
import numpy as np, pandas as pd
from poslayers import simulate_genome
from poslayers.decompose import gc_correct, lag_profile_linear
from poslayers.config import OUTDIR


def fit_decay(p):
    lags = np.arange(1, len(p)); v = p[1:]; ok = (v > 0.05 * v[0]) & (v > 1e-4)
    if ok.sum() < 4: return np.nan
    s = np.polyfit(lags[ok], np.log(v[ok]), 1)[0]
    return -1 / s if s < 0 else np.nan


rows = []
for lam in (2, 3, 5, 8, 12):
    for rep in range(3):
        s = simulate_genome(n_genes=1000, n_chrom=4, n_samples=200, decay_len=lam, gc_phi=0.97, gc_bias_sd=0.25, seed=100 * lam + rep)
        Y = s['X'] - s['X'].mean(0)
        rows.append(dict(true_lambda=lam, rep=rep, naive=fit_decay(lag_profile_linear(Y, 1000, 4, 60)),
                         corrected=fit_decay(lag_profile_linear(gc_correct(Y, s['gc']), 1000, 4, 60))))
D = pd.DataFrame(rows); D.to_csv(OUTDIR + 'sim_decay_recovery.csv', index=False)
print(D.groupby('true_lambda')[['naive', 'corrected']].median().round(2).to_string())
