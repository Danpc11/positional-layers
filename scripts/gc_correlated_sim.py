"""GC correction when part of the real cis signal tracks GC (Supplementary Table 9C).
A per-sample regulatory factor acts on each gene in proportion to its locally smoothed GC content, making up 0-75% of the cis variance;
we report the fraction of the TRUE cis lag correlation that survives GC correction."""
import numpy as np, pandas as pd
from poslayers import simulate_genome
from poslayers.decompose import gc_correct, lag_profile_linear
from poslayers.config import OUTDIR
rows = []
for share in (0.0, 0.25, 0.5, 0.75):
    out = []
    for rep in range(3):
        s = simulate_genome(n_genes=1000, n_chrom=4, n_samples=200, decay_len=5, gc_bias_sd=0.25, cis_sd=0.6, seed=50 + rep)
        rng = np.random.default_rng(90 + rep); gc = s['gc']; cis = s['cis']
        bio = np.outer(rng.normal(0, 1, cis.shape[0]), np.convolve(gc, np.ones(9) / 9, 'same')); bio *= cis.std() / bio.std()
        new = np.sqrt(1 - share) * cis + np.sqrt(share) * bio
        truth = lag_profile_linear(new - new.mean(0), 1000, 4, 10); kept = lag_profile_linear(gc_correct(new - new.mean(0), gc), 1000, 4, 10)
        out.append((truth[1], truth[10], kept[1] / truth[1], kept[10] / truth[10]))
    rows.append((share, *np.mean(out, 0)))
R = pd.DataFrame(rows, columns=['gc_tracking_share', 'true_lag1', 'true_lag10', 'surviving_lag1', 'surviving_lag10'])
R.to_csv(OUTDIR + 'gc_correlated_cis_sim.csv', index=False); print(R.round(3).to_string(index=False))
