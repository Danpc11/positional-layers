"""Out-of-sample prediction of the copy-number layer from the shared-source model (Fig. 3b; Supplementary Table 9A).

Theory: the dosage component of expression covariance between genes i and j is beta_i * beta_j * Cov(CN_i, CN_j), where
beta is each gene's response to its own copy number (its loading on the copy-number source) and Cov(CN_i, CN_j) comes
from the segment structure alone. Out-of-sample test: beta is estimated in a random half of the tumours; in the other
half the prediction is built from copy-number data only and compared with the observed dosage component
(covariance before minus after each gene's own copy-number correction), averaged over gene pairs at each lag.
Stability: the prediction is repeated for 10 random splits of tumours (the first split, seed 5, is the one shown in Fig. 3b).
The copy-number correction that defines the observed component (linear and quadratic in each gene's own copy number,
lib_tumour.prepare) is fitted on all tumours; the footprints beta are linear slopes fitted on the training half only.
Outputs: OUTDIR/dosage_prediction.csv (seed 5), OUTDIR/dosage_prediction_splits.csv (ratio by cohort, split and lag)
"""
import sys, numpy as np, pandas as pd
from poslayers.config import OUTDIR
from lib_tumour import prepare

LAGS = [1, 2, 3, 5, 8, 12, 20, 30, 50, 80]
rows = []
for c in (sys.argv[1].split(',') if len(sys.argv) > 1 else ['BLCA', 'UCEC']):
    D, D_raw, chrs, samp, cov, gids = prepare(c, cn='continuous', return_raw=True, return_gids=True)
    Z = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CN = pd.DataFrame(Z['CN'], index=Z['genes'], columns=Z['samples']).reindex(index=gids, columns=samp).fillna(0).values
    for split, seed in enumerate([5] + list(range(101, 110))):
        rng = np.random.default_rng(seed); half = rng.permutation(len(samp)); A, B = half[: len(samp) // 2], half[len(samp) // 2:]
        xa = CN[:, A] - CN[:, A].mean(1, keepdims=True); ya = D_raw[:, A] - D_raw[:, A].mean(1, keepdims=True)
        beta = np.sum(xa * ya, 1) / np.maximum(np.sum(xa * xa, 1), 1e-9)                 # footprints from half A
        cB = lambda M: M[:, B] - M[:, B].mean(1, keepdims=True)
        Rb, Cb, Kb = cB(D_raw), cB(D), cB(CN)
        for L in LAGS:
            o, p, n = 0.0, 0.0, 0
            for ch in np.unique(chrs):
                x = np.where(chrs == ch)[0]
                if len(x) <= L + 5: continue
                i, j = x[:-L], x[L:]
                o += np.sum(np.mean(Rb[i] * Rb[j], 1) - np.mean(Cb[i] * Cb[j], 1))
                p += np.sum(beta[i] * beta[j] * np.mean(Kb[i] * Kb[j], 1)); n += len(i)
            rows.append({'cohort': c, 'split': split, 'seed': seed, 'lag': L, 'observed_dosage_cov': o / n, 'predicted_dosage_cov': p / n})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); R['ratio'] = R.observed_dosage_cov / R.predicted_dosage_cov
R.to_csv(OUTDIR + 'dosage_prediction_splits.csv', index=False)
P = R[R.split == 0].drop(columns=['split', 'seed']); P.to_csv(OUTDIR + 'dosage_prediction.csv', index=False); print(P.round(4).to_string(index=False))
for c, d in R.groupby('cohort'):
    print(c, 'primary split: ratio %.3f-%.3f, r across lags %.3f; 10 splits: ratio %.3f-%.3f' % (d[d.split == 0].ratio.min(), d[d.split == 0].ratio.max(),
          np.corrcoef(d[d.split == 0].observed_dosage_cov, d[d.split == 0].predicted_dosage_cov)[0, 1], d.ratio.min(), d.ratio.max()))
