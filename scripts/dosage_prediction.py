"""Out-of-sample prediction of the copy-number layer from the shared-source model (Fig. 3b; Supplementary Table 9A,F).

For each tumour cohort and each of ten random splits into a training and a validation half:
  1. In the TRAINING half only, each gene's GC-corrected expression is regressed on its own continuous copy number with
     an intercept, a linear and a quadratic term: f_j(CN) = b0_j + b1_j CN + b2_j CN^2 (the dosage model).
  2. In the VALIDATION half, the trained model is applied without refitting:
       predicted dosage covariance  = Cov(f_i(CN_i), f_j(CN_j))   -- uses validation copy number and training coefficients
                                                                     only, no validation expression;
       observed dosage covariance   = Cov(x_i, x_j) - Cov(e_i, e_j), with e = x - f(CN)  -- an operational estimand that
                                                                     also contains cross terms between dosage and residual.
     Both are averaged over gene pairs at ten lags of 1-80 genes.
  3. Sensitivity: the linear-only prediction beta_i beta_j Cov(CN_i, CN_j), with beta the training slope.
Nothing is fitted or recalibrated in the validation half. The first split (seed 5) is the one drawn in Fig. 3b.
Usage: python dosage_prediction.py BLCA,UCEC
Outputs: OUTDIR/dosage_prediction_splits.csv (all splits; rows replaced by cohort, split and lag on re-runs) and
         OUTDIR/dosage_prediction.csv (split 0)
"""
import sys, numpy as np, pandas as pd
from poslayers.config import OUTDIR, upsert_csv
from lib_tumour import prepare

LAGS = [1, 2, 3, 5, 8, 12, 20, 30, 50, 80]
SEEDS = [5] + list(range(101, 110))

def fit_dosage(Y, C):
    """Per-gene least squares of Y on [1, C, C^2] (genes x samples); returns coefficients (genes x 3)."""
    X = np.stack([np.ones_like(C), C, C ** 2], -1)                       # genes x samples x 3
    XtX = np.einsum('gsi,gsj->gij', X, X) + 1e-9 * np.eye(3); XtY = np.einsum('gsi,gs->gi', X, Y)
    return np.linalg.solve(XtX, XtY[..., None])[..., 0]
def lagged_cov(M, chrs, L):
    M = M - M.mean(1, keepdims=True); tot, n = 0.0, 0
    for ch in np.unique(chrs):
        x = np.where(chrs == ch)[0]
        if len(x) <= L + 5: continue
        tot += np.sum(np.mean(M[x[:-L]] * M[x[L:]], 1)); n += len(x) - L
    return tot / n

rows = []
for c in sys.argv[1].split(','):
    D, D_raw, chrs, samp, cov, gids = prepare(c, cn='continuous', return_raw=True, return_gids=True)
    Z = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CN = pd.DataFrame(Z['CN'], index=Z['genes'], columns=Z['samples']).reindex(index=gids, columns=samp).fillna(0).values
    for split, seed in enumerate(SEEDS):
        half = np.random.default_rng(seed).permutation(len(samp)); A, B = half[: len(samp) // 2], half[len(samp) // 2:]
        coef = fit_dosage(D_raw[:, A], CN[:, A])                                          # training half only
        CB = CN[:, B]; F = coef[:, [0]] + coef[:, [1]] * CB + coef[:, [2]] * CB ** 2       # trained model applied to validation
        XB = D_raw[:, B]; E = XB - F
        xa = CN[:, A] - CN[:, A].mean(1, keepdims=True); ya = D_raw[:, A] - D_raw[:, A].mean(1, keepdims=True)
        beta = np.sum(xa * ya, 1) / np.maximum(np.sum(xa * xa, 1), 1e-9)
        for L in LAGS:
            obs = lagged_cov(XB, chrs, L) - lagged_cov(E, chrs, L)
            pred = lagged_cov(F, chrs, L)
            Kc = CB - CB.mean(1, keepdims=True); lin, n = 0.0, 0
            for ch in np.unique(chrs):
                x = np.where(chrs == ch)[0]
                if len(x) <= L + 5: continue
                lin += np.sum(beta[x[:-L]] * beta[x[L:]] * np.mean(Kc[x[:-L]] * Kc[x[L:]], 1)); n += len(x) - L
            rows.append({'cohort': c, 'split': split, 'seed': seed, 'lag': L, 'observed_dosage_cov': obs, 'predicted_dosage_cov': pred, 'predicted_linear_only': lin / n})
    print(c, 'done', flush=True)
R = pd.DataFrame(rows); R['ratio'] = R.observed_dosage_cov / R.predicted_dosage_cov; R['ratio_linear_only'] = R.observed_dosage_cov / R.predicted_linear_only
R = upsert_csv(R, OUTDIR + 'dosage_prediction_splits.csv', ['cohort', 'split', 'lag'])
P = R[R.split == 0].drop(columns=['split', 'seed']); P.to_csv(OUTDIR + 'dosage_prediction.csv', index=False)
for c, d in R.groupby('cohort'):
    p0 = d[d.split == 0]
    print(f'{c}: split 0 ratio {p0.ratio.min():.3f}-{p0.ratio.max():.3f} (r across lags {np.corrcoef(p0.observed_dosage_cov, p0.predicted_dosage_cov)[0, 1]:.3f}); '
          f'10 splits {d.ratio.min():.3f}-{d.ratio.max():.3f}; linear-only sensitivity {d.ratio_linear_only.min():.3f}-{d.ratio_linear_only.max():.3f}')
