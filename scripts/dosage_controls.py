"""Controls for the copy-number layer (Fig. 3b): null models and adjustment for tumour purity.

Split 0 of dosage_prediction.py (seed 5), same dosage model fitted in the training half only.
  true      : each gene's own copy number;
  permuted  : copy-number profiles permuted across tumours (the link between a tumour's copy number and its
              expression is broken; segment structure is kept);
  shifted   : each chromosome's copy-number profile shifted by half its length (circularly) in every tumour, so that
              segments are kept but no longer aligned with their genes;
  purity    : true copy number, with expression-based immune and stromal scores (ESTIMATE-like purity proxies)
              added to the training model; their fitted contribution is removed in both halves before the dosage terms.
For each: observed dosage covariance (covariance of expression minus that of its residual) and predicted dosage
covariance (covariance of the fitted dosage components), averaged over the ten lags of dosage_prediction.py.
Usage: python dosage_controls.py BLCA,UCEC    Output: OUTDIR/dosage_controls.csv
"""
import sys
import numpy as np, pandas as pd
from poslayers.config import OUTDIR, replace_groups_csv
from lib_tumour import prepare

LAGS = [1, 2, 3, 5, 8, 12, 20, 30, 50, 80]
def fit(Y, Xd):                                   # per-gene least squares; Xd: genes x samples x k
    XtX = np.einsum('gsi,gsj->gij', Xd, Xd) + 1e-9 * np.eye(Xd.shape[2]); return np.linalg.solve(XtX, np.einsum('gsi,gs->gi', Xd, Y)[..., None])[..., 0]
def lagged_cov(M, chrs, L):
    M = M - M.mean(1, keepdims=True); tot, n = 0.0, 0
    for ch in np.unique(chrs):
        x = np.where(chrs == ch)[0]
        if len(x) <= L + 5: continue
        tot += np.sum(np.mean(M[x[:-L]] * M[x[L:]], 1)); n += len(x) - L
    return tot / n
for c in sys.argv[1].split(','):
    D, D_raw, chrs, samp, cov, gids = prepare(c, cn='continuous', return_raw=True, return_gids=True)
    Z = np.load(OUTDIR + f'cn_cont_{c}.npz', allow_pickle=True); CN = pd.DataFrame(Z['CN'], index=Z['genes'], columns=Z['samples']).reindex(index=gids, columns=samp).fillna(0).values
    half = np.random.default_rng(5).permutation(len(samp)); A, B = half[: len(samp) // 2], half[len(samp) // 2:]
    pur = cov[['immune', 'stromal']].values; pur = (pur - pur.mean(0)) / pur.std(0)
    CNp = CN[:, np.random.default_rng(9).permutation(len(samp))]
    CNs = CN.copy()
    for ch in np.unique(chrs):
        x = np.where(chrs == ch)[0]; CNs[x] = np.roll(CN[x], len(x) // 2, axis=0)
    rows = []
    for label, C, use_pur in (('true', CN, False), ('permuted across tumours', CNp, False), ('segments shifted', CNs, False), ('purity-adjusted', CN, True)):
        G, S = D_raw.shape
        base = [np.ones((G, S)), C, C ** 2] + ([np.broadcast_to(pur[:, k], (G, S)) for k in range(2)] if use_pur else [])
        Xd = np.stack(base, -1); coef = fit(D_raw[:, A], Xd[:, A])
        F = coef[:, [1]] * C[:, B] + coef[:, [2]] * C[:, B] ** 2                         # dosage component only
        P = sum(coef[:, [3 + k]] * pur[B, k][None, :] for k in range(2)) if use_pur else 0
        XB = D_raw[:, B] - P; E = XB - F
        obs = np.mean([lagged_cov(XB, chrs, L) - lagged_cov(E, chrs, L) for L in LAGS]); pred = np.mean([lagged_cov(F, chrs, L) for L in LAGS])
        rows.append({'cohort': c, 'copy_number': label, 'observed_dosage_cov': obs, 'predicted_dosage_cov': pred, 'ratio': obs / pred if pred else np.nan})
    R = pd.DataFrame(rows); R['observed_relative_to_true'] = R.observed_dosage_cov / R.observed_dosage_cov.iloc[0]
    print(R.round(4).to_string(index=False), flush=True)
    replace_groups_csv(R, OUTDIR + 'dosage_controls.csv', ['cohort'], [c], key=['cohort', 'copy_number'])
