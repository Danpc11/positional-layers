"""Simulations for the positional-layer theory.
x[s,j] = mu[j] + beta[j]*z[s] + c[s,j] + b[s]*gc[j] + e[s,j]
  mu: tissue landscape (no positional structure); c: cis co-regulation (moving average of local regulators, decay length lam);
  b*gc: per-sample technical GC bias acting on isochore-clustered GC; beta*z: stage effect, locally smoothed.
Checks: (1) spectral decomposition identity; (2) 'universal peaks' with real vs random gene order; (3) GC creates a false
domain scale; (4) recovery of the cis decay length by naive vs corrected estimators."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import numpy as np, pandas as pd
rng = np.random.default_rng(42)
def ar1(n, phi, rng): x = np.zeros(n); x[0] = rng.normal()
for _ in range(0): pass
def ar1(n, phi, rng):
    x = np.empty(n); x[0] = rng.normal()
    for i in range(1, n): x[i] = phi * x[i - 1] + np.sqrt(1 - phi ** 2) * rng.normal()
    return x
def simulate(n_chr=4, n=1000, S=200, lam=3.0, cis_sd=0.6, gc_sd=0.35, stage_sd=0.15, land_sd=2.0, seed=0, cis=True, gc=True, stage=True):
    r = np.random.default_rng(seed); chrs = np.repeat(np.arange(n_chr), n); N = n_chr * n
    mu = r.normal(0, land_sd, N)                                    # landscape: independent across positions
    GC = np.concatenate([ar1(n, 0.97, r) for _ in range(n_chr)])     # isochore-clustered GC
    z = r.integers(0, 6, S).astype(float)
    X = np.tile(mu, (S, 1)) + r.normal(0, 1, (S, N))
    if cis:
        k = np.arange(-int(6 * lam), int(6 * lam) + 1); w = np.exp(-np.abs(k) / lam); w /= np.sqrt(np.sum(w ** 2))
        for c in range(n_chr):
            sl = slice(c * n, (c + 1) * n); U = r.normal(0, 1, (S, n + 2 * len(k)))
            X[:, sl] += cis_sd * np.array([np.convolve(u, w, mode='same')[len(k):len(k) + n] for u in U])
    if gc: X += np.outer(r.normal(0, gc_sd, S), GC)
    if stage:
        beta = np.concatenate([np.convolve(r.normal(0, 1, n + 20), np.ones(3) / np.sqrt(3), mode='same')[10:10 + n] for _ in range(n_chr)]) * stage_sd
        X += np.outer(z, beta)
    return X, chrs, GC, mu, z
def pgram(x): n = len(x); return np.abs(np.fft.rfft(x)[1:n // 2 + 1]) ** 2
def whiten(P, n):
    k = np.arange(1, len(P) + 1); lx = np.log10(k / n); b = np.polyfit(lx, np.log10(P + 1e-12), 1); w = P / 10 ** (b[0] * lx + b[1]); return w / w.mean()
def universal(X, chrs, perm_seed=None):
    tot = 0
    for c in np.unique(chrs):
        idx = np.where(chrs == c)[0]
        if perm_seed is not None: idx = np.random.default_rng(perm_seed).permutation(idx)
        Xc = X[:, idx] - X[:, idx].mean(1, keepdims=True); W = np.array([whiten(pgram(x), len(idx)) for x in Xc]); tot += int((np.mean(W > 3, 0) >= 0.9).sum())
    return tot
LAGS = np.array([1, 2, 3, 5, 7, 10, 15, 20, 30, 50])
def lagprof(D, chrs):
    out = []
    for L in LAGS:
        v = []
        for c in np.unique(chrs):
            Z = D[:, chrs == c]; Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9); v.append(np.mean(Z[:, :-L] * Z[:, L:]))
        out.append(np.mean(v))
    return np.array(out)
def deviations(X, GC=None, z=None):
    D = X - X.mean(0)                                                 # remove the landscape (per-gene mean across samples)
    if z is not None: Zd = np.column_stack([np.ones(len(z)), z]); D = D - Zd @ np.linalg.lstsq(Zd, D, rcond=None)[0]
    if GC is not None: G = np.column_stack([np.ones(len(GC)), GC, GC ** 2]); D = D - (G @ np.linalg.lstsq(G, D.T, rcond=None)[0]).T
    return D
def decay_length(prof):
    m = prof > 0.01; b = np.polyfit(LAGS[m], np.log(prof[m]), 1) if m.sum() >= 3 else [np.nan]; return -1 / b[0]
if __name__ == '__main__':
    out = {}
    # (1) identity: mean periodogram = |FT(mu_hat)|^2 + FT of the lag-summed (circular) sample covariance
    X, chrs, GC, mu, z = simulate(seed=1); c0 = chrs == 0; Xc = X[:, c0]; n = Xc.shape[1]
    mbar = Xc.mean(0); Dv = Xc - mbar
    lhs = np.mean([np.abs(np.fft.fft(x)) ** 2 for x in Xc], 0)
    C = np.array([np.mean(np.sum(Dv * np.roll(Dv, -L, axis=1), 1)) for L in range(n)])
    rhs = np.abs(np.fft.fft(mbar)) ** 2 + np.real(np.fft.fft(C))
    out['identity_max_rel_error'] = float(np.max(np.abs(lhs - rhs) / lhs))
    out['landscape_share_of_power'] = float(np.sum(np.abs(np.fft.fft(mbar)) ** 2) / np.sum(lhs))
    # (2) universal peaks with real vs random order, with and without positional structure
    for name, kw in [('landscape only', dict(cis=False, gc=False, stage=False)), ('full model', {})]:
        Xs, ch, *_ = simulate(seed=2, **kw); out[f'universal_{name}_real'] = universal(Xs, ch); out[f'universal_{name}_random'] = universal(Xs, ch, perm_seed=5)
    # (3) GC creates a false domain scale; correction recovers the cis-only truth
    Xg, ch, GCg, _, zg = simulate(seed=3, cis=True, gc=True); Xt, cht, *_ = simulate(seed=3, cis=True, gc=False)
    prof_raw = lagprof(deviations(Xg, z=zg), ch); prof_corr = lagprof(deviations(Xg, GC=GCg, z=zg), ch); prof_truth = lagprof(deviations(Xt, z=zg), cht)
    pd.DataFrame({'lag': LAGS, 'raw': prof_raw, 'GC_corrected': prof_corr, 'truth_no_GC': prof_truth}).to_csv(OUTDIR + 'sim_lag_profiles.csv', index=False)
    # (4) recovery of the true cis decay length: see decay_recovery.py (Fig. 1h). The fixed-window estimator that
    #     used to live here was replaced after review; it let the noisy tail dominate when the decay is short.
    pd.Series(out).to_csv(OUTDIR + 'sim_summary.csv')
    print(pd.Series(out).round(4).to_string()); print(pd.DataFrame({'lag': LAGS, 'raw': prof_raw, 'GC_corr': prof_corr, 'truth': prof_truth}).round(3).to_string(index=False))
