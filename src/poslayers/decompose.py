"""Exact decomposition of a positional expression matrix (Supplementary Note 1, sections 1-2)."""
from __future__ import annotations
import numpy as np


def lag_covariance(Y: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Circular lag covariance C(L) of deviations Y (samples x genes), averaged over samples.

    C(L) = (1/S) sum_s sum_j y(s, j) y(s, j + L), the quantity whose Fourier transform
    is the covariance spectrum S(f).

    The sum is CIRCULAR: index j + L wraps around the end of the array, which is what
    makes the Wiener-Khinchin identity in periodogram_identity() exact. Pass one
    chromosome at a time. For a linear (non-wrapping) lag profile at a few short lags,
    use lag_profile_linear() instead; the two agree to O(L / n_genes).
    """
    Y = np.asarray(Y, float)
    n = Y.shape[1]
    max_lag = n if max_lag is None else min(max_lag, n)
    F = np.fft.rfft(Y, axis=1)
    C = np.fft.irfft((F * np.conj(F)).real, n=n, axis=1).mean(0)
    return C[:max_lag]


def periodogram_identity(X: np.ndarray) -> dict:
    """Verify  mean_s |FT(x_s)|^2  ==  |FT(xbar)|^2 + FT(C).

    X is samples x genes, genes ordered along a chromosome. Returns the three spectra and
    the maximum relative error of the identity, which is an algebraic identity and so
    should be at machine precision.
    """
    X = np.asarray(X, float)
    xbar = X.mean(0)
    Y = X - xbar
    lhs = (np.abs(np.fft.rfft(X, axis=1)) ** 2).mean(0)
    land = np.abs(np.fft.rfft(xbar)) ** 2
    cov = (np.abs(np.fft.rfft(Y, axis=1)) ** 2).mean(0)
    err = np.max(np.abs(lhs - land - cov) / (np.abs(lhs) + 1e-12))
    return {"mean_periodogram": lhs, "landscape": land, "covariance": cov, "max_relative_error": float(err)}


def landscape_share(X: np.ndarray) -> float:
    """Fraction of the spectral power of a single sample carried by the tissue landscape."""
    X = np.asarray(X, float)
    xbar = X.mean(0)
    land = (xbar ** 2).sum() * X.shape[0]
    resid = ((X - xbar) ** 2).sum()
    return float(land / (land + resid))



def _gc_basis(gc, degree):
    """Orthonormal basis of [1, g, ..., g^degree] with its EFFECTIVE rank.

    With constant GC (or fewer distinct GC values than the degree), the polynomial columns are collinear.
    A plain QR still returns `degree + 1` columns, and projecting on the spurious ones removes signal unrelated to GC.
    We keep only directions whose singular value is non-negligible."""
    g = np.asarray(gc, float)
    sd = g.std()
    g = (g - g.mean()) / sd if sd > 0 else np.zeros_like(g)
    D = np.vstack([g ** k for k in range(degree + 1)]).T
    U, s, _ = np.linalg.svd(D, full_matrices=False)
    keep = s > s[0] * 1e-10 * max(D.shape)
    return U[:, keep]


def gc_correct(X: np.ndarray, gc: np.ndarray, degree: int = 2) -> np.ndarray:
    """Remove a per-sample polynomial trend on gene GC content (the technical isochore layer).

    X is samples x genes. The basis includes the constant column, so each sample's mean across genes is also removed: a
    per-sample offset is not positional information, and rows come out centred. When GC is constant the basis has rank 1
    and only the per-sample mean is removed. On simulated data with no GC bias the correction costs 0.3% of the
    adjacent-gene correlation; when real regulation tracks GC it can remove much more at domain scale (see
    scripts/gc_correlated_sim.py), so it cannot separate technical bias from GC-associated biology.
    """
    X = np.asarray(X, float)
    if X.ndim != 2 or X.shape[1] != len(gc):
        raise ValueError(f'X must be samples x genes with {len(gc)} genes; got shape {X.shape}')
    Q = _gc_basis(gc, degree)
    return X - (X @ Q) @ Q.T


def gc_slopes(X: np.ndarray, gc: np.ndarray) -> np.ndarray:
    """Per-sample GC slope b_s (OLS on standardised GC with an intercept), whose variance drives the GC-associated layer.
    Returns zeros when GC is constant, where the slope is not identifiable."""
    X = np.asarray(X, float)
    g = np.asarray(gc, float)
    if X.shape[1] != len(g):
        raise ValueError(f'X must be samples x genes with {len(g)} genes; got shape {X.shape}')
    if g.std() == 0:
        return np.zeros(X.shape[0])
    g = (g - g.mean()) / g.std()
    Y = X - X.mean(1, keepdims=True)
    return (Y @ g) / (g @ g)


def lag_profile_linear(Y: np.ndarray, n_genes: int, n_chrom: int, max_lag: int = 60) -> np.ndarray:
    """Non-circular lag correlation profile, averaged over chromosomes and samples.

    Y is samples x genes with the genes of n_chrom chromosomes of n_genes each, concatenated. This is the estimator used for
    the coupling atlas: genes standardised across samples, pairs taken within a chromosome without wrapping. Lags with no
    pairs (max_lag >= n_genes) are returned as NaN instead of failing.
    """
    Y = np.asarray(Y, float)
    if Y.ndim != 2 or Y.shape[1] != n_genes * n_chrom:
        raise ValueError(f'Y must be samples x (n_genes * n_chrom) = {n_genes * n_chrom} columns; got shape {Y.shape}')
    Z = (Y - Y.mean(0)) / (Y.std(0) + 1e-12)
    out = np.full(max_lag + 1, np.nan)
    for L in range(min(max_lag, n_genes - 1) + 1):
        acc, cnt = 0.0, 0
        for c in range(n_chrom):
            s = slice(c * n_genes, (c + 1) * n_genes)
            A, B = Z[:, s][:, :n_genes - L], Z[:, s][:, L:]
            acc += float(np.mean(A * B)) * A.shape[1]; cnt += A.shape[1]
        out[L] = acc / cnt
    return out
