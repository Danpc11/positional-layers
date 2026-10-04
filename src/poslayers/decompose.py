"""Exact decomposition of a positional expression matrix (Supplementary Note 1, sections 1-2)."""
from __future__ import annotations
import numpy as np


def lag_covariance(Y: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Circular lag covariance C(L) of deviations Y (samples x genes), averaged over samples.

    C(L) = (1/S) sum_s sum_j y(s, j) y(s, j + L), the quantity whose Fourier transform
    is the covariance spectrum S(f).
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
    X = np.asarray(float_or(X), float)
    xbar = X.mean(0)
    land = (xbar ** 2).sum() * X.shape[0]
    resid = ((X - xbar) ** 2).sum()
    return float(land / (land + resid))


def float_or(x):
    return x


def gc_correct(X: np.ndarray, gc: np.ndarray, degree: int = 2) -> np.ndarray:
    """Remove a per-sample polynomial trend on gene GC content (the technical isochore layer)."""
    X = np.asarray(X, float)
    g = (np.asarray(gc, float) - np.mean(gc)) / (np.std(gc) + 1e-12)
    D = np.vstack([g ** k for k in range(degree + 1)]).T
    Q, _ = np.linalg.qr(D)
    return X - (X @ Q) @ Q.T


def gc_slopes(X: np.ndarray, gc: np.ndarray) -> np.ndarray:
    """Per-sample GC slope b_s, whose variance drives the isochore law."""
    X = np.asarray(X, float)
    g = (np.asarray(gc, float) - np.mean(gc)) / (np.std(gc) + 1e-12)
    Y = X - X.mean(0)
    return (Y @ g) / (g @ g)
