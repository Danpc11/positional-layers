"""Regression tests for the problems found in the audit of the repository."""
import numpy as np
from poslayers import simulate_genome, lag_covariance, lag_profile_linear
from poslayers.decompose import gc_correct


def _fit_decay(prof):
    """Log-linear fit over the lags where the profile is still above 5% of its lag-1 value.

    A fixed lag window lets the noisy tail dominate when the decay is short: with a fixed
    window of 30 lags, a true length of 2 genes came back as 6. The relative window
    recovers lambda to within 3% for lambda >= 2, and returns nan below that, which is the
    honest answer for a decay of about one gene.
    """
    lags = np.arange(1, len(prof))
    v = prof[lags]
    ok = (v > 0.05 * v[0]) & (v > 1e-4)
    if ok.sum() < 4:
        return np.nan
    slope = np.polyfit(lags[ok], np.log(v[ok]), 1)[0]
    return -1 / slope if slope < 0 else np.nan


def test_kernel_autocorrelation_decays_with_lambda():
    """No calibration constant: the decay length of the autocorrelation IS lambda.

    This is the bug the audit found. The code used to divide by 1.56, which made the
    GC-corrected estimator undershoot by about a third.
    """
    for lam in (2.0, 3.0, 5.0, 8.0, 12.0):
        s = simulate_genome(n_genes=800, n_chrom=4, n_samples=300, decay_len=lam,
                            gc_bias_sd=0.0, noise_sd=0.0, cis_sd=1.0, seed=2)
        Y = s["X"] - s["X"].mean(0)
        fitted = _fit_decay(lag_profile_linear(Y, s["n_genes_per_chrom"], s["n_chrom"], 40))
        assert 0.85 * lam < fitted < 1.15 * lam, (lam, fitted)


def test_gc_correction_does_not_eat_real_signal():
    """With no GC bias at all, correcting costs almost nothing."""
    s = simulate_genome(n_genes=800, n_chrom=4, n_samples=200, decay_len=5.0,
                        gc_bias_sd=0.0, seed=6)
    Y = s["X"] - s["X"].mean(0)
    raw = lag_profile_linear(Y, 800, 4, 5)[1]
    cor = lag_profile_linear(gc_correct(Y, s["gc"]), 800, 4, 5)[1]
    assert cor > 0.95 * raw


def test_gc_correction_removes_a_pure_artefact():
    """With a GC bias and no cis layer, correction should flatten the lag profile."""
    s = simulate_genome(n_genes=800, n_chrom=4, n_samples=200, decay_len=3.0,
                        gc_bias_sd=0.6, cis_sd=0.0, noise_sd=0.3, seed=11)
    Y = s["X"] - s["X"].mean(0)
    raw = lag_profile_linear(Y, 800, 4, 20)
    cor = lag_profile_linear(gc_correct(Y, s["gc"]), 800, 4, 20)
    assert abs(cor[10]) < 0.5 * abs(raw[10])


def test_lag_covariance_is_circular_as_documented():
    """lag_covariance wraps; lag_profile_linear does not. Both are intended."""
    rng = np.random.default_rng(0)
    Y = rng.normal(size=(40, 120))
    n = Y.shape[1]
    circ = [float(np.mean([np.sum(Y[s] * np.roll(Y[s], -L)) for s in range(40)])) for L in range(4)]
    assert np.allclose(lag_covariance(Y, max_lag=4), circ)


def test_gc_correct_centres_each_sample():
    """The polynomial basis includes the constant column, so rows come out centred."""
    rng = np.random.default_rng(1)
    g = rng.normal(size=200)
    X = rng.normal(size=(20, 200)) + 7.0
    assert np.allclose(gc_correct(X, g).mean(1), 0.0, atol=1e-9)
