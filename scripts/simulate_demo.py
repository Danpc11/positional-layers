"""The web simulator, in Python. Reproduces Fig. 1b and Fig. 1h.

    python scripts/simulate_demo.py --decay 5 --gc-bias 0.25
"""
import argparse
import numpy as np
from poslayers import simulate_genome, periodogram_identity
from poslayers.decompose import gc_correct, gc_slopes
from poslayers.laws import gc_layer_correlation

# The cis kernel is exp(-i/lambda), whose normalised autocorrelation is exp(-L/lambda).
# The decay length of that autocorrelation is lambda itself, so no calibration is needed.
CAL = 1.0


def lag_profile(Y, n_genes, n_chrom, max_lag=60):
    sd = Y.std(0) + 1e-12
    Z = Y / sd
    out = np.zeros(max_lag + 1)
    for L in range(max_lag + 1):
        acc, cnt = 0.0, 0
        for c in range(n_chrom):
            s = slice(c * n_genes, (c + 1) * n_genes)
            A, B = Z[:, s][:, :n_genes - L], Z[:, s][:, L:]
            acc += float(np.mean(A * B) * A.shape[1]); cnt += A.shape[1]
        out[L] = acc / cnt
    return out


def fit_decay(prof):
    lags = np.arange(1, len(prof))
    v = prof[lags]
    ok = (v > 0.05 * v[0]) & (v > 1e-4)   # relative window: a fixed window lets the noisy tail dominate short decays
    if ok.sum() < 4:
        return np.nan                      # a decay of ~1 gene cannot be fitted this way
    slope = np.polyfit(lags[ok], np.log(v[ok]), 1)[0]
    return -1 / slope / CAL if slope < 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--decay", type=float, default=5.0)
    ap.add_argument("--gc-bias", type=float, default=0.25)
    ap.add_argument("--phi", type=float, default=0.97)
    ap.add_argument("--cis", type=float, default=0.6)
    ap.add_argument("--samples", type=int, default=150)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    s = simulate_genome(n_genes=400, n_chrom=2, n_samples=a.samples, decay_len=a.decay,
                        gc_phi=a.phi, gc_bias_sd=a.gc_bias, cis_sd=a.cis, seed=a.seed)
    X, gc = s["X"], s["gc"]
    Y = X - X.mean(0)

    print(f"identity error      {periodogram_identity(X)['max_relative_error']:.2e}")
    naive = lag_profile(Y, s["n_genes_per_chrom"], s["n_chrom"])
    corr = lag_profile(gc_correct(Y, gc), s["n_genes_per_chrom"], s["n_chrom"])
    print(f"true decay length   {a.decay:.1f} genes")
    print(f"naive estimate      {fit_decay(naive):.1f} genes")
    print(f"GC-corrected        {fit_decay(corr):.1f} genes")
    print(f"adjacent, naive     {naive[1]:.3f}")
    print(f"adjacent, corrected {corr[1]:.3f}")

    var_b = float(np.var(gc_slopes(X, gc)))
    sd = Y.std(0)
    print("\nisochore law, predicted vs observed GC component")
    for L in (1, 2, 5, 10, 20, 30):
        i = np.arange(len(gc) - L)
        pred = float(np.mean(gc_layer_correlation(var_b, gc[i], gc[i + L], sd[i], sd[i + L])))
        print(f"  {L:>2} genes   predicted {pred:+.4f}   observed {naive[L] - corr[L]:+.4f}")


if __name__ == "__main__":
    main()
