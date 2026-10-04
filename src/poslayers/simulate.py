"""Simulated genomes with a known cis structure, used to validate the estimators."""
from __future__ import annotations
import numpy as np


def simulate_genome(n_genes=1000, n_chrom=4, n_samples=200, decay_len=5.0,
                    gc_phi=0.97, gc_bias_sd=0.25, cis_sd=0.6, noise_sd=0.4, seed=0) -> dict:
    """Generate an expression matrix with a landscape, a GC-bias layer and a known cis layer.

    decay_len is the true exponential decay length of the cis kernel, in genes. The naive
    estimator of that length is biased upward by the GC layer; the GC-corrected estimator
    recovers it. Returns samples x genes matrices and the ground truth.
    """
    rng = np.random.default_rng(seed)
    n = n_genes * n_chrom
    landscape = rng.gamma(1.3, 1.0, n) * 3

    gc = np.zeros(n)
    for c in range(n_chrom):
        sl = slice(c * n_genes, (c + 1) * n_genes)
        e = rng.normal(0, 1, n_genes)
        x = np.zeros(n_genes)
        for i in range(1, n_genes):
            x[i] = gc_phi * x[i - 1] + e[i]
        gc[sl] = (x - x.mean()) / (x.std() + 1e-12)

    L = max(1, int(np.ceil(6 * decay_len)))
    kern = np.exp(-np.arange(L) / decay_len)
    kern /= np.sqrt((kern ** 2).sum())
    cis = np.zeros((n_samples, n))
    for c in range(n_chrom):
        sl = slice(c * n_genes, (c + 1) * n_genes)
        drive = rng.normal(0, 1, (n_samples, n_genes + L))
        cis[:, sl] = np.apply_along_axis(lambda v: np.convolve(v, kern, "valid")[:n_genes], 1, drive) * cis_sd

    b = rng.normal(0, gc_bias_sd, n_samples)
    X = landscape[None, :] + b[:, None] * gc[None, :] + cis + rng.normal(0, noise_sd, (n_samples, n))
    return {"X": X, "gc": gc, "landscape": landscape, "cis": cis, "gc_slope": b,
            "true_decay_len": decay_len, "n_genes_per_chrom": n_genes, "n_chrom": n_chrom}
