"""The three quantitative laws (Supplementary Note 1, sections 3-5)."""
from __future__ import annotations
import numpy as np
from scipy.optimize import curve_fit


def gc_layer_correlation(var_b: float, gc_i, gc_j, sd_i, sd_j) -> np.ndarray:
    """Correlation between two genes induced by per-sample GC bias: the GC-associated layer (Supplementary Note 1,
    section 2). Called isochore_law in release 0.1.0; that name is kept as an alias.

    rho_GC(i, j) = Var(b) * g_i * g_j / (sigma_i * sigma_j)

    var_b   variance across samples of the per-sample GC slope (see decompose.gc_slopes)
    gc_i/j  standardised GC content of each gene
    sd_i/j  standard deviation across samples of each gene's deviations
    """
    return var_b * np.asarray(gc_i) * np.asarray(gc_j) / (np.asarray(sd_i) * np.asarray(sd_j))


def predicted_genetic_covariance(freq, beta1, beta2, p_same_causal=1.0) -> np.ndarray:
    """Covariance (correlation, for unit-variance expression) between two genes induced by a shared causal variant,
    2p(1 - p) beta1 beta2 (Supplementary Note 1, section 4). Called eqtl_law in release 0.1.0; kept as an alias.

    rho_eQTL = P(same variant) * 2p(1 - p) * beta1 * beta2

    Effects must be on expression scaled to unit variance (GTEx slopes are estimated on
    inverse-normal-transformed expression, so they can be used directly). The sign of the
    product sets the sign of the coupling: this is the law's falsifiable part.
    """
    p = np.asarray(freq, float)
    return np.asarray(p_same_causal, float) * 2 * p * (1 - p) * np.asarray(beta1, float) * np.asarray(beta2, float)


def saturation_exponent(contact, K: float) -> np.ndarray:
    """Local exponent of coupling on contact, k = 1 - occupancy = K / (c + K).

    Not used in the article; kept for compatibility with release 0.1.0."""
    c = np.asarray(contact, float)
    return K / (c + K)


def _power(c, r0, A, k):
    return r0 + A * c ** k


def _hub(c, r0, A, cstar):
    return r0 + A * (1 + (cstar / c) ** (2 / 3)) ** -1.5


def fit_contact_law(contact, coupling, sigma=None) -> dict:
    """Fit the power law and the fixed-size hub model to mean coupling versus mean contact.
    The article reports only the power-law exponent, as a description of the Hi-C maps; the hub model is kept for
    compatibility with release 0.1.0.

    Returns both fits with their weighted chi-square, so the hub model can be rejected as it
    is in the paper. The power-law exponent k is the quantity that equals one minus the
    mean regulatory occupancy of the pairs in the bin.
    """
    c = np.asarray(contact, float)
    r = np.asarray(coupling, float)
    kw = dict(sigma=sigma, maxfev=50000)
    p_pw, _ = curve_fit(_power, c, r, p0=[0.0, 0.002, 0.5], bounds=([-0.05, 0, 0.01], [0.05, 5, 3]), **kw)
    p_hb, _ = curve_fit(_hub, c, r, p0=[0.005, 0.2, float(np.median(c))], bounds=([-0.05, 0, 1], [0.05, 5, 1e6]), **kw)
    w = np.ones_like(r) if sigma is None else 1 / np.asarray(sigma, float) ** 2
    return {
        "power_law": {"r0": p_pw[0], "A": p_pw[1], "k": p_pw[2], "chi2": float(np.sum(w * (r - _power(c, *p_pw)) ** 2))},
        "hub": {"r0": p_hb[0], "A": p_hb[1], "c_star": p_hb[2], "chi2": float(np.sum(w * (r - _hub(c, *p_hb)) ** 2))},
    }


def source_covariance(source_variance, footprint_i, footprint_j):
    """Covariance of two genes under the shared-source model (Supplementary Note 1, section 3).

    Each source q (regulatory element, genotype, copy-number segment) has variance Var(A_q) across samples and reaches
    gene i with footprint a_iq. If sources are uncorrelated with each other and with gene-private noise,
        Cov(x_i, x_j) = sum_q Var(A_q) * a_iq * a_jq.
    The GC layer (footprint = gene GC, source = per-sample GC slope) and the copy-number layer (footprint = dosage
    response, source = copy number) are special cases.
    """
    import numpy as np
    v, a, b = (np.asarray(x, dtype=float) for x in (source_variance, footprint_i, footprint_j))
    if not (v.shape == a.shape == b.shape):
        raise ValueError('source_variance and the two footprints must have the same shape')
    return float(np.sum(v * a * b))


# Names used in release 0.1.0, kept so that existing code keeps working.
isochore_law = gc_layer_correlation
eqtl_law = predicted_genetic_covariance
