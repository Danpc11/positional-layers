"""poslayers — decomposition of positional signals in transcriptomic data.

Reference implementation of the three laws described in
"Common cis-regulatory inputs shape local gene co-expression across human tissues".
"""
from .decompose import (periodogram_identity, landscape_share, lag_covariance,
                        lag_profile_linear, gc_correct, gc_slopes)
from .laws import isochore_law, eqtl_law, saturation_exponent, fit_contact_law, source_covariance
from .simulate import simulate_genome

__version__ = "1.3.1"
__all__ = ["source_covariance", "periodogram_identity", "landscape_share", "lag_covariance",
           "lag_profile_linear", "gc_correct", "gc_slopes",
           "lag_profile_linear", "gc_correct", "gc_slopes",
           "isochore_law", "eqtl_law", "saturation_exponent", "fit_contact_law",
           "simulate_genome"]
