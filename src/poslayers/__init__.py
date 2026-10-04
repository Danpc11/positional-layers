"""poslayers — decomposition of positional signals in transcriptomic data.

Reference implementation of the three laws described in
"Shared cis-regulation drives local gene co-expression across human tissues".
"""
from .decompose import periodogram_identity, landscape_share, lag_covariance
from .laws import isochore_law, eqtl_law, saturation_exponent, fit_contact_law
from .simulate import simulate_genome

__version__ = "1.0.0"
__all__ = ["periodogram_identity", "landscape_share", "lag_covariance",
           "isochore_law", "eqtl_law", "saturation_exponent", "fit_contact_law",
           "simulate_genome"]
