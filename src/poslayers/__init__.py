"""poslayers — decomposition of positional signals in transcriptomic data.

Reference implementation of the decomposition and the model relations described in
"Common cis-regulatory inputs shape local gene co-expression across human tissues".
"""
from .decompose import (periodogram_identity, landscape_share, lag_covariance,
                        lag_profile_linear, gc_correct, gc_slopes)
from .laws import (gc_layer_correlation, predicted_genetic_covariance, source_covariance,
                   isochore_law, eqtl_law, saturation_exponent, fit_contact_law)   # second line: release 0.1.0 names
from .simulate import simulate_genome

__version__ = "0.1.0"
__all__ = ["periodogram_identity", "landscape_share", "lag_covariance", "lag_profile_linear", "gc_correct", "gc_slopes",
           "gc_layer_correlation", "predicted_genetic_covariance", "source_covariance", "simulate_genome",
           "isochore_law", "eqtl_law", "saturation_exponent", "fit_contact_law"]
