"""Shared helper: per-tumour positional score at given lags (used by cohesin_freedman_lane.py and tumour_baseline_sensitivity.py)."""
import numpy as np
from lib_tumour import CHR


def lagscore(D, chrs, lags):
    Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
    return np.mean([np.concatenate([Z[x[:-L]] * Z[x[L:]] for c in CHR for x in [np.where(chrs == c)[0]] if len(x) > L + 5], 0).mean(0) for L in lags], 0)
