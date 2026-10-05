"""Shared-source covariance (Supplementary Note 1, section 3)."""
import numpy as np, pytest
from poslayers import source_covariance


def test_two_shared_sources_add():
    assert source_covariance([2.0, 0.5], [1.0, 0.4], [0.5, -1.0]) == pytest.approx(2.0 * 0.5 + 0.5 * 0.4 * -1.0)


def test_matches_simulation():
    rng = np.random.default_rng(0); n = 200_000
    A = rng.normal(size=(n, 2)) * np.sqrt([1.5, 0.6]); fi, fj = np.array([0.8, 0.3]), np.array([0.5, -0.7])
    xi = A @ fi + rng.normal(size=n); xj = A @ fj + rng.normal(size=n)
    assert np.cov(xi, xj)[0, 1] == pytest.approx(source_covariance([1.5, 0.6], fi, fj), abs=0.01)


def test_shape_check():
    with pytest.raises(ValueError):
        source_covariance([1.0, 2.0], [1.0], [1.0, 1.0])
