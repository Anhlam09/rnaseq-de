import numpy as np
import pandas as pd

from rnaseq_de import filter_low_counts, log_normalize, size_factors


def test_size_factors_recover_known_scaling():
    rng = np.random.default_rng(0)
    base = rng.integers(20, 2000, size=(500, 1))
    scale = np.array([0.5, 1.0, 2.0, 4.0])
    counts = pd.DataFrame(base * scale, columns=list("abcd"))
    sf = size_factors(counts)
    np.testing.assert_allclose(sf / sf["b"], scale, rtol=1e-9)


def test_log_normalize_removes_library_size():
    counts = pd.DataFrame({"a": [10, 100, 1000], "b": [20, 200, 2000]})
    log_expr = log_normalize(counts, pseudocount=0)
    np.testing.assert_allclose(log_expr["a"], log_expr["b"])


def test_filter_low_counts():
    counts = pd.DataFrame(
        {"a": [0, 50, 9], "b": [0, 50, 50], "c": [1, 50, 50]}, index=["x", "y", "z"]
    )
    assert list(filter_low_counts(counts, min_count=10, min_samples=2).index) == ["y", "z"]
