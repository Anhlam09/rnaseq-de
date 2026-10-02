import numpy as np
import pandas as pd


def filter_low_counts(
    counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3
) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2-style median-of-ratios size factors."""
    expressed = counts.loc[(counts > 0).all(axis=1)]
    log_counts = np.log(expressed)
    log_geo_means = log_counts.mean(axis=1)
    return np.exp(log_counts.sub(log_geo_means, axis=0).median(axis=0))


def log_normalize(counts: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    return np.log2(counts / size_factors(counts) + pseudocount)
