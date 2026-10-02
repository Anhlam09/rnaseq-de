import numpy as np
import pandas as pd

from rnaseq_de import differential_expression, evaluate_calls, log_normalize, simulate_counts
from rnaseq_de.normalize import filter_low_counts


def test_detects_planted_signal():
    conditions = pd.Series(["control"] * 3 + ["treated"] * 3, index=[f"s{i}" for i in range(6)])
    expr = pd.DataFrame(
        [[5.0, 5.1, 4.9, 8.0, 8.1, 7.9], [5.0, 5.2, 4.8, 5.1, 4.9, 5.0]],
        index=["up", "flat"],
        columns=conditions.index,
    )
    res = differential_expression(expr, conditions)
    assert res.index[0] == "up"
    assert res.loc["up", "log2fc"] == np.float64(3.0)
    assert res.loc["up", "p_value"] < 0.001
    assert res.loc["flat", "p_value"] > 0.5


def test_end_to_end_controls_fdr_on_simulation():
    counts, metadata, truth = simulate_counts(n_genes=3000, n_per_group=5, seed=42)
    results = differential_expression(
        log_normalize(filter_low_counts(counts)), metadata["condition"]
    )
    metrics = evaluate_calls(results, truth, alpha=0.05)
    assert metrics["called"] > 50
    assert metrics["observed_fdr"] < 0.15
    assert results["padj"].between(0, 1).all()
