import numpy as np
import pandas as pd
from scipy import stats


def differential_expression(
    log_expr: pd.DataFrame,
    conditions: pd.Series,
    reference: str = "control",
    treatment: str = "treated",
) -> pd.DataFrame:
    """Per-gene Welch t-test on log2 normalized expression with BH-adjusted p-values."""
    ref = log_expr.loc[:, conditions[log_expr.columns] == reference]
    trt = log_expr.loc[:, conditions[log_expr.columns] == treatment]

    t_stat, p_value = stats.ttest_ind(trt, ref, axis=1, equal_var=False)
    # Genes with zero variance in both groups produce NaN; treat them as uninformative.
    p_value = np.nan_to_num(p_value, nan=1.0)

    result = pd.DataFrame(
        {
            "log2fc": trt.mean(axis=1) - ref.mean(axis=1),
            "mean_log_expr": log_expr.mean(axis=1),
            "t_stat": t_stat,
            "p_value": p_value,
            "padj": stats.false_discovery_control(p_value, method="bh"),
        },
        index=log_expr.index,
    )
    return result.sort_values("p_value")


def evaluate_calls(
    results: pd.DataFrame, truth: pd.DataFrame, alpha: float = 0.05, min_abs_lfc: float = 0.0
) -> dict[str, float]:
    called = (results["padj"] < alpha) & (results["log2fc"].abs() >= min_abs_lfc)
    actual = truth.loc[results.index, "is_de"]
    tp = int((called & actual).sum())
    fp = int((called & ~actual).sum())
    fn = int((~called & actual).sum())
    return {
        "called": int(called.sum()),
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "precision": tp / (tp + fp) if tp + fp else float("nan"),
        "recall": tp / (tp + fn) if tp + fn else float("nan"),
        "observed_fdr": fp / (tp + fp) if tp + fp else 0.0,
    }
