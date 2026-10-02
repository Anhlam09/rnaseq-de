from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SIG_COLOR = "#c2410c"
NS_COLOR = "#9ca3af"


def volcano(results: pd.DataFrame, path: Path, alpha: float = 0.05, lfc: float = 1.0) -> None:
    sig = (results["padj"] < alpha) & (results["log2fc"].abs() >= lfc)
    fig, ax = plt.subplots(figsize=(6, 5))
    y = -np.log10(results["p_value"].clip(lower=1e-300))
    ax.scatter(results.loc[~sig, "log2fc"], y[~sig], s=4, c=NS_COLOR, alpha=0.5, label="n.s.")
    ax.scatter(
        results.loc[sig, "log2fc"], y[sig], s=6, c=SIG_COLOR, label=f"padj<{alpha}, |LFC|≥{lfc}"
    )
    for x in (-lfc, lfc):
        ax.axvline(x, ls="--", lw=0.8, c="k", alpha=0.4)
    ax.set_xlabel("log2 fold change (treated / control)")
    ax.set_ylabel("-log10 p-value")
    ax.set_title("Volcano plot")
    ax.legend(frameon=False, loc="upper left")
    _save(fig, path)


def ma_plot(results: pd.DataFrame, path: Path, alpha: float = 0.05) -> None:
    sig = results["padj"] < alpha
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        results.loc[~sig, "mean_log_expr"], results.loc[~sig, "log2fc"], s=4, c=NS_COLOR, alpha=0.5
    )
    ax.scatter(results.loc[sig, "mean_log_expr"], results.loc[sig, "log2fc"], s=6, c=SIG_COLOR)
    ax.axhline(0, lw=0.8, c="k")
    ax.set_xlabel("mean log2 normalized expression")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot")
    _save(fig, path)


def pca(log_expr: pd.DataFrame, conditions: pd.Series, path: Path, top_n: int = 500) -> None:
    top = log_expr.loc[log_expr.var(axis=1).nlargest(top_n).index]
    centered = (top.T - top.T.mean()).to_numpy()
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    explained = s**2 / (s**2).sum()

    fig, ax = plt.subplots(figsize=(5.5, 5))
    for cond, color in zip(sorted(conditions.unique()), ("#2563eb", SIG_COLOR), strict=False):
        mask = (conditions[log_expr.columns] == cond).to_numpy()
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=60, c=color, label=cond, edgecolor="k")
    ax.set_xlabel(f"PC1 ({explained[0]:.0%})")
    ax.set_ylabel(f"PC2 ({explained[1]:.0%})")
    ax.set_title(f"PCA, top {top_n} variable genes")
    ax.legend(frameon=False)
    _save(fig, path)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
