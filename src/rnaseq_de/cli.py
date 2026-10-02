import argparse
from pathlib import Path

import pandas as pd

from rnaseq_de import plots
from rnaseq_de.normalize import filter_low_counts, log_normalize
from rnaseq_de.simulate import simulate_counts
from rnaseq_de.stats import differential_expression, evaluate_calls


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Differential expression on RNA-seq counts.")
    parser.add_argument(
        "--counts", type=Path, help="CSV of counts (genes x samples). Simulated if omitted."
    )
    parser.add_argument(
        "--metadata", type=Path, help="CSV with sample index and a 'condition' column."
    )
    parser.add_argument("--n-genes", type=int, default=5000)
    parser.add_argument("--n-per-group", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--out", type=Path, default=Path("results"))
    args = parser.parse_args(argv)

    truth = None
    if args.counts:
        if not args.metadata:
            parser.error("--metadata is required with --counts")
        counts = pd.read_csv(args.counts, index_col=0)
        metadata = pd.read_csv(args.metadata, index_col=0)
    else:
        counts, metadata, truth = simulate_counts(
            n_genes=args.n_genes, n_per_group=args.n_per_group, seed=args.seed
        )

    filtered = filter_low_counts(counts)
    log_expr = log_normalize(filtered)
    results = differential_expression(log_expr, metadata["condition"])

    args.out.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out / "de_results.csv")
    plots.volcano(results, args.out / "volcano.png", alpha=args.alpha)
    plots.ma_plot(results, args.out / "ma_plot.png", alpha=args.alpha)
    plots.pca(log_expr, metadata["condition"], args.out / "pca.png")

    print(f"genes: {len(counts)} total, {len(filtered)} after filtering")
    print(f"significant at padj < {args.alpha}: {(results['padj'] < args.alpha).sum()}")
    print(results.head(10).round(4).to_string())
    if truth is not None:
        metrics = evaluate_calls(results, truth, alpha=args.alpha)
        print("\nagainst simulated ground truth:")
        for k, v in metrics.items():
            print(f"  {k:>16}: {v:.3f}" if isinstance(v, float) else f"  {k:>16}: {v}")
    print(f"\nwrote results to {args.out}/")


if __name__ == "__main__":
    main()
