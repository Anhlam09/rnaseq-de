import numpy as np
import pandas as pd


def simulate_counts(
    n_genes: int = 5000,
    n_per_group: int = 4,
    frac_de: float = 0.1,
    lfc_sd: float = 1.5,
    dispersion: float = 0.1,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Simulate a two-condition RNA-seq experiment with negative-binomial counts.

    Returns (counts genes x samples, sample metadata, per-gene ground truth).
    """
    rng = np.random.default_rng(seed)

    base_mean = rng.lognormal(mean=4.0, sigma=1.8, size=n_genes)
    is_de = np.zeros(n_genes, dtype=bool)
    is_de[rng.choice(n_genes, int(frac_de * n_genes), replace=False)] = True
    true_lfc = np.where(is_de, rng.normal(0.0, lfc_sd, n_genes), 0.0)

    conditions = np.array(["control"] * n_per_group + ["treated"] * n_per_group)
    samples = [f"{c}_{i % n_per_group + 1}" for i, c in enumerate(conditions)]
    # Library-size variation is what normalization has to undo.
    library_scale = rng.uniform(0.6, 1.6, size=len(samples))

    treated = conditions == "treated"
    mu = (
        base_mean[:, None]
        * library_scale[None, :]
        * np.where(treated, 2.0 ** true_lfc[:, None], 1.0)
    )

    # numpy's NB is parameterized by (r, p); var = mu + dispersion * mu^2 when r = 1/dispersion.
    r = 1.0 / dispersion
    counts = rng.negative_binomial(r, r / (r + mu))

    genes = [f"gene_{i:05d}" for i in range(n_genes)]
    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    metadata = pd.DataFrame(
        {"condition": conditions, "library_scale": library_scale}, index=samples
    )
    truth = pd.DataFrame(
        {"is_de": is_de, "true_lfc": true_lfc, "base_mean": base_mean}, index=genes
    )
    return counts_df, metadata, truth
