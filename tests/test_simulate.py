import numpy as np

from rnaseq_de import simulate_counts


def test_shapes_and_labels():
    counts, metadata, truth = simulate_counts(n_genes=300, n_per_group=3, frac_de=0.2, seed=1)
    assert counts.shape == (300, 6)
    assert list(metadata.index) == list(counts.columns)
    assert (metadata["condition"] == "treated").sum() == 3
    assert truth["is_de"].sum() == 60
    assert (truth.loc[~truth["is_de"], "true_lfc"] == 0).all()


def test_reproducible_with_seed():
    a, _, _ = simulate_counts(n_genes=100, seed=7)
    b, _, _ = simulate_counts(n_genes=100, seed=7)
    c, _, _ = simulate_counts(n_genes=100, seed=8)
    assert a.equals(b)
    assert not np.array_equal(a.to_numpy(), c.to_numpy())
