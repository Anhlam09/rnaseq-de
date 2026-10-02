from rnaseq_de.cli import main


def test_cli_writes_outputs(tmp_path, capsys):
    main(["--n-genes", "500", "--seed", "3", "--out", str(tmp_path)])
    for name in ("de_results.csv", "volcano.png", "ma_plot.png", "pca.png"):
        assert (tmp_path / name).stat().st_size > 0
    assert "against simulated ground truth" in capsys.readouterr().out
