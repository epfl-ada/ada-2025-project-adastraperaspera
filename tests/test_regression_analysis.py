import numpy as np
import pandas as pd
import pytest

from src.scripts.analysis.regression_analysis import (
    compute_gene_spatial_stats,
)


@pytest.fixture
def fake_data():
    np.random.seed(0)
    n = 50
    distance = np.linspace(0, 100, n)
    df = pd.DataFrame(
        {
            "distance_to_plaque": distance,
            "GeneA": distance * 0.5 + np.random.randn(n) * 2,  # positive trend
            "GeneB": -distance * 0.3 + np.random.randn(n) * 2,  # negative trend
            "GeneC": np.random.randn(n),  # random noise
        }
    )
    return df


def test_compute_gene_spatial_stats(fake_data):
    gene_cols = ["GeneA", "GeneB", "GeneC"]
    stats_df = compute_gene_spatial_stats(fake_data, gene_cols)
    assert {"gene", "spearman_r", "spearman_p", "slope"}.issubset(stats_df.columns)
    assert "GeneA" in stats_df["gene"].values

    # sanity checks: GeneA slope positive, GeneB negative
    slope_dict = dict(zip(stats_df["gene"], stats_df["slope"], strict=False))
    assert slope_dict["GeneA"] > 0
    assert slope_dict["GeneB"] < 0


def test_compute_gene_spatial_stats_includes_fdr(fake_data):
    stats_df = compute_gene_spatial_stats(fake_data, ["GeneA", "GeneB", "GeneC"])
    assert "fdr_pval" in stats_df.columns
    assert (stats_df["fdr_pval"] >= 0).all()
    assert (stats_df["fdr_pval"] <= 1).all()
