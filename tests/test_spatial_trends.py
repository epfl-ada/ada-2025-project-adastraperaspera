import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from src.scripts.preprocessing.spatial_trends import (
    assign_distance_bins,
    mean_expression_by_bin,
    plot_gene_trends,
    plot_mean_heatmap,
)


@pytest.fixture
def fake_cells_df():
    """Simple fake cell dataset for testing."""
    np.random.seed(0)
    return pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(20)],
            "distance_to_plaque": np.linspace(0, 500, 20),
            "GeneA": np.random.rand(20),
            "GeneB": np.random.rand(20),
            "GeneC": np.random.rand(20),
        }
    )


def test_assign_distance_bins_creates_bins(fake_cells_df):
    df_binned = assign_distance_bins(fake_cells_df, n_bins=4)
    assert "distance_bin" in df_binned.columns
    assert isinstance(df_binned["distance_bin"].dtype, pd.CategoricalDtype)
    assert df_binned["distance_bin"].nunique() <= 4


def test_mean_expression_by_bin_computes_means(fake_cells_df):
    df_binned = assign_distance_bins(fake_cells_df, n_bins=4)
    grouped = mean_expression_by_bin(df_binned, ["GeneA", "GeneB"])
    assert isinstance(grouped, pd.DataFrame)
    assert {"GeneA", "GeneB"}.issubset(grouped.columns)

    assert grouped.shape[0] <= 4


def test_plot_gene_trends_executes(fake_cells_df):
    df_binned = assign_distance_bins(fake_cells_df, n_bins=3)
    grouped = mean_expression_by_bin(df_binned, ["GeneA", "GeneB"])

    plot_gene_trends(grouped, genes=["GeneA", "GeneB"])


def test_plot_mean_heatmap_executes(fake_cells_df):
    df_binned = assign_distance_bins(fake_cells_df, n_bins=3)
    grouped = mean_expression_by_bin(df_binned, ["GeneA", "GeneB", "GeneC"])

    plot_mean_heatmap(grouped, top_n=2)


def test_missing_distance_column_raises():
    df = pd.DataFrame({"GeneA": [1, 2, 3]})
    with pytest.raises(ValueError, match="distance_to_plaque"):
        assign_distance_bins(df)


def test_missing_distance_bin_for_mean_expr_raises(fake_cells_df):
    with pytest.raises(ValueError, match="distance_bin"):
        mean_expression_by_bin(fake_cells_df, ["GeneA"])
