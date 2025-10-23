"""
Unit tests for src.data.loaders module.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from src.scripts.preprocessing.loaders import (
    expression_to_dataframe,
    load_cells_table,
    load_expression_matrix,
    merge_expression_with_cells,
)


def test_load_cells_table_reads_parquet(tmp_path: Path) -> None:
    """Ensure load_cells_table reads a parquet file correctly."""

    df = pd.DataFrame({"cell_id": ["A", "B"], "transcript_counts": [100, 200]})
    p = tmp_path / "cells.parquet"
    df.to_parquet(p)

    loaded = load_cells_table(str(p))
    assert isinstance(loaded, pd.DataFrame)
    assert loaded.equals(df)
    assert "cell_id" in loaded.columns


def test_load_cells_table_raises_for_missing(tmp_path: Path) -> None:
    """Ensure missing parquet file raises an informative error."""
    with pytest.raises(FileNotFoundError, match="Cells table not found"):
        load_cells_table(str(tmp_path / "missing.parquet"))


def test_load_expression_matrix_reads_h5(tmp_path: Path) -> None:
    """Ensure load_expression_matrix handles non-10x H5 files gracefully."""

    # Create minimal synthetic AnnData and save as .h5ad (not a valid 10x-format file)
    adata = sc.AnnData(
        X=np.array([[1, 0], [0, 1]]),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["geneA", "geneB"]),
    )
    test_path = tmp_path / "matrix.h5ad"
    adata.write_h5ad(test_path)

    # Expect a ValueError or OSError due to invalid 10x H5 structure
    with pytest.raises((ValueError, OSError)):
        _ = load_expression_matrix(str(test_path))


def _fake_adata() -> sc.AnnData:
    x = np.array([[1, 2], [3, 4], [5, 6]])
    obs = pd.DataFrame(index=["A", "B", "C"])
    var = pd.DataFrame(index=["Gene1", "Gene2"])
    return sc.AnnData(X=x, obs=obs, var=var)


def test_expression_to_dataframe_returns_clean_df():
    adata = _fake_adata()
    df = expression_to_dataframe(adata)
    assert isinstance(df, pd.DataFrame)
    assert "Gene1" in df.columns
    assert df.shape == (3, 2)
    assert df.index.name == "cell_id"


def test_merge_expression_with_cells_inner_join(tmp_path):
    adata = _fake_adata()
    cells_df = pd.DataFrame(
        {
            "cell_id": ["A", "B", "X"],
            "x_centroid": [1, 2, 3],
            "distance_to_plaque": [10, 20, 30],
        }
    )

    merged = merge_expression_with_cells(cells_df, adata)
    assert "Gene1" in merged.columns
    # Only common IDs (A,B)
    assert set(merged["cell_id"]) == {"A", "B"}
    assert merged.shape[0] == 2
