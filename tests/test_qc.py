"""
Unit tests for src.preprocessing.qc module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from src.scripts.preprocessing.qc import (
    compute_nonzero_fraction,
    filter_cells,
    normalize_expression,
)


def test_filter_cells_relative_quantiles() -> None:
    """Check that percentile-based filtering adapts dynamically."""
    df = pd.DataFrame(
        {
            "cell_id": [f"C{i}" for i in range(10)],
            "transcript_counts": np.linspace(0, 100, 10),
            "cell_area": np.linspace(10, 100, 10),
        }
    )
    # Keep top 80% → should drop bottom 2 cells
    filtered = filter_cells(df, transcript_q=0.2, area_q=0.2)
    assert len(filtered) == 8
    # Ensure the lowest transcript_counts and areas are gone
    assert filtered["transcript_counts"].min() > df["transcript_counts"].quantile(0.2)


def test_normalize_expression_and_nonzero_fraction() -> None:
    """Check normalization and sparsity computation."""
    x = np.array([[1, 0], [3, 4]], dtype=float)
    adata = sc.AnnData(X=x)
    normed = normalize_expression(adata)
    # After normalization and log1p, should remain nonnegative
    assert (normed.X >= 0).all()
    frac = compute_nonzero_fraction(normed)
    assert 0 < frac <= 1


def test_filter_cells_all_fail_safe() -> None:
    """Ensure filtering works safely even if few cells exist."""
    df = pd.DataFrame(
        {
            "cell_id": ["A", "B"],
            "transcript_counts": [1, 2],
            "cell_area": [10, 15],
        }
    )
    filtered = filter_cells(df, transcript_q=0.5, area_q=0.5)
    assert isinstance(filtered, pd.DataFrame)
    assert 0 <= len(filtered) <= len(df)
