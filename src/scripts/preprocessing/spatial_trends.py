"""
Functions to quantify and visualize spatial trends of gene expression
relative to amyloid-beta plaque distance.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src.utils.logging_utils import logger


def assign_distance_bins(df: pd.DataFrame, n_bins: int = 5) -> pd.DataFrame:
    """Discretize 'distance_to_plaque' into quantile bins."""
    if "distance_to_plaque" not in df.columns:
        raise ValueError("Missing required column: distance_to_plaque")

    logger.info(f"Assigning {n_bins} quantile bins by distance_to_plaque")

    df = df.copy()
    df["distance_bin"] = pd.qcut(df["distance_to_plaque"], q=n_bins, duplicates="drop")
    return df


def mean_expression_by_bin(df: pd.DataFrame, gene_cols: list[str]) -> pd.DataFrame:
    """Compute mean expression per gene across distance bins."""
    if "distance_bin" not in df.columns:
        raise ValueError("distance_bin column missing; run assign_distance_bins first")

    logger.info("Computing mean expression per bin")

    grouped = df.groupby("distance_bin", observed=False)[gene_cols].mean()
    grouped = grouped.sort_index(key=lambda x: x.map(lambda i: i.mid))
    return grouped


def sem_expression_by_bin(df: pd.DataFrame, gene_cols: list[str]) -> pd.DataFrame:
    """
    Compute SEM (standard error of the mean) per gene across distance bins.

    Requires:
        - df['distance_bin'] created by assign_distance_bins
        - gene_cols: list of gene columns (log1p-normalized if that is the convention)

    Returns:
        DataFrame indexed by ordered distance_bin, columns = gene_cols, values = SEM.
    """
    if "distance_bin" not in df.columns:
        raise ValueError("distance_bin column missing; run assign_distance_bins first")

    logger.info("Computing SEM per bin")

    grouped = df.groupby("distance_bin", observed=False)[gene_cols].sem()

    grouped = grouped.sort_index(key=lambda x: x.map(lambda i: i.mid))
    return grouped
