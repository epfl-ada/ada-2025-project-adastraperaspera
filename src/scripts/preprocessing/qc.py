"""
Quality control (QC) and normalization functions for Xenium single-cell data.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from src.utils.logging_utils import logger


def filter_cells(
    cells_df: pd.DataFrame,
    transcript_q: float = 0.05,
    area_q: float = 0.05,
    nucleus_q: float = 0.05,
    genes_q: float = 0.05,
    gene_cols: list[str] = [],
) -> pd.DataFrame:
    """
    Filter out low-quality cells using relative thresholds (percentile-based).

    Removes cells in the bottom `transcript_q` quantile for transcript counts
    or the bottom `area_q` quantile for cell area.

    Args:
        cells_df (pd.DataFrame): Input cell metadata table.
        transcript_q (float): Quantile for transcript cutoff (e.g., 0.05 keeps top 95%).
        area_q (float): Quantile for area cutoff (e.g., 0.05 keeps top 95%).

    Returns:
        pd.DataFrame: Filtered DataFrame of high-quality cells.
    """
    logger.info(
        f"Filtering cells below {transcript_q*100:.1f}th percentile for transcript_counts\n"
        f"and {area_q*100:.1f}th percentile for cell_area\n"
        f"and {nucleus_q*100:.1f}th percentile for nucleus_area\n"
        f"and {genes_q*100:.1f}th percentile for n_genes\n"
        f"and removing cells with <=0 cell area / nucleus area / total number of nonzero genes\n"
    )

    if (
        "transcript_counts" not in cells_df.columns
        or "cell_area" not in cells_df.columns
        or "nucleus_area" not in cells_df.columns
    ):
        raise KeyError(
            "Required columns 'transcript_counts', 'cell_area' and 'nucleus_area' not found."
        )
    cells_df = cells_df.loc[
        (cells_df["nucleus_area"] > 0) & (cells_df["cell_area"] > 0) & (cells_df["n_genes"] > 0)
    ].copy()

    min_transcripts = cells_df["transcript_counts"].quantile(transcript_q)
    min_area = cells_df["cell_area"].quantile(area_q)
    min_nucleus = cells_df["nucleus_area"].quantile(nucleus_q)
    min_genes = cells_df["n_genes"].quantile(genes_q)

    logger.info(
        f"Computed thresholds:\n"
        f"transcript_counts > {min_transcripts:.2f}\n"
        f"cell_area > {min_area:.2f}\n"
        f"nucleus_area > {min_nucleus:.2f}\n"
        f"n_genes > {min_genes:.2f}\n"
    )

    mask = (
        (cells_df["transcript_counts"] > min_transcripts)
        & (cells_df["cell_area"] > min_area)
        & (cells_df["nucleus_area"] > min_nucleus)
        & (cells_df["n_genes"] > min_genes)
    )
    filtered = cells_df.loc[mask].copy()

    logger.info(
        f"Retained {filtered.shape[0]} / {cells_df.shape[0]} cells "
        f"({filtered.shape[0] / cells_df.shape[0] * 100:.1f}%)"
    )
    return filtered


def normalize_expression(adata: sc.AnnData) -> sc.AnnData:
    """
    Normalize and log-transform expression counts in an AnnData object.

    Steps:
      1. Normalize total counts per cell to 1e4.
      2. Apply log1p transformation.
      3. Return a copy with normalized data.

    Args:
        adata (sc.AnnData): Raw expression data.

    Returns:
        sc.AnnData: Normalized AnnData object.
    """
    logger.info("⚙️ Normalizing expression matrix (target_sum=1e4)...")
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logger.info("Normalization complete.")
    return adata


def compute_nonzero_fraction(adata: sc.AnnData) -> float:
    """
    Compute the fraction of non-zero entries in the expression matrix.

    Args:
        adata (sc.AnnData): AnnData object (raw or normalized).

    Returns:
        float: Fraction of non-zero values.
    """
    logger.info("Computing nonzero fraction in expression matrix...")
    x = adata.X
    if hasattr(x, "nnz"):
        fraction = x.nnz / (x.shape[0] * x.shape[1])
    else:
        fraction = np.count_nonzero(x) / x.size
    logger.info(f"Nonzero fraction: {fraction:.3f}")
    return float(fraction)


def plot_qc_distributions(
    cells_df: pd.DataFrame,
    transcript_q: float = 0.05,
    area_q: float = 0.05,
    figsize: tuple[int, int] = (10, 4),
) -> None:
    """
    Plot QC histograms for transcript counts and cell area with percentile thresholds.

    Args:
        cells_df (pd.DataFrame): DataFrame containing 'transcript_counts' and 'cell_area'.
        transcript_q (float): Quantile for transcript cutoff (e.g., 0.05 keeps top 95%).
        area_q (float): Quantile for area cutoff (e.g., 0.05 keeps top 95%).
        figsize (tuple[int, int]): Figure size.
    """
    if not {"transcript_counts", "cell_area"}.issubset(cells_df.columns):
        raise KeyError("Missing required columns in cells_df")

    thr_trans = cells_df["transcript_counts"].quantile(transcript_q)
    thr_area = cells_df["cell_area"].quantile(area_q)
    logger.info(
        f"Plotting QC distributions with cutoffs: "
        f"transcripts > {thr_trans:.2f}, area > {thr_area:.2f}"
    )

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Transcript counts
    sns.histplot(cells_df["transcript_counts"], bins=100, ax=axes[0], color="lightgray")
    axes[0].axvline(thr_trans, color="red", linestyle="--", label=f"{transcript_q*100:.0f}th pct")
    axes[0].set_title("Transcript counts per cell")
    axes[0].set_xlabel("transcript_counts")
    axes[0].legend()
    axes[0].set_xscale("log")
    axes[1].set_xscale("log")

    # Cell area
    sns.histplot(cells_df["cell_area"], bins=100, ax=axes[1], color="lightgray")
    axes[1].axvline(thr_area, color="red", linestyle="--", label=f"{area_q*100:.0f}th pct")
    axes[1].set_title("Cell area distribution")
    axes[1].set_xlabel("cell_area")
    axes[1].legend()

    plt.tight_layout()
    plt.show()
