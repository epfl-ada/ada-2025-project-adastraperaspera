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
    grouped = grouped.sort_index(key=lambda x: x.map(lambda i: i.mid))  # order by bin midpoint
    return grouped


def plot_gene_trends(
    mean_expr: pd.DataFrame,
    genes: list[str],
    ylabel: str = "Mean expression (log1p normalized)",
    xlabel: str = "Distance to plaque (µm, binned)",
    figsize: tuple = (7, 4),
    title: str = "Spatial gene expression gradients",
) -> None:
    """Plot mean gene expression trends across distance bins."""
    if mean_expr.empty:
        raise ValueError("mean_expr is empty; check your inputs")

    plt.figure(figsize=figsize)
    x = range(len(mean_expr))
    bin_labels = [f"{b.left:.0f}-{b.right:.0f}" for b in mean_expr.index]

    found = False
    for gene in genes:
        if gene not in mean_expr.columns:
            logger.warning(f"[WARN] Gene '{gene}' not found, skipping.")
            continue
        plt.plot(x, mean_expr[gene], marker="o", linewidth=2, label=gene)
        found = True

    if not found:
        plt.text(0.5, 0.5, "No valid genes found", ha="center", va="center")
        plt.axis("off")
        plt.show()
        return

    plt.xticks(x, bin_labels, rotation=45, ha="right")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_mean_heatmap(
    mean_expr: pd.DataFrame,
    top_n: int = 25,
    zscore: bool = True,
) -> None:
    """
    Visualize top N genes with strongest spatial variation across plaque distance.
    """
    import matplotlib.pyplot as plt
    import pandas as pd

    non_gene_cols = {
        "cell_id",
        "x_centroid",
        "y_centroid",
        "cell_area",
        "nucleus_area",
        "total_counts",
        "transcript_counts",
        "distance_to_plaque",
        "distance_bin",
    }
    gene_cols = [c for c in mean_expr.columns if c not in non_gene_cols]
    mean_expr = mean_expr[gene_cols]

    grad = mean_expr.diff().abs().sum().sort_values(ascending=False)
    top_genes = grad.head(top_n).index
    sub_df = mean_expr[top_genes]

    # --- convert sparse columns to dense if needed ---
    if zscore:
        sub_df = sub_df.apply(lambda x: x.sparse.to_dense() if pd.api.types.is_sparse(x) else x)
        sub_df = (sub_df - sub_df.mean()) / sub_df.std(ddof=0)

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        sub_df.T,
        cmap="vlag" if zscore else "magma",
        center=0 if zscore else None,
        cbar_kws={"label": "Z-scored mean expression" if zscore else "Mean log1p expression"},
    )
    plt.title(f"Top {top_n} genes varying with plaque distance")
    plt.xlabel("Distance bin")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.show()
