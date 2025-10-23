"""
Gene-wise regression and correlation analyses versus plaque distance.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests

from src.utils.logging_utils import logger


def compute_gene_spatial_stats(
    df: pd.DataFrame,
    gene_cols: list[str],
    distance_col: str = "distance_to_plaque",
) -> pd.DataFrame:
    """
    For each gene, compute correlation and regression slope vs. plaque distance,
    plus FDR-corrected p-values.
    """
    if distance_col not in df.columns:
        raise ValueError(f"Missing distance column: {distance_col}")

    x = df[[distance_col]].to_numpy()
    results = []

    for g in gene_cols:
        y = df[g].to_numpy()
        if np.allclose(y, 0):
            continue  # skip constant genes

        # Spearman correlation
        r, p = spearmanr(x.ravel(), y)
        if np.isnan(r):
            continue

        # Linear regression slope
        model = LinearRegression().fit(x, y)
        slope = float(model.coef_[0])

        results.append((g, r, p, slope))

    res_df = pd.DataFrame(results, columns=["gene", "spearman_r", "spearman_p", "slope"])

    # --- Apply FDR correction ---
    if not res_df.empty:
        _, fdr, _, _ = multipletests(res_df["spearman_p"], method="fdr_bh")
        res_df["fdr_pval"] = fdr

    res_df.sort_values("spearman_r", ascending=False, inplace=True)
    logger.info(f"Computed spatial trends for {len(res_df)} genes (FDR-corrected).")
    return res_df


def plot_top_spatial_genes(
    stats_df: pd.DataFrame,
    top_n: int = 20,
    metric: str = "spearman_r",
    figsize: tuple = (6, 6),
) -> None:
    """Bar plot of genes most correlated with plaque distance."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    top = stats_df.nlargest(top_n, metric)
    plt.figure(figsize=figsize)
    sns.barplot(
        data=top,
        x=metric,
        y="gene",
        palette="vlag" if metric == "spearman_r" else "crest",
    )
    plt.title(f"Top {top_n} genes by {metric}")
    plt.xlabel(metric)
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.show()
