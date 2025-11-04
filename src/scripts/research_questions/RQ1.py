from __future__ import annotations

from collections.abc import Sequence
import logging
import math
import re
from typing import Any

from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


def summarize_gene_stats(
    reg_df: pd.DataFrame,
    anova_df: pd.DataFrame,
    *,
    sort_column: str = "slope",
    ascending: bool = True,
    top_n: int = 10,
    bottom_n: int = 5,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Merge regression (continuous) and ANOVA (categorical) gene statistics, sort, and
    optionally display the top/bottom genes.

    This function replicates and generalizes the provided snippet:
      1) Renames overlapping p/q-value columns to avoid collisions.
      2) Left-joins ANOVA stats onto regression stats by the "gene" key.
      3) Sorts by the specified column (default: "slope").
      4) Selects a compact view of key columns if present.
      5) Optionally logs/displays the top and bottom rows.

    Parameters
    ----------
    reg_df : pd.DataFrame
        DataFrame containing continuous regression results. Must include columns:
        "gene" and the sort column (default "slope"). If present, "pval" and "qval"
        are renamed to "pval_continuous" and "qval_continuous".
    anova_df : pd.DataFrame
        DataFrame containing ANOVA results keyed by "gene". If present, "pval" and
        "qval" are renamed to "anova_pval" and "qval_anova" before merging.
    sort_column : str, optional
        Column name in the merged DataFrame to sort by (default "slope").
    ascending : bool, optional
        Sort order for `sort_column`. True for ascending (default), False for descending.
    top_n : int, optional
        Number of top rows to show/return from the sorted DataFrame (default 10).
    bottom_n : int, optional
        Number of bottom rows to show/return from the sorted DataFrame (default 5).
    logger : Optional[logging.Logger], optional
        Logger to use for info messages. If None, messages are not logged (default).

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        A tuple of:
        - summary_stats : The merged and sorted DataFrame.
        - top_df        : The first `top_n` rows of `summary_stats` with key columns.
        - bottom_df     : The last `bottom_n` rows of `summary_stats` with key columns.

    Raises
    ------
    ValueError
        If required columns ("gene" and `sort_column`) are missing from inputs.
    """
    # --- Validate required columns
    if "gene" not in reg_df.columns:
        raise ValueError("`reg_df` must contain a 'gene' column.")
    if sort_column not in reg_df.columns and sort_column not in anova_df.columns:
        raise ValueError(
            f"`sort_column='{sort_column}'` not found in inputs. "
            "It must exist in either reg_df or anova_df (preferably reg_df)."
        )

    # --- Rename overlapping columns to avoid collisions
    reg_renamed = reg_df.rename(
        columns={"pval": "pval_continuous", "qval": "qval_continuous"}
    )
    anova_renamed = anova_df.rename(
        columns={"pval": "anova_pval", "qval": "qval_anova"}
    )

    # --- Merge
    summary_stats = reg_renamed.merge(anova_renamed, on="gene", how="left")

    # --- Sort
    if sort_column not in summary_stats.columns:
        # If the requested sort column isn't in merged (unlikely), raise a clearer error.
        raise ValueError(
            f"`sort_column='{sort_column}'` not present after merge. "
            "Ensure it exists in reg_df or choose another column."
        )
    summary_stats = summary_stats.sort_values(sort_column, ascending=ascending)

    # --- Key columns to show (only those present)
    key_cols = [
        c
        for c in [
            "gene",
            "slope",
            "pval_continuous",
            "qval_continuous",
            "anova_pval",
            "qval_anova",
        ]
        if c in summary_stats.columns
    ]

    # --- Slices to return
    top_df = summary_stats.head(top_n)[key_cols]
    bottom_df = summary_stats.tail(bottom_n)[key_cols]

    # --- Optional logging/display
    if logger is not None:
        logger.info(
            "=== Top plaque-proximal genes (negative slope, smallest q-values) ==="
        )
        logger.info(top_df.head())

    if logger is not None:
        logger.info("\n=== Least plaque-responsive genes (weakest or flat slopes) ===")
        logger.info(bottom_df.head())

    return summary_stats, top_df, bottom_df


def compute_genewise_categorical_anova(
    data: pd.DataFrame,
    genes: Sequence[str],
    group_col: str = "distance_bin",
    *,
    anova_type: int = 2,
    multitest_method: str = "fdr_bh",
    min_groups: int = 2,
    raise_on_missing_gene: bool = False,
) -> pd.DataFrame:
    """
    Run one-way (categorical) ANOVA for each gene against a grouping variable.

    This function fits an OLS model for each gene using the formula
    ``gene_expression ~ C(group_col)`` and obtains a Type-II ANOVA table to
    extract the p-value for the group effect. Gene names are sanitized to be
    valid Patsy identifiers by replacing non-alphanumeric/underscore
    characters with underscores. Multiple testing correction is applied across
    genes using the specified method.

    Parameters
    ----------
    data : pd.DataFrame
        Input dataframe containing gene expression columns and the grouping column.
    genes : Sequence[str]
        List/sequence of gene column names to test.
    group_col : str, default "distance_bin"
        Name of the categorical grouping column in `data`.
    anova_type : int, default 2
        Type of sums of squares for ANOVA (passed to `sm.stats.anova_lm`).
        Common choices are 1, 2, or 3.
    multitest_method : str, default "fdr_bh"
        Multiple testing correction method passed to `statsmodels.stats.multitest.multipletests`.
        Examples: "bonferroni", "holm", "fdr_bh", "fdr_by".
    min_groups : int, default 2
        Minimum number of distinct groups required to attempt the ANOVA for a gene.
        Genes with fewer groups after dropping NA rows will receive NaN p-values.
    raise_on_missing_gene : bool, default False
        If True, raise a KeyError when a gene is not found in `data`. If False, skip it.

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per gene and the following columns:
        - ``gene``: original gene name
        - ``anova_pval``: p-value for the group effect from the ANOVA
        - ``qval``: multiple-testing–adjusted p-value (NaN if `anova_pval` is NaN)
        The rows are sorted by ``anova_pval`` ascending (NaNs last).
    """
    if group_col not in data.columns:
        raise KeyError(f"`group_col` '{group_col}' not found in `data`.")

    results: list[dict[str, float | str | None]] = []

    for gene in genes:
        if gene not in data.columns:
            if raise_on_missing_gene:
                raise KeyError(f"Gene '{gene}' not found in `data`.")
            # Skip missing genes silently if not raising
            continue

        # Subset and drop rows with NA in either column
        df_sub = data[[gene, group_col]].dropna().copy()

        # Require at least `min_groups` distinct categories to attempt ANOVA
        if df_sub[group_col].nunique(dropna=True) < min_groups:
            results.append({"gene": gene, "anova_pval": np.nan})
            continue

        # Sanitize the gene name for Patsy/formula usage
        safe_gene = re.sub(r"[^0-9a-zA-Z_]", "_", gene)
        df_sub = df_sub.rename(columns={gene: safe_gene})

        try:
            formula = f"{safe_gene} ~ C({group_col})"
            model = smf.ols(formula, data=df_sub).fit()
            aov_table = sm.stats.anova_lm(model, typ=anova_type)

            term = f"C({group_col})"
            if term in aov_table.index:
                pval = float(aov_table.loc[term, "PR(>F)"])
            else:
                # Fallback: take the first row if the expected term name isn't present
                pval = float(aov_table["PR(>F)"].iloc[0])

        except Exception:
            pval = np.nan

        results.append({"gene": gene, "anova_pval": pval})

    # Assemble results
    anova_df = pd.DataFrame(results, columns=["gene", "anova_pval"])

    if anova_df.empty:
        return anova_df

    # Multiple testing correction (handle NaNs safely)
    pvals = anova_df["anova_pval"].to_numpy(dtype=float)
    mask = np.isfinite(pvals)
    qvals = np.full_like(pvals, np.nan, dtype=float)
    if mask.any():
        _, q_corrected, _, _ = multipletests(pvals[mask], method=multitest_method)
        qvals[mask] = q_corrected

    anova_df["qval"] = qvals

    # Sort by raw p-value, placing NaNs at the end
    anova_df = anova_df.sort_values("anova_pval", na_position="last").reset_index(
        drop=True
    )

    return anova_df


def regress_expression_vs_distance(
    data: pd.DataFrame,
    gene_list: Sequence[str],
    distance_col: str = "distance_to_plaque",
    top_n: int = 10,
    fdr_method: str = "fdr_bh",
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run per-gene linear regressions of expression on distance and summarize results.

    For each gene present in `data`, fits the model:
        expression_gene ~ 1 + distance

    Collects the slope and p-value for the distance term, adjusts p-values
    via multiple testing correction, and returns a full results table as well
    as the top plaque-proximal genes (negative slope) ranked by adjusted p-value.

    Parameters
    ----------
    data : pd.DataFrame
        Input dataframe containing gene expression columns and a distance column.
    gene_list : Sequence[str]
        Candidate gene names (columns) to test; only those found in `data` are used.
    distance_col : str, optional
        Name of the distance column in `data`. Defaults to "distance_to_plaque".
    top_n : int, optional
        Number of top plaque-proximal genes (negative slope) to return. Defaults to 10.
    fdr_method : str, optional
        Multiple testing correction method passed to `statsmodels.stats.multitest.multipletests`.
        Defaults to "fdr_bh".
    logger : Optional[logging.Logger], optional
        If provided, brief result summaries are logged with `logger.info`.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - reg_df: Full results table with columns ["gene", "slope", "pval", "qval"],
          sorted by "slope" (ascending).
        - proximal_genes: Subset of reg_df with negative slopes, sorted by "qval",
          containing up to `top_n` rows.

    Raises
    ------
    ValueError
        If `distance_col` is not present in `data`.
    """
    if distance_col not in data.columns:
        raise ValueError(f"'{distance_col}' not found in data columns.")

    # Use only genes present in the dataframe
    pig_cols = [g for g in gene_list if g in data.columns]

    results = []
    for gene in pig_cols:
        df = data[[gene, distance_col]].dropna()
        # Need at least 2 observations to fit OLS with intercept and one predictor
        if len(df) < 2:
            continue

        # Design matrix with intercept
        X = sm.add_constant(df[distance_col])
        y = df[gene]

        try:
            model = sm.OLS(y, X).fit()
            slope = float(model.params[distance_col])
            pval = float(model.pvalues[distance_col])
            results.append({"gene": gene, "slope": slope, "pval": pval})
        except Exception:
            # Skip genes where the regression fails (e.g., singular matrix)
            continue

    # Convert to DataFrame
    reg_df = pd.DataFrame(results)
    if not reg_df.empty:
        reg_df["qval"] = multipletests(reg_df["pval"], method=fdr_method)[1]
        reg_df = reg_df.sort_values("slope", ascending=True, kind="mergesort")
    else:
        reg_df["qval"] = pd.Series(dtype=float)

    if logger is not None:
        logger.info("Continuous regression results (expression ~ distance).")
        logger.info("Computed %d gene regressions.", len(reg_df))

    # Highlight top plaque-proximal genes (negative slope)
    proximal_genes = (
        reg_df.loc[reg_df["slope"] < 0]
        .sort_values("qval", ascending=True, kind="mergesort")
        .head(top_n)
        .copy()
    )

    if logger is not None:
        logger.info(
            "Top plaque-proximal genes (negative slope, smallest q-values): %d returned.",
            len(proximal_genes),
        )

    return reg_df, proximal_genes


def summarize_pig_expression_by_distance(
    df: pd.DataFrame,
    pig_genes: Sequence[str],
    distance_col: str = "distance_to_plaque",
    q: int = 5,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Bin cells by distance to plaque and compute mean ± SEM of PIG gene expression per bin.

    This function:
      1) Assigns each row to a distance bin (in µm) using `pandas.cut`.
      2) Computes per-bin means and standard errors of the mean (SEM) for the provided PIG genes.
      3) Returns a long-form summary DataFrame and the count of cells in each distance bin.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing a distance column and gene expression columns.
    pig_genes : Sequence[str]
        Collection of gene names to summarize (columns expected to exist in `df`).
    distance_col : str, default "distance_to_plaque"
        Name of the column in `df` with distances (in µm) to be binned.
    bins : Sequence[float], default (0, 20, 50, 100, 200, np.inf)
        Monotonic sequence of bin edges (µm) passed to `pandas.cut`.
    labels : Sequence[str], default ("0–20 µm", "20–50 µm", "50–100 µm", "100–200 µm", ">200 µm")
        Labels for the resulting bins. Must be length `len(bins) - 1`.
    include_lowest : bool, default True
        Whether the first interval should be left-inclusive.

    Returns
    -------
    summary_df : pd.DataFrame
        Long-form table with columns:
        - 'distance_bin' (Categorical): the distance bin label
        - 'gene' (str): gene name
        - 'mean_expr' (float): mean expression within the bin
        - 'sem_expr' (float): standard error of the mean within the bin
    counts : pd.Series
        Count of rows (cells) in each distance bin, indexed by the bin labels and
        sorted by bin order.

    Raises
    ------
    KeyError
        If `distance_col` is not present in `df`.
    ValueError
        If none of the requested `pig_genes` are present in `df`.
        If `labels` length does not match `len(bins) - 1`.
    """
    # Basic validations
    if distance_col not in df.columns:
        raise KeyError(f"Column '{distance_col}' not found in DataFrame.")

    pig_cols = [g for g in pig_genes if g in df.columns]
    if not pig_cols:
        raise ValueError(
            "None of the specified `pig_genes` are present in the DataFrame."
        )

    # Work on a copy to avoid mutating the caller's DataFrame
    tmp = df.copy()

    # Assign distance bins via equal quantiles
    tmp["distance_bin"] = pd.qcut(
        tmp[distance_col].astype(float), q=q, duplicates="drop"
    )

    # Compute mean and SEM for each distance bin
    mean_expr = tmp.groupby("distance_bin")[pig_cols].mean().reset_index()
    sem_expr = tmp.groupby("distance_bin")[pig_cols].sem().reset_index()

    # Melt to long-form and merge mean/sem
    summary_df = mean_expr.melt(
        id_vars="distance_bin", var_name="gene", value_name="mean_expr"
    )
    sem_melted = sem_expr.melt(
        id_vars="distance_bin", var_name="gene", value_name="sem_expr"
    )
    summary_df = summary_df.merge(sem_melted, on=["distance_bin", "gene"], how="left")

    # Cell counts per bin (sorted by categorical order)
    counts = tmp["distance_bin"].value_counts().sort_index()

    return summary_df, counts
