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
    reg_renamed = reg_df.rename(columns={"pval": "pval_continuous", "qval": "qval_continuous"})
    anova_renamed = anova_df.rename(columns={"pval": "anova_pval", "qval": "qval_anova"})

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
        logger.info("=== Top plaque-proximal genes (negative slope, smallest q-values) ===")
        logger.info(top_df.head())

    if logger is not None:
        logger.info("\n=== Least plaque-responsive genes (weakest or flat slopes) ===")
        logger.info(bottom_df.head())

    return summary_stats, top_df, bottom_df


def plot_top_genes_by_distance(
    anova_df: pd.DataFrame,
    combined_df_normalized: pd.DataFrame,
    *,
    top_k: int = 8,
    n_cols: int = 4,
    gene_col: str = "gene",
    qval_col: str = "qval",
    distance_col: str = "distance_bin",
    y_label: str = "Expression (log1p)",
    palette: str = "viridis",
    max_points: int = 1000,
    jitter: bool = True,
    jitter_size: float = 2.0,
    jitter_alpha: float = 0.3,
    random_state: int = 0,
    sharey: bool = False,
    fig_width: float = 14.0,
    fig_height_per_row: float = 3.0,
    suptitle: str = "Expression by Distance Bin for Top Plaque-Responsive Genes",
    show: bool = True,
    preselected_genes: list[str] | None = None,
) -> tuple[Figure, list[Axes], list[str]]:
    """
    Plot expression distributions by distance bins for the top genes from an ANOVA table.

    This function:
      1) Selects the top `top_k` genes from `anova_df` by ascending `qval_col`
         (or uses `preselected_genes` if provided),
      2) Creates a grid of boxplots across `distance_col` for each gene,
      3) Overlays a jittered strip of up to `max_points` nonzero observations per gene.

    Parameters
    ----------
    anova_df : pd.DataFrame
        DataFrame containing at least the columns specified by `gene_col` and `qval_col`.
    combined_df_normalized : pd.DataFrame
        Long or wide expression matrix that includes `distance_col` and one column per gene to plot.
        Each gene column should be numeric and represent log1p-normalized expression values.
    top_k : int, optional
        Number of top genes (by `qval_col`) to plot. Ignored if `preselected_genes` is provided.
    n_cols : int, optional
        Number of columns in the subplot grid.
    gene_col : str, optional
        Column in `anova_df` that contains gene names.
    qval_col : str, optional
        Column in `anova_df` that contains q-values for ranking genes.
    distance_col : str, optional
        Column in `combined_df_normalized` indicating distance bins (categorical or discrete).
    y_label : str, optional
        Y-axis label for expression.
    palette : str, optional
        Seaborn palette for the boxplots (e.g., "viridis").
    max_points : int, optional
        Maximum number of nonzero points to overlay per gene for the stripplot.
    jitter : bool, optional
        Whether to jitter the overlaid points.
    jitter_size : float, optional
        Marker size for the overlaid points.
    jitter_alpha : float, optional
        Transparency for the overlaid points.
    random_state : int, optional
        Random seed for sampling overlaid points.
    sharey : bool, optional
        Whether subplots share the y-axis scale.
    fig_width : float, optional
        Total figure width in inches.
    fig_height_per_row : float, optional
        Figure height per row in inches.
    suptitle : str, optional
        Figure title shown above all subplots.
    show : bool, optional
        If True, calls `plt.show()` before returning.
    preselected_genes : list of str, optional
        If provided, use this exact list of genes instead of picking from `anova_df`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    axes_list : list of matplotlib.axes.Axes
        Flattened list of axes corresponding to the subplots (unused axes are hidden).
    top_genes : list of str
        The list of genes that were plotted, in order.

    Raises
    ------
    ValueError
        If required columns are missing, or if `top_k` <= 0 and no `preselected_genes` provided,
        or if no genes are available to plot.
    """
    # ---- Basic validation
    for col, df_name, df in [
        (gene_col, "anova_df", anova_df),
        (qval_col, "anova_df", anova_df),
        (distance_col, "combined_df_normalized", combined_df_normalized),
    ]:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in {df_name}.")

    if preselected_genes is not None:
        top_genes = list(preselected_genes)
    else:
        if top_k <= 0:
            raise ValueError("`top_k` must be > 0 when `preselected_genes` is not provided.")
        if gene_col not in anova_df.columns or qval_col not in anova_df.columns:
            raise ValueError(f"`anova_df` must contain '{gene_col}' and '{qval_col}'.")
        top_genes = (
            anova_df.sort_values(qval_col, ascending=True)
            .head(top_k)[gene_col]
            .astype(str)
            .tolist()
        )

    if len(top_genes) == 0:
        raise ValueError("No genes available to plot. Check inputs or `preselected_genes`.")

    # Ensure all genes exist in data
    missing = [g for g in top_genes if g not in combined_df_normalized.columns]
    if missing:
        raise ValueError(
            f"The following genes are missing from `combined_df_normalized`: {missing}"
        )

    # ---- Layout
    n_rows = max(1, math.ceil(len(top_genes) / n_cols))
    fig_height = n_rows * float(fig_height_per_row)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(float(fig_width), fig_height), sharey=sharey)

    # Normalize axes to a flat list
    if isinstance(axes, np.ndarray):
        axes_flat = axes.ravel().tolist()
    else:
        axes_flat = [axes]  # single subplot case

    # ---- Plot each gene
    for i, g in enumerate(top_genes):
        ax = axes_flat[i]
        # Base boxplot
        sns.boxplot(
            data=combined_df_normalized,
            x=distance_col,
            y=g,
            ax=ax,
            showfliers=False,
            palette=palette,
        )

        # Overlay nonzero points (sampled)
        nonzero = combined_df_normalized.loc[combined_df_normalized[g] > 0, [distance_col, g]]
        if len(nonzero) > 0 and max_points > 0:
            sample_n = min(int(max_points), len(nonzero))
            sample = nonzero.sample(sample_n, random_state=random_state)
            sns.stripplot(
                data=sample,
                x=distance_col,
                y=g,
                ax=ax,
                color="black",
                size=jitter_size,
                alpha=jitter_alpha,
                jitter=jitter,
            )

        ax.set_title(str(g))
        ax.set_xlabel("")
        ax.set_ylabel(y_label)
        ax.tick_params(axis="x", rotation=45)

    # Hide any unused axes
    for j in range(len(top_genes), len(axes_flat)):
        axes_flat[j].axis("off")

    # ---- Figure formatting
    fig.suptitle(suptitle, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if show:
        plt.show()

    return fig, axes_flat, top_genes


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
    anova_df = anova_df.sort_values("anova_pval", na_position="last").reset_index(drop=True)

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


def plot_expression_heatmap(
    summary_df: pd.DataFrame,
    *,
    index_col: str = "gene",
    column_col: str = "distance_bin",
    value_col: str = "mean_expr",
    title: str = "PIG Expression by Distance to Plaque (TgCRND8 17.9m)",
    figsize: tuple[float, float] = (8, 6),
    cmap: str = "mako_r",
    cbar_label: str = "Mean log-expression",
    xlabel: str | None = None,
    ylabel: str | None = None,
    ax: Axes | None = None,
    show: bool = True,
    tight_layout: bool = True,
    cbar_kws: dict[str, Any] | None = None,
    **heatmap_kwargs: Any,
) -> Axes:
    """
    Plot a heatmap of mean expression by distance bin using a pivot of `summary_df`.

    This refactors:
        plt.figure(figsize=(8, 6))
        heatmap_df = summary_df.pivot(index="gene", columns="distance_bin", values="mean_expr")
        sns.heatmap(heatmap_df, cmap="mako_r", cbar_kws={"label": "Mean log-expression"})
        plt.title("PIG Expression by Distance to Plaque (TgCRND8 17.9m)")
        plt.xlabel("Distance bin")
        plt.ylabel("Gene")
        plt.tight_layout()
        plt.show()

    Args:
        summary_df: Input DataFrame containing at least `index_col`, `column_col`, and `value_col`.
        index_col: Column to use for the heatmap's y-axis (rows of the pivot).
        column_col: Column to use for the heatmap's x-axis (columns of the pivot).
        value_col: Column providing cell values in the heatmap.
        title: Figure title.
        figsize: Figure size when creating a new Axes.
        cmap: Colormap for the heatmap.
        cbar_label: Label for the colorbar.
        xlabel: Optional custom x-axis label. Defaults to a title-cased version of `column_col`.
        ylabel: Optional custom y-axis label. Defaults to a title-cased version of `index_col`.
        ax: Optional existing Matplotlib Axes to draw on. If None, a new figure/Axes is created.
        show: If True, calls `plt.show()` at the end.
        tight_layout: If True and a new figure is created, applies `plt.tight_layout()`.
        cbar_kws: Extra kwargs for the colorbar; merged with the label provided in `cbar_label`.
        **heatmap_kwargs: Additional keyword arguments forwarded to `sns.heatmap`.

    Returns:
        The Matplotlib Axes containing the heatmap.

    Raises:
        ValueError: If required columns are missing from `summary_df`.
    """
    required = {index_col, column_col, value_col}
    missing = required - set(summary_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in `summary_df`: {sorted(missing)}")

    # Prepare pivot table for the heatmap
    heatmap_df = summary_df.pivot(index=index_col, columns=column_col, values=value_col)

    created_ax = ax is None
    if created_ax:
        _, ax = plt.subplots(figsize=figsize)

    # Merge colorbar kwargs, ensuring label is set unless explicitly overridden
    merged_cbar_kws: dict[str, Any] = {"label": cbar_label}
    if cbar_kws:
        merged_cbar_kws.update(cbar_kws)

    sns.heatmap(
        heatmap_df,
        ax=ax,
        cmap=cmap,
        cbar_kws=merged_cbar_kws,
        **heatmap_kwargs,
    )

    ax.set_title(title)
    ax.set_xlabel(xlabel if xlabel is not None else column_col.replace("_", " ").title())
    ax.set_ylabel(ylabel if ylabel is not None else index_col.replace("_", " ").title())

    if created_ax and tight_layout:
        plt.tight_layout()
    if show:
        plt.show()

    return ax


def plot_gene_expression_by_distance(
    summary_df: pd.DataFrame,
    pig_genes: Sequence[str],
    *,
    n_cols: int = 4,
    title: str = "Expression of Plaque-Induced Genes vs Distance to Plaque (TgCRND8 17.9m)",
    x_label: str = "Distance bin",
    y_label: str = "Mean log-expression (± SEM)",
    row_height: float = 2.5,
    base_width: float = 14.0,
    sharey: bool = True,
    style: str = "whitegrid",
    distance_col: str = "distance_bin",
    mean_col: str = "mean_expr",
    sem_col: str = "sem_expr",
    gene_col: str = "gene",
    marker: str = "o",
    capsize: float = 3.0,
    linewidth: float = 1.0,
    rotate_xticks: int = 45,
    tight_rect: Sequence[float] = (0.03, 0.03, 1, 0.95),
    show: bool = True,
) -> tuple[Figure, NDArray[Axes]]:
    """
    Plot mean gene expression (± SEM) against distance-to-plaque bins for a set of genes.

    This function filters the provided `pig_genes` to those present in `summary_df`,
    lays out small multiples in a grid, and draws errorbar plots for each gene.
    Global titles and axis labels are added at the figure level.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Long-form dataframe containing at least the columns specified by
        `gene_col`, `distance_col`, `mean_col`, and `sem_col`.
    pig_genes : Sequence[str]
        List or sequence of gene names to include. Only genes present in
        `summary_df[gene_col]` will be plotted.
    n_cols : int, optional
        Number of subplot columns. Defaults to 4.
    title : str, optional
        Figure-level title. Defaults to a descriptive title.
    x_label : str, optional
        Global x-axis label. Defaults to "Distance bin".
    y_label : str, optional
        Global y-axis label. Defaults to "Mean log-expression (± SEM)".
    row_height : float, optional
        Height (inches) of each subplot row. Defaults to 2.5.
    base_width : float, optional
        Figure width in inches. Defaults to 14.0.
    sharey : bool, optional
        Whether to share the y-axis across subplots. Defaults to True.
    style : str, optional
        Seaborn style to apply. Defaults to "whitegrid".
    distance_col : str, optional
        Column name for the distance/bin x-values. Defaults to "distance_bin".
    mean_col : str, optional
        Column name for the mean expression values. Defaults to "mean_expr".
    sem_col : str, optional
        Column name for the SEM values. Defaults to "sem_expr".
    gene_col : str, optional
        Column name for gene identifiers. Defaults to "gene".
    marker : str, optional
        Marker style for errorbar points. Defaults to "o".
    capsize : float, optional
        Capsize for error bars. Defaults to 3.0.
    linewidth : float, optional
        Line width for error bars. Defaults to 1.0.
    rotate_xticks : int, optional
        Rotation (degrees) for x-tick labels in each subplot. Defaults to 45.
    tight_rect : Sequence[float], optional
        Rect parameter for `plt.tight_layout`. Defaults to (0.03, 0.03, 1, 0.95).
    show : bool, optional
        Whether to call `plt.show()` at the end. Defaults to True.

    Returns
    -------
    (Figure, np.ndarray[Axes])
        The Matplotlib figure and a flattened NumPy array of Axes.

    Raises
    ------
    ValueError
        If none of the requested genes are present in `summary_df`.
    """
    # Filter genes to those present in the dataframe
    available_genes = set(summary_df[gene_col].unique())
    genes_to_plot = [g for g in pig_genes if g in available_genes]

    if not genes_to_plot:
        raise ValueError("None of the requested genes were found in the dataframe.")

    n_genes = len(genes_to_plot)
    n_rows = math.ceil(n_genes / n_cols)

    # Styling
    sns.set(style=style)

    # Create subplots
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(base_width, n_rows * row_height),
        sharey=sharey,
    )

    # Normalize axes to a flat array for easy indexing
    axes = np.atleast_1d(axes).flatten()

    # Plot per-gene panels
    for i, g in enumerate(genes_to_plot):
        df_g = summary_df[summary_df[gene_col] == g]
        ax = axes[i]
        ax.errorbar(
            df_g[distance_col],
            df_g[mean_col],
            yerr=df_g[sem_col],
            marker=marker,
            capsize=capsize,
            linewidth=linewidth,
        )
        ax.set_title(g)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=rotate_xticks)

    # Hide any unused axes if grid is larger than number of genes
    for j in range(n_genes, len(axes)):
        axes[j].set_visible(False)

    # Global labels and title
    fig.suptitle(title, fontsize=14)
    fig.text(0.5, 0.04, x_label, ha="center", fontsize=12)
    fig.text(0.04, 0.5, y_label, va="center", rotation="vertical", fontsize=12)

    plt.tight_layout(rect=tight_rect)

    if show:
        plt.show()

    return fig, axes


def summarize_pig_expression_by_distance(
    df: pd.DataFrame,
    pig_genes: Sequence[str],
    distance_col: str = "distance_to_plaque",
    bins: Sequence[float] = (0, 20, 50, 100, 200, np.inf),
    labels: Sequence[str] = (
        "0–20 µm",
        "20–50 µm",
        "50–100 µm",
        "100–200 µm",
        ">200 µm",
    ),
    include_lowest: bool = True,
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

    if len(labels) != (len(bins) - 1):
        raise ValueError(
            f"`labels` must have length len(bins) - 1 = {len(bins) - 1}, got {len(labels)}."
        )

    pig_cols = [g for g in pig_genes if g in df.columns]
    if not pig_cols:
        raise ValueError("None of the specified `pig_genes` are present in the DataFrame.")

    # Work on a copy to avoid mutating the caller's DataFrame
    tmp = df.copy()

    # Assign distance bins
    tmp["distance_bin"] = pd.cut(
        tmp[distance_col],
        bins=bins,
        labels=labels,
        include_lowest=include_lowest,
    )

    # Compute mean and SEM for each distance bin
    mean_expr = tmp.groupby("distance_bin")[pig_cols].mean().reset_index()
    sem_expr = tmp.groupby("distance_bin")[pig_cols].sem().reset_index()

    # Melt to long-form and merge mean/sem
    summary_df = mean_expr.melt(id_vars="distance_bin", var_name="gene", value_name="mean_expr")
    sem_melted = sem_expr.melt(id_vars="distance_bin", var_name="gene", value_name="sem_expr")
    summary_df = summary_df.merge(sem_melted, on=["distance_bin", "gene"], how="left")

    # Cell counts per bin (sorted by categorical order)
    counts = tmp["distance_bin"].value_counts().sort_index()

    return summary_df, counts
