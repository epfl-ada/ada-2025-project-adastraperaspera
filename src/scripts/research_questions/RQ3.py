from __future__ import annotations

from collections.abc import Iterable, Sequence
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import copy
import matplotlib.pyplot as plt
from typing import Iterable, Optional, Tuple


def plot_boxgrid(
    df,
    cols: Iterable[str],
    *,
    figsize_per_plot: Tuple[float, float] = (4.5, 4.0),
    showfliers: bool = True,
    log_scale_cols: Optional[Iterable[str]] = None,
    title: Optional[str] = None,
    annotate_median: bool = True,
    median_fmt: str = "{:.2f}",
    median_fontsize: int = 9,
):
    """
    Plot a 1xN grid of box plots for selected dataframe columns, optionally annotating medians.
    """
    cols = list(cols)
    if not cols:
        raise ValueError("`cols` must contain at least one column name.")

    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing columns in df: {missing}")

    log_set = set(log_scale_cols) if log_scale_cols is not None else set()

    w, h = figsize_per_plot
    fig, axes = plt.subplots(1, len(cols), figsize=(w * len(cols), h), squeeze=False)
    axes = axes[0]

    for ax, col in zip(axes, cols):
        s = df[col].dropna()
        if s.empty:
            ax.set_title(col, fontsize=10)
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes
            )
            ax.set_xticks([])
            ax.grid(axis="y", alpha=0.3)
            continue

        data = s.to_numpy()
        bp = ax.boxplot(data, vert=True, showfliers=showfliers)

        if col in log_set:
            ax.set_yscale("log")

        ax.set_title(col, fontsize=10)
        ax.set_xticks([])
        ax.grid(axis="y", alpha=0.3)

        if annotate_median:
            median_val = float(s.median())
            x = 1

            if ax.get_yscale() == "log":
                y = median_val * 1.08
            else:
                y_min, y_max = ax.get_ylim()
                y = median_val + 0.02 * (y_max - y_min)

            ax.text(
                x,
                y,
                median_fmt.format(median_val),
                ha="center",
                va="bottom",
                fontsize=median_fontsize,
            )

    if title:
        fig.suptitle(title)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()

    plt.show()
    return fig, axes


def compute_pig_neighbor_corr_matrix(
    cells_df: pd.DataFrame,
    pig_genes: list[str],
    use_spearman: bool = False,
    neighbor_prefix: str = "neigh_mean_",
) -> pd.DataFrame:
    """
    Returns a (len(pig_genes) x len(pig_genes)) correlation matrix where:
      rows = target gene expression columns (e.g., 'GeneA')
      cols = neighbor mean feature columns (e.g., 'neigh_mean_GeneB')
    """
    method = "spearman" if use_spearman else "pearson"

    missing_targets = [g for g in pig_genes if g not in cells_df.columns]
    missing_neighbors = [
        f"{neighbor_prefix}{g}"
        for g in pig_genes
        if f"{neighbor_prefix}{g}" not in cells_df.columns
    ]

    if missing_targets or missing_neighbors:
        msg = []
        if missing_targets:
            msg.append(f"Missing target gene columns: {missing_targets}")
        if missing_neighbors:
            msg.append(f"Missing neighbor feature columns: {missing_neighbors}")
        raise KeyError("\n".join(msg))

    corr = pd.DataFrame(index=pig_genes, columns=pig_genes, dtype=float)

    for tg in pig_genes:
        x = cells_df[tg]
        for ng in pig_genes:
            y = cells_df[f"{neighbor_prefix}{ng}"]
            corr.loc[tg, ng] = x.corr(y, method=method)

    return corr


def plot_corr_matrix(
    corr_df: pd.DataFrame,
    title: str = "PIG target vs neighbor-mean correlation",
    mask_upper_triangle: bool = True,
    cmap: str = "RdBu_r",
    mask_color: str = "black",
    color_scale: str = "zscore",
    z_clip: float | None = 2.5,
    raw_vmin: float | None = None,
    raw_vmax: float | None = None,
    annotate: bool = True,
    fmt: str = ".2f",
    fontsize: int = 9,
    mask_diagonal: bool = True,
):
    """
    Matplotlib correlation matrix plot with:
      - per-cell annotations
      - optional upper-triangle masking (blacked out)
      - optional diagonal masking (blacked out), independent of upper-triangle masking
      - diverging colormap centered at 0
      - optional z-score scaling for colors (annotations remain raw correlations)
    """
    data = corr_df.to_numpy(dtype=float)
    nrows, ncols = data.shape

    structural_mask = np.zeros_like(data, dtype=bool)

    if mask_upper_triangle:
        structural_mask |= np.triu(np.ones_like(data, dtype=bool), k=1)

    if mask_diagonal:
        structural_mask |= np.eye(nrows, ncols, dtype=bool)

    invalid_mask = ~np.isfinite(data)
    full_mask = structural_mask | invalid_mask

    if color_scale.lower() == "zscore":
        vals = data[~full_mask]
        if vals.size == 0:
            z = np.zeros_like(data)
            sigma = 1.0
        else:
            mu = vals.mean()
            sigma = vals.std(ddof=0)
            z = np.zeros_like(data) if sigma == 0 else (data - mu) / sigma

        plot_arr = np.ma.array(z, mask=full_mask)

        if z_clip is not None:
            plot_arr = np.ma.clip(plot_arr, -z_clip, z_clip)
            norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-z_clip, vmax=z_clip)
            cbar_label = f"Correlation (z-score, clipped to ±{z_clip})"
        else:
            maxabs = (
                float(np.nanmax(np.abs(z[~full_mask]))) if np.any(~full_mask) else 1.0
            )
            maxabs = max(maxabs, 1e-12)
            norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-maxabs, vmax=maxabs)
            cbar_label = "Correlation (z-score)"
    else:
        plot_arr = np.ma.array(data, mask=full_mask)

        if raw_vmin is not None and raw_vmax is not None:
            norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=raw_vmin, vmax=raw_vmax)
        else:
            vals = data[~full_mask]
            maxabs = float(np.nanmax(np.abs(vals))) if vals.size else 1.0
            maxabs = max(maxabs, 1e-12)
            norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-maxabs, vmax=maxabs)

        cbar_label = "Correlation"

    cmap_obj = plt.get_cmap(cmap)
    try:
        cmap_obj = cmap_obj.copy()
    except AttributeError:
        cmap_obj = copy.copy(cmap_obj)
    cmap_obj.set_bad(color=mask_color)

    figsize = (max(8, 0.55 * ncols + 4), max(6, 0.55 * nrows + 3))
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(plot_arr, cmap=cmap_obj, norm=norm)

    ax.set_title(title)
    ax.set_xticks(np.arange(ncols))
    ax.set_yticks(np.arange(nrows))
    ax.set_xticklabels(corr_df.columns, rotation=45, ha="right")
    ax.set_yticklabels(corr_df.index)

    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which="minor", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    if annotate:
        for i in range(nrows):
            for j in range(ncols):
                if full_mask[i, j]:
                    continue
                val = data[i, j]
                if not np.isfinite(val):
                    continue

                bg = plot_arr[i, j]
                bg_norm = (
                    0.5
                    if (bg is np.ma.masked or not np.isfinite(float(bg)))
                    else float(norm(float(bg)))
                )
                text_color = "white" if (bg_norm < 0.25 or bg_norm > 0.75) else "black"

                ax.text(
                    j,
                    i,
                    format(val, fmt),
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    color=text_color,
                )

    fig.tight_layout()
    plt.show()


def analyze_plaque_distance_effects(
    cells: pd.DataFrame,
    pig_genes: Iterable[str],
    *,
    apoe_gene: str = "Apoe",
    baseline_broad_type: str | None = "Unlabeled",
    predicted_label_col: str = "predicted_label",
    unknown_label_value: int = 34,
    cell_type_col: str = "cell_type",
    distance_col: str = "distance_to_plaque",
    distance_bin_col: str = "distance_bin",
    broad_type_col: str = "broad_type",
    create_distance_bins: bool = True,
    logger: logging.Logger | None = None,
) -> tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run OLS models relating gene expression to distance from plaque while controlling
    for broad cell type (with a configurable baseline level), and compute binned
    summary stats for PIG genes.

    Steps performed:
      1) If `predicted_label_col` exists, mark cells with `unknown_label_value` as unlabeled
         by setting their `cell_type_col` to NA.
      2) Ensure `broad_type_col` exists by simplifying `cell_type_col` into:
         {"Glutamatergic", "GABAergic", "Glial", "Unlabeled", "Other"}.
      3) Fit an OLS for `apoe_gene` (treatment coding for broad_type with chosen baseline):
             Q(apoe_gene) ~ distance + C(broad_type, Treatment(reference=...))
         and log the model summary.
      4) For each gene in `pig_genes` present in the data, fit:
             Q(gene) ~ distance + C(broad_type, Treatment(reference=...))
         and collect coefficients and p-values into a tidy DataFrame.
      5) Produce mean, SEM, and counts per (gene, broad_type, distance_bin). If
         `distance_bin_col` is missing and `create_distance_bins=True`, create it from
         `distance_col` using bins: [0,20), [20,50), [50,100), [100,200), [200, ∞).

    Parameters
    ----------
    cells : pd.DataFrame
        DataFrame with at least the columns for distance (`distance_col`) and cell type
        (`cell_type_col`). Optionally contains `predicted_label_col` and `distance_bin_col`.
    pig_genes : Iterable[str]
        Collection of PIG gene names to analyze.
    apoe_gene : str, optional
        Name of the Apoe column (default "Apoe").
    baseline_broad_type : str | None, optional
        Desired baseline level for treatment coding of `broad_type`.
        Defaults to "Unlabeled". If not present, falls back to a sensible
        available level (tries "Other" next, then the most frequent level).
    predicted_label_col : str, optional
        Column indicating predicted label class (default "predicted_label").
    unknown_label_value : int, optional
        Value in `predicted_label_col` that means "unknown" (default 34).
    cell_type_col : str, optional
        Column with detailed cell types to be simplified (default "cell_type").
    distance_col : str, optional
        Column with distance to plaque (default "distance_to_plaque").
    distance_bin_col : str, optional
        Column with categorical distance bins (default "distance_bin").
    broad_type_col : str, optional
        Output column for simplified cell types (default "broad_type").
    create_distance_bins : bool, optional
        If True and `distance_bin_col` is missing, compute bins from `distance_col`.
    logger : logging.Logger, optional
        Logger to use for info messages. Defaults to a module-level logger.

    Returns
    -------
    apoe_summary_text : str
        Text representation of the OLS summary for the Apoe model.
    gene_coefs : pd.DataFrame
        Tidy coefficients for each gene with columns:
        ["gene", "variable", "coef", "pval"].
    binned_stats : pd.DataFrame
        Aggregated statistics with columns:
        ["gene", "broad_type", "distance_bin", "mean", "sem", "n"].
    cells_out : pd.DataFrame
        A copy of the input DataFrame with ensured `broad_type_col`
        (and `distance_bin_col` if created).

    Raises
    ------
    ValueError
        If required columns are missing.
    """
    logger = logger or logging.getLogger(__name__)
    pig_genes = list(pig_genes)

    required_cols = {cell_type_col, distance_col}
    missing = required_cols - set(cells.columns)
    if missing:
        raise ValueError(f"Missing required columns in `cells`: {sorted(missing)}")

    cells = cells.copy()

    if predicted_label_col in cells.columns:
        mask_unknown = cells[predicted_label_col] == unknown_label_value
        n_unknown = int(mask_unknown.sum())
        if n_unknown:
            logger.info(
                "Setting %d rows with %s == %s to NA in %s.",
                n_unknown,
                predicted_label_col,
                unknown_label_value,
                cell_type_col,
            )
            cells.loc[mask_unknown, cell_type_col] = pd.NA

    def _simplify_celltype(ct: object) -> str:
        if pd.isna(ct):
            return "Unlabeled"
        celltype_to_broad = {
            "Corticothalamic, Gluta": "Neuron_Glutamatergic",
            "Dentate, Gluta": "Neuron_Glutamatergic",
            "Hypothalamic Gnrh1, Gluta": "Neuron_Glutamatergic",
            "Hypothalamic medial, Gluta": "Neuron_Glutamatergic",
            "Intra/Extratelencephalic, Gluta": "Neuron_Glutamatergic",
            "Olfactory bulb, Gluta": "Neuron_Glutamatergic",
            "Pineal, Gluta": "Neuron_Glutamatergic",
            "Pons, Gluta": "Neuron_Glutamatergic",
            "Thalamic, Gluta": "Neuron_Glutamatergic",
            "Cerebellar, GABA": "Neuron_GABAergic",
            "Cerebral LGE, GABA": "Neuron_GABAergic",
            "Cortex caudal, GABA": "Neuron_GABAergic",
            "Cortex medial, GABA": "Neuron_GABAergic",
            "Hypothalamic GABA": "Neuron_GABAergic",
            "Medulla, GABA": "Neuron_GABAergic",
            "Astrocyte": "Glia_Astrocyte_Ependymal",
            "Oligodendrocyte": "Glia_Oligodendrocyte_Lineage",
            "Immune": "Immune_Microglia_Macrophage",
            "Vascular": "Vascular_Endothelial_Pericyte",
        }
        return celltype_to_broad.get(ct, "Unlabeled")

    if broad_type_col not in cells.columns:
        logger.info("Creating '%s' by simplifying '%s'.", broad_type_col, cell_type_col)
        cells[broad_type_col] = cells[cell_type_col].apply(_simplify_celltype)

    if distance_bin_col not in cells.columns and create_distance_bins:
        logger.info(
            "Creating '%s' from '%s' using 5 equal quantile bins.",
            distance_bin_col,
            distance_col,
        )
        cells[distance_bin_col] = pd.qcut(
            cells[distance_col].astype(float), q=5, duplicates="drop"
        )

    observed_levels = set(cells[broad_type_col].dropna().astype(str).unique())
    ref = baseline_broad_type
    if ref is None or ref not in observed_levels:
        for candidate in ("Unlabeled", "Other"):
            if candidate in observed_levels:
                ref = candidate
                break
        else:

            ref = cells[broad_type_col].dropna().astype(str).value_counts().idxmax()

    if apoe_gene not in cells.columns:
        raise ValueError(
            f"Column '{apoe_gene}' not found in `cells` for the Apoe model."
        )
    apoe_df = cells[[apoe_gene, distance_col, broad_type_col]].dropna()
    if apoe_df.empty:
        raise ValueError("No rows available for the Apoe model after dropping NAs.")
    apoe_formula = f"Q('{apoe_gene}') ~ {distance_col} + C({broad_type_col}, Treatment(reference='{ref}'))"
    apoe_res = smf.ols(apoe_formula, data=apoe_df).fit()
    apoe_summary_text = apoe_res.summary().as_text()
    logger.info("\n%s", apoe_summary_text)

    results: list[dict] = []
    present_genes = [g for g in pig_genes if g in cells.columns]
    missing_genes = [g for g in pig_genes if g not in cells.columns]
    if missing_genes:
        logger.info(
            "Skipping %d genes not found in data: %s",
            len(missing_genes),
            ", ".join(missing_genes[:10]) + ("..." if len(missing_genes) > 10 else ""),
        )

    for gene in present_genes:
        df_gene = cells[[distance_col, gene, broad_type_col]].dropna()
        if df_gene.empty:
            continue
        formula = f"Q('{gene}') ~ {distance_col} + C({broad_type_col}, Treatment(reference='{ref}'))"
        res = smf.ols(formula, data=df_gene).fit()
        for var in res.params.index:
            results.append(
                {
                    "gene": gene,
                    "variable": var,
                    "coef": float(res.params[var]),
                    "pval": float(res.pvalues[var]),
                }
            )
    gene_coefs = pd.DataFrame(results)

    if gene_coefs.empty:
        logger.info(
            "No coefficients produced (no qualifying PIG gene models could be fit)."
        )
    else:
        logger.info(
            "Fitted %d gene models. Example rows:\n%s",
            gene_coefs["gene"].nunique(),
            gene_coefs.head().to_string(index=False),
        )

    if distance_bin_col not in cells.columns:

        logger.info(
            "Column '%s' not available; binned statistics will be empty.",
            distance_bin_col,
        )
        binned_stats = pd.DataFrame(
            columns=["gene", "broad_type", "distance_bin", "mean", "sem", "n"]
        )
    else:
        long = cells.melt(
            id_vars=[distance_bin_col, broad_type_col],
            value_vars=present_genes,
            var_name="gene",
            value_name="expression",
        ).dropna(subset=["expression"])

        def _sem(x: pd.Series) -> float:
            n = len(x)
            if n <= 1:
                return 0.0
            return float(x.std(ddof=1) / np.sqrt(n))

        agg = (
            long.groupby(["gene", broad_type_col, distance_bin_col], observed=True)
            .agg(
                mean=("expression", "mean"),
                sem=("expression", _sem),
                n=("expression", "size"),
            )
            .reset_index()
        )

        if hasattr(agg[distance_bin_col].dtype, "ordered") and getattr(
            agg[distance_bin_col].dtype, "ordered", False
        ):
            agg = agg.sort_values(["gene", broad_type_col, distance_bin_col])

        binned_stats = agg.rename(
            columns={broad_type_col: "broad_type", distance_bin_col: "distance_bin"}
        )

    if distance_bin_col in cells.columns and hasattr(
        cells[distance_bin_col].dtype, "categories"
    ):
        bin_order = list(map(str, cells[distance_bin_col].cat.categories))
    else:
        bin_order = (
            list(map(str, sorted(cells[distance_bin_col].unique())))
            if distance_bin_col in cells.columns
            else []
        )

    return apoe_summary_text, gene_coefs, binned_stats, cells, agg, bin_order


def _clean_cell_id(x: Any) -> str:
    """
    Normalize a cell identifier to a clean UTF-8 string.

    Handles bytes-like inputs and strings that look like ``"b'xxxx'"``.

    Parameters
    ----------
    x : Any
        Original cell identifier (bytes, bytearray, numpy bytes, or any object).

    Returns
    -------
    str
        A normalized string representation of the cell ID.
    """
    if isinstance(x, (bytes, bytearray, np.bytes_)):
        return x.decode("utf-8", errors="ignore")
    s = str(x)
    if s.startswith("b'") and s.endswith("'"):
        return s[2:-1]
    return s


def process_cell_annotations(
    combined_df_normalized: pd.DataFrame,
    annotation_csv: str | Path,
    *,
    distance_col: str = "distance_to_plaque",
    distance_bins: Sequence[float] = (0, 20, 50, 100, 200, 1e9),
    distance_labels: Sequence[str] = (
        "0-20",
        "20-50",
        "50-100",
        "100-200",
        ">200",
    ),
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    """
    Merge cell annotations into a normalized dataframe, build distance bins,
    and compute per-bin cell type proportions.

    This function refactors and generalizes the provided script into a reusable
    pipeline. It reads the annotation CSV, infers/renames the cell ID column,
    normalizes the ID format, merges with the provided dataframe, constructs
    distance bins (if missing), and returns diagnostic metrics along with the
    merged data and proportions.

    Parameters
    ----------
    combined_df_normalized : pandas.DataFrame
        Input dataframe containing at least a ``"cell_id"`` column (or a column
        that will be merged to ``"cell_id"`` from the annotation file), and
        typically a distance column (default: ``"distance_to_plaque"``).
    annotation_csv : str | pathlib.Path
        Path to the annotation CSV file. Must contain a cell identifier column
        (e.g., ``cell_id`` or similar containing both "cell" and "id" in its name).
    distance_col : str, optional
        Name of the column with distances used to create bins if ``"distance_bin"``
        is not already present, by default ``"distance_to_plaque"``.
    distance_bins : Sequence[float], optional
        Bin edges passed to ``pandas.cut`` for distance binning, by default
        ``(0, 20, 50, 100, 200, 1e9)``.
    distance_labels : Sequence[str], optional
        Labels corresponding to ``distance_bins`` intervals, by default
        ``("0-20", "20-50", "50-100", "100-200", ">200")``.
    logger : logging.Logger, optional
        Logger to use for info messages. If ``None``, uses ``logging.getLogger(__name__)``.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing:
        - ``"cells"`` : pandas.DataFrame
            The merged dataframe, including (optional) ``"distance_bin"`` and annotations.
        - ``"props"`` : pandas.DataFrame
            A dataframe with columns ``["distance_bin", "cell_type", "n", "proportion"]``.
        - ``"match_rate"`` : float
            The fraction of rows with a non-null ``"cell_type"`` after merge.
        - ``"unknown_count"`` : int
            Count of cells with ``predicted_label == 34`` (if present).
        - ``"merged_rows_before"`` : int
            Row count of the input dataframe prior to merging.
        - ``"merged_rows_after"`` : int
            Row count of the merged dataframe.
        - ``"bin_proportion_sums"`` : pandas.Series
            Sanity check: sum of proportions per ``distance_bin`` (should be ~1.0).

    Raises
    ------
    AssertionError
        If no cell ID-like column can be identified in the annotation CSV.
    FileNotFoundError
        If ``annotation_csv`` does not exist.
    ValueError
        If distance binning is needed but ``distance_col`` is missing from the data.
    """
    logger = logger or logging.getLogger(__name__)

    annotation_csv = Path(annotation_csv)
    if not annotation_csv.exists():
        raise FileNotFoundError(f"Annotation CSV not found: {annotation_csv}")

    annot = pd.read_csv(annotation_csv)

    id_col = next(
        (c for c in annot.columns if "cell" in c.lower() and "id" in c.lower()),
        None,
    )
    assert (
        id_col is not None
    ), f"No 'cell_id' column detected in annotations. Available columns: {annot.columns.tolist()}"

    annot = annot.rename(columns={id_col: "cell_id"})
    annot["cell_id"] = annot["cell_id"].apply(_clean_cell_id)

    keep_cols = ["cell_id", "cell_type", "coord_X", "coord_Y", "predicted_label"]
    keep_cols = [c for c in keep_cols if c in annot.columns]
    annot = annot[keep_cols]

    logger.info(annot.head())

    before = combined_df_normalized.shape[0]
    cells = combined_df_normalized.copy()
    cells["cell_id"] = (
        cells["cell_id"].apply(_clean_cell_id)
        if "cell_id" in cells.columns
        else cells["cell_id"]
    )
    cells = cells.merge(annot, on="cell_id", how="left")
    logger.info("Merged: %d -> %d rows", before, cells.shape[0])

    match_rate = (
        cells["cell_type"].notna().mean() if "cell_type" in cells.columns else 0.0
    )
    logger.info("Annotated cells rate: %.1f%%", 100.0 * match_rate)

    unknown_count = 0
    if "predicted_label" in cells.columns:
        unknown_count = int((cells["predicted_label"] == 34).sum())
        if unknown_count:
            logger.info(
                "Attention: %d cells with predicted_label==34 (not a valid tag).",
                unknown_count,
            )

    if "distance_bin" not in cells.columns:
        if distance_col not in cells.columns:
            raise ValueError(
                f"'distance_bin' is missing and '{distance_col}' not found to create it."
            )
        cells["distance_bin"] = pd.qcut(
            cells[distance_col].astype(float), q=5, duplicates="drop"
        )

    if hasattr(cells["distance_bin"].dtype, "ordered"):
        cells["distance_bin"] = cells["distance_bin"].cat.as_ordered()

    props = (
        cells.groupby(["distance_bin", "cell_type"]).size().rename("n").reset_index()
    )
    props["proportion"] = props.groupby("distance_bin")["n"].transform(
        lambda x: x / x.sum()
    )

    bin_proportion_sums = props.groupby("distance_bin")["proportion"].sum().round(6)
    logger.info("Sum of proportions per bin: %s", dict(bin_proportion_sums))

    return {
        "cells": cells,
        "props": props,
        "match_rate": float(match_rate),
        "unknown_count": int(unknown_count),
        "merged_rows_before": int(before),
        "merged_rows_after": int(cells.shape[0]),
        "bin_proportion_sums": bin_proportion_sums,
    }
