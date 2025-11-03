from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import logging
from typing import Any

from anndata import AnnData
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.formula.api as smf


def analyze_plaque_distance(
    df: pd.DataFrame,
    *,
    column: str = "distance_to_plaque",
    prox_thresh: float = 30.0,
    distal_thresh: float = 100.0,
    n_bins: int = 60,
    kde: bool = True,
    figsize: tuple[float, float] = (6.0, 4.0),
    category_col: str = "dist_bin_simple",
    palette: str = "viridis",
    show_plots: bool = True,
    in_place: bool = True,
    logger: logging.Logger | None = None,
) -> dict[str, Any] | tuple[dict[str, Any], pd.DataFrame]:
    """
    Summarize and visualize distances to the nearest plaque, and classify cells
    into proximity bins.

    This function replicates and generalizes the original analysis:
      1) Logs basic summary statistics for the distance column.
      2) Plots a histogram on linear scale.
      3) Plots a histogram of log10(distance + 1).
      4) Computes key quantiles (5, 25, 50, 75, 95%).
      5) Categorizes cells into proximal / intermediate / distal and plots counts.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing the distance column.
    column : str, optional
        Name of the column with distances (in µm). Default is "distance_to_plaque".
    prox_thresh : float, optional
        Threshold (in µm) for the proximal class upper bound. Default is 30.0.
    distal_thresh : float, optional
        Threshold (in µm) for the intermediate class upper bound
        (distal starts above this). Default is 100.0.
    n_bins : int, optional
        Number of bins for histograms. Default is 60.
    kde : bool, optional
        Whether to overlay a KDE on histograms. Default is True.
    figsize : tuple[float, float], optional
        Figure size for each plot. Default is (6.0, 4.0).
    category_col : str, optional
        Name of the output categorical column to add. Default is "dist_bin_simple".
    palette : str, optional
        Seaborn palette for the countplot. Default is "viridis".
    show_plots : bool, optional
        If True, shows the generated plots. Default is True.
    in_place : bool, optional
        If True, the new category column is added to `df`. If False, a copy is made
        and returned alongside the results. Default is True.
    logger : logging.Logger | None, optional
        Logger to use for informational output. If None, a module-level logger is used.

    Returns
    -------
    results : dict
        A dictionary with keys:
            - "describe": pd.Series of basic summary statistics.
            - "quantiles": dict with keys "5%", "25%", "50% (median)", "75%", "95%".
            - "proximity_counts": pd.Series of counts per proximity class (ordered).
            - "labels": list of class labels used.
            - "category_col": name of the column added.
    (results, df_out) : tuple
        If `in_place` is False, returns a tuple where `df_out` is the modified copy.

    Raises
    ------
    ValueError
        If `column` is not present in `df`, or if thresholds are invalid.
    """
    _logger = logger or logging.getLogger(__name__)

    # --- Validation ---
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")
    if prox_thresh <= 0 or distal_thresh <= prox_thresh:
        raise ValueError("Thresholds must satisfy: 0 < prox_thresh < distal_thresh.")

    # Choose working DataFrame
    df_out = df if in_place else df.copy()

    # Extract data vector, drop NaNs
    x = df_out[column].astype(float)
    x_clean = x.dropna()

    # --- 1. Basic summary statistics ---
    _logger.info("=== Distance to nearest plaque (µm) summary ===")
    desc = x.describe()
    _logger.info(desc)

    # --- 2. Histogram (linear scale) ---
    if show_plots:
        plt.figure(figsize=figsize)
        sns.histplot(x_clean, bins=n_bins, kde=kde)
        plt.xlabel("Distance to nearest plaque (µm)")
        plt.ylabel("Cell count")
        plt.title("Distribution of cell distances to nearest plaque (linear scale)")
        plt.tight_layout()
        plt.show()

    # --- 3. Histogram (log scale) ---
    if show_plots:
        plt.figure(figsize=figsize)
        # Use +1 to remain defined at zero
        sns.histplot(np.log10(x_clean + 1.0), bins=n_bins, kde=kde)
        plt.xlabel("log10(Distance + 1)")
        plt.ylabel("Cell count")
        plt.title("Distribution of distances (log10 scale)")
        plt.tight_layout()
        plt.show()

    # --- 4. Quantiles and thresholds ---
    q5, q25, q50, q75, q95 = np.percentile(x_clean, [5, 25, 50, 75, 95])
    quantiles: dict[str, float] = {
        "5%": float(q5),
        "25%": float(q25),
        "50% (median)": float(q50),
        "75%": float(q75),
        "95%": float(q95),
    }
    _logger.info(quantiles)

    # --- 5. Proximity classes ---
    bins = [0.0, float(prox_thresh), float(distal_thresh), np.inf]
    labels: list[str] = [
        f"proximal (<{prox_thresh:g} µm)",
        f"intermediate ({prox_thresh:g}–{distal_thresh:g} µm)",
        f"distal (>{distal_thresh:g} µm)",
    ]
    df_out[category_col] = pd.cut(x, bins=bins, labels=labels, include_lowest=True, right=False)

    _logger.info("\n=== Cell counts by proximity class ===")
    counts = df_out[category_col].value_counts().reindex(labels, fill_value=0)
    _logger.info(counts)

    if show_plots:
        sns.countplot(y=category_col, data=df_out, order=labels, palette=palette)
        plt.title("Cell counts by plaque proximity class")
        plt.xlabel("Count")
        plt.ylabel("Proximity class")
        plt.tight_layout()
        plt.show()

    results: dict[str, Any] = {
        "describe": desc,
        "quantiles": quantiles,
        "proximity_counts": counts,
        "labels": labels,
        "category_col": category_col,
    }

    return results if in_place else (results, df_out)


@dataclass
class LeidenAnalysisResult:
    """
    Container for outputs of `analyze_leiden_spatial`.
    """

    df: pd.DataFrame  # Copy of the input with added columns (cluster, bins, etc.)
    adata: AnnData  # AnnData object used by Scanpy
    gene_cols: list[str]  # Detected numeric gene-expression columns
    freq_df: pd.DataFrame  # Cluster frequency by distance bin (with midpoints & percentages)
    logit_df: (
        pd.DataFrame
    )  # Logistic regression summary per cluster (with Bonferroni adj. p-values)
    expr_z: pd.DataFrame  # Z-scored marker expression (clusters x genes)
    top_z: pd.Series  # For each marker gene, the cluster with the highest z-score


def analyze_leiden_spatial(
    combined_df: pd.DataFrame,
    *,
    meta_columns: Iterable[str] | None = None,
    leiden_resolution: float = 0.6,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    bin_width_um: float = 25.0,
    proximal_threshold_um: float = 30.0,
    sample_for_scatter: int | None = 20_000,
    random_state: int = 0,
    marker_genes: Mapping[str, Sequence[str]] | None = None,
    plot: bool = True,
    logger: logging.Logger | None = None,
) -> LeidenAnalysisResult:
    """
    Run a Leiden clustering + spatial analysis pipeline on a single-cell dataframe.

    This function:
      1) Detects numeric gene-expression columns (excluding metadata).
      2) Builds an AnnData object and runs Scanpy: scale → PCA → neighbors → Leiden → UMAP.
      3) Adds the Leiden cluster labels back to the dataframe.
      4) Visualizes UMAP and spatial distribution (optional).
      5) Bins distances to plaque and computes cluster frequency trends vs. distance.
      6) Runs per-cluster logistic regression of membership vs. distance.
      7) Computes z-scored marker enrichment per cluster and reports top clusters per marker gene.

    Parameters
    ----------
    combined_df : pd.DataFrame
        Input dataframe with at least the following columns:
        'cell_id', 'distance_to_plaque', and (for spatial scatter) 'x_centroid', 'y_centroid'.
        Gene-expression columns should be numeric.
    meta_columns : Iterable[str], optional
        Column names to exclude from gene detection. If None, a sensible default set is used.
    leiden_resolution : float, default=0.6
        Resolution parameter for `sc.tl.leiden`.
    n_neighbors : int, default=15
        Number of neighbors for `sc.pp.neighbors`.
    n_pcs : int, default=30
        Number of principal components to use for neighbors/UMAP.
    bin_width_um : float, default=25.0
        Width of distance bins (in µm) for frequency curves.
    proximal_threshold_um : float, default=30.0
        Threshold (in µm) for creating a binary 'is_proximal' label.
    sample_for_scatter : int or None, default=20000
        If set, subsample this many rows for the spatial scatterplot (speeds up plotting).
        If None, use all rows.
    random_state : int, default=0
        Random seed for subsampling.
    marker_genes : Mapping[str, Sequence[str]], optional
        Dict of marker lists by cell type to evaluate enrichment. If None, a default set is used.
    plot : bool, default=True
        If True, produce UMAP and spatial plots.
    logger : logging.Logger, optional
        Logger for progress messages. If None, a module logger is used.

    Returns
    -------
    LeidenAnalysisResult
        Dataclass holding the augmented dataframe, AnnData, gene columns, frequency table,
        logistic regression table, z-scored marker matrix, and per-gene top cluster.
    """
    log = logger or logging.getLogger(__name__)
    df = combined_df.copy()

    # -----------------------------
    # 0) Defaults & validations
    # -----------------------------
    required_cols = {"cell_id", "distance_to_plaque"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input dataframe is missing required columns: {missing}")

    default_meta = {
        "cell_id",
        "x_centroid",
        "y_centroid",
        "transcript_counts",
        "control_probe_counts",
        "control_codeword_counts",
        "unassigned_codeword_counts",
        "total_counts",
        "cell_area",
        "nucleus_area",
        "distance_to_plaque",
        "nearest_plaque_center_dist",
        "inside_any_plaque",
        "nearest_plaque_id",
        "nearest_plaque_area",
        "distance_bin",
        "dist_bin_simple",
    }
    meta_cols = set(meta_columns) if meta_columns is not None else default_meta

    # -----------------------------
    # 1) Detect numeric gene columns
    # -----------------------------
    gene_cols: list[str] = [
        c for c in df.columns if c not in meta_cols and pd.api.types.is_numeric_dtype(df[c])
    ]
    log.info("Detected %d gene columns", len(gene_cols))

    if len(gene_cols) == 0:
        raise ValueError("No numeric gene columns detected. Check `meta_columns` or input dtypes.")

    # -----------------------------
    # 2) Build AnnData for Scanpy
    # -----------------------------
    adata = sc.AnnData(df[gene_cols].to_numpy())
    # Minimal obs; 'dist_bin_simple' optional
    adata.obs = pd.DataFrame(
        {
            "cell_id": df["cell_id"].astype(str).values,
            "distance_to_plaque": pd.to_numeric(df["distance_to_plaque"], errors="coerce").values,
            "dist_bin_simple": (
                df["dist_bin_simple"].astype(str).values
                if "dist_bin_simple" in df.columns
                else pd.Series([""] * len(df)).values
            ),
        },
        index=df["cell_id"].astype(str).values,
    )
    adata.obs_names = df["cell_id"].astype(str).values
    adata.var_names = pd.Index(gene_cols, name="genes")

    log.info("AnnData created with shape %s", adata.shape)
    log.info("=== PCA → Neighbors → Leiden Clustering ===")

    # -----------------------------
    # 3) Scanpy pipeline
    # -----------------------------
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden")

    log.info("Leiden cluster counts:\n%s", adata.obs["leiden"].value_counts())

    # UMAP (for visualization)
    sc.tl.umap(adata)
    if plot:
        sc.pl.umap(
            adata,
            color=["leiden"],
            legend_loc="on data",
            frameon=False,
            title="Leiden clusters",
            show=True,
        )

    # -----------------------------
    # 4) Add cluster labels to df
    # -----------------------------
    df["cluster_leiden"] = adata.obs["leiden"].reindex(df["cell_id"].astype(str)).values
    log.info("Added 'cluster_leiden' column to dataframe.")

    # -----------------------------
    # 5) Spatial distribution plot
    # -----------------------------
    if plot and {"x_centroid", "y_centroid"}.issubset(df.columns):
        plot_df = df
        if sample_for_scatter is not None and sample_for_scatter < len(df):
            plot_df = df.sample(sample_for_scatter, random_state=random_state)

        plt.figure(figsize=(8, 8))
        sns.scatterplot(
            data=plot_df,
            x="x_centroid",
            y="y_centroid",
            hue="cluster_leiden",
            palette="tab20",
            s=6,
            linewidth=0,
            alpha=0.7,
        )
        plt.gca().invert_yaxis()  # match microscopy orientation
        plt.title("Spatial map of Leiden clusters")
        plt.xlabel("X coordinate (µm)")
        plt.ylabel("Y coordinate (µm)")
        plt.legend(bbox_to_anchor=(1.05, 1), title="Cluster", ncol=1)
        plt.tight_layout()
        plt.show()
    elif plot:
        log.warning("Skipping spatial scatterplot: 'x_centroid' and/or 'y_centroid' not found.")

    if plot:
        sc.pl.umap(
            adata,
            color=["distance_to_plaque"],
            color_map="viridis",
            size=15,
            alpha=0.8,
            frameon=False,
            title="UMAP colored by distance to plaque",
            show=True,
        )

    # -----------------------------------------------
    # 6) Cluster frequency vs. continuous distance
    # -----------------------------------------------
    # Smooth binning up to ~99th percentile to reduce outliers
    max_dist = np.nanpercentile(pd.to_numeric(df["distance_to_plaque"], errors="coerce"), 99)
    bins = np.arange(0, max_dist + bin_width_um, bin_width_um)
    df["distance_bin_cont"] = pd.cut(
        pd.to_numeric(df["distance_to_plaque"], errors="coerce"), bins=bins
    )

    freq_df = df.groupby(["distance_bin_cont", "cluster_leiden"]).size().reset_index(name="n_cells")
    totals = (
        df["distance_bin_cont"]
        .value_counts()
        .rename_axis("distance_bin_cont")
        .reset_index(name="total_cells")
    )
    freq_df = freq_df.merge(totals, on="distance_bin_cont", how="left")
    freq_df["pct"] = 100 * freq_df["n_cells"] / freq_df["total_cells"]
    freq_df["bin_mid"] = freq_df["distance_bin_cont"].apply(
        lambda x: x.mid if pd.notna(x) else np.nan
    )

    if plot and not freq_df.empty:
        plt.figure(figsize=(10, 6))
        sns.lineplot(
            data=freq_df,
            x="bin_mid",
            y="pct",
            hue="cluster_leiden",
            marker="o",
            alpha=0.7,
        )
        plt.title("Cluster frequency (%) vs. distance to nearest plaque")
        plt.xlabel(f"Distance to plaque (µm, binned every {int(bin_width_um)} µm)")
        plt.ylabel("% of cells in each bin")
        plt.legend(bbox_to_anchor=(1.05, 1), title="Cluster", ncol=1)
        plt.tight_layout()
        plt.show()

    # -----------------------------------------------
    # 7) Logistic regression: in_cluster ~ distance
    # -----------------------------------------------
    df["is_proximal"] = (
        pd.to_numeric(df["distance_to_plaque"], errors="coerce") <= proximal_threshold_um
    ).astype(int)

    results: list[dict] = []
    for cl in sorted(df["cluster_leiden"].dropna().unique(), key=str):
        df_cl = df.copy()
        df_cl["in_cluster"] = (df_cl["cluster_leiden"] == cl).astype(int)

        try:
            model = smf.logit("in_cluster ~ distance_to_plaque", data=df_cl).fit(disp=0)
            slope = float(model.params.get("distance_to_plaque", np.nan))
            pval = float(model.pvalues.get("distance_to_plaque", np.nan))
        except Exception as e:
            log.warning("Logit failed for cluster %s (%s). Setting slope/pval = NaN.", cl, e)
            slope, pval = np.nan, np.nan

        results.append({"cluster": str(cl), "slope": slope, "pval": pval})

    logit_df = pd.DataFrame(results)
    if not logit_df.empty and "pval" in logit_df:
        # Bonferroni
        m = logit_df["pval"].notna().sum()
        logit_df["adj_pval"] = np.minimum(logit_df["pval"] * max(m, 1), 1.0)
    else:
        logit_df["adj_pval"] = np.nan

    # -----------------------------------------------
    # 8) Marker enrichment heatmap (z-scored)
    # -----------------------------------------------
    if marker_genes is None:
        marker_genes = {
            "microglia": ["Apoe", "Hexb", "Cst3", "C4b"],
            "astrocyte": ["Gfap", "Serpina3n", "Vim", "S100a6"],
            "neuron": ["Nrep", "Snap25", "Rbfox3"],
            "oligodendrocyte": ["Mbp", "Plp1", "Cldn11"],
        }

    genes_to_check = [g for genes in marker_genes.values() for g in genes if g in df.columns]
    expr_z = pd.DataFrame()
    top_z = pd.Series(dtype=object)

    if genes_to_check:
        expr = df.groupby("cluster_leiden")[genes_to_check].mean(numeric_only=True)
        # z-score across clusters per gene
        expr_z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0)

        if plot and not expr_z.empty:
            plt.figure(figsize=(9, 5))
            sns.heatmap(
                expr_z,
                cmap="vlag",
                center=0,
                cbar_kws={"label": "Z-score (relative expression)"},
            )
            plt.title("Marker gene enrichment (z-scored across clusters)")
            plt.xlabel("Gene")
            plt.ylabel("Cluster ID")
            plt.tight_layout()
            plt.show()

        top_z = expr_z.idxmax(axis=0).rename("max_in_cluster")

    else:
        log.warning("No marker genes found in dataframe columns; skipping enrichment heatmap.")

    return LeidenAnalysisResult(
        df=df,
        adata=adata,
        gene_cols=gene_cols,
        freq_df=freq_df,
        logit_df=logit_df,
        expr_z=expr_z,
        top_z=top_z,
    )
