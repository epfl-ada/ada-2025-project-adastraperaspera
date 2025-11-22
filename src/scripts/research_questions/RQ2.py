from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import logging
from typing import Any

from anndata import AnnData
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import statsmodels.formula.api as smf


@dataclass
class LeidenClusteringResult:
    """
    Outputs for clustering-only pipeline.
    """

    df: pd.DataFrame  # Copy of the input with added 'cluster_leiden'
    adata: AnnData  # AnnData object used by Scanpy (with UMAP computed)
    gene_cols: list[str]  # Detected numeric gene-expression columns


def cluster_cells_leiden(
    combined_df: pd.DataFrame,
    *,
    meta_columns: Iterable[str] | None = None,
    leiden_resolution: float = 0.6,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    logger: logging.Logger | None = None,
) -> LeidenClusteringResult:
    """
    Run Scanpy clustering (scale → PCA → neighbors → Leiden → UMAP) on a single-cell dataframe.
    Only clustering is performed; no plaque analysis or plotting.

    Parameters
    ----------
    combined_df : pd.DataFrame
        Input dataframe with at least a 'cell_id' column and numeric gene-expression columns.
    meta_columns : Iterable[str], optional
        Column names to exclude from gene detection. If None, a sensible default set is used.
    leiden_resolution : float, default=0.6
        Resolution parameter for `sc.tl.leiden`.
    n_neighbors : int, default=15
        Number of neighbors for `sc.pp.neighbors`.
    n_pcs : int, default=30
        Number of principal components to use for neighbors/UMAP.
    logger : logging.Logger, optional
        Logger for progress messages. If None, a module logger is used.

    Returns
    -------
    LeidenClusteringResult
        Dataclass holding the augmented dataframe, AnnData, and detected gene columns.
    """
    log = logger or logging.getLogger(__name__)
    df = combined_df.copy()

    # Defaults for metadata (exclude from gene detection)
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

    # 1) Detect numeric gene columns
    gene_cols: list[str] = [
        c for c in df.columns if c not in meta_cols and pd.api.types.is_numeric_dtype(df[c])
    ]
    if len(gene_cols) == 0:
        raise ValueError("No numeric gene columns detected. Check `meta_columns` or input dtypes.")

    log.info("Detected %d gene columns", len(gene_cols))

    # 2) Build AnnData for Scanpy
    adata = sc.AnnData(df[gene_cols].to_numpy())
    # Minimal obs
    adata.obs = pd.DataFrame(
        {
            "cell_id": df["cell_id"].astype(str).values,
        },
        index=df["cell_id"].astype(str).values,
    )
    adata.obs_names = df["cell_id"].astype(str).values
    adata.var_names = pd.Index(gene_cols, name="genes")

    log.info("AnnData created with shape %s", adata.shape)
    log.info("=== PCA → Neighbors → Leiden Clustering ===")

    # 3) Scanpy pipeline
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden")
    sc.tl.umap(adata)

    log.info("Leiden cluster counts:\n%s", adata.obs["leiden"].value_counts())

    # 4) Add cluster labels back to df
    df["cluster_leiden"] = adata.obs["leiden"].reindex(df["cell_id"].astype(str)).values
    log.info("Added 'cluster_leiden' column to dataframe.")

    return LeidenClusteringResult(df=df, adata=adata, gene_cols=gene_cols)


def plot_leiden_umap_grid(
    adata_by_mouse: Mapping[str, AnnData],
    order: Sequence[str],
    *,
    n_cols: int = 3,
    figsize: tuple[float, float] = (12.0, 8.0),
    legend_loc: str | None = "on data",
    suptitle: str | None = "Leiden clusters across mice",
) -> None:
    """
    Plot a grid of UMAPs colored by Leiden clusters for a given mouse order.
    """
    n = len(order)
    n_cols = max(1, n_cols)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    for i, mouse in enumerate(order):
        r, c = divmod(i, n_cols)
        ax = axes[r][c]
        if mouse not in adata_by_mouse:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{mouse}\n(no data)", ha="center", va="center", fontsize=10)
            continue
        adata = adata_by_mouse[mouse]
        sc.pl.umap(
            adata,
            color=["leiden"],
            legend_loc=legend_loc,
            frameon=False,
            title=str(mouse),
            show=False,
            ax=ax,
        )
    # Hide any unused axes
    for j in range(n, n_rows * n_cols):
        r, c = divmod(j, n_cols)
        axes[r][c].axis("off")
    if suptitle:
        fig.suptitle(suptitle)
    fig.tight_layout()
    plt.show()


def plot_leiden_spatial_grid(
    df_by_mouse: Mapping[str, pd.DataFrame],
    order: Sequence[str],
    *,
    n_cols: int = 3,
    figsize: tuple[float, float] = (12.0, 8.0),
    sample_for_scatter: int | None = 20_000,
    random_state: int = 0,
    suptitle: str | None = "Spatial map of Leiden clusters across mice",
) -> None:
    """
    Plot a grid of spatial scatterplots colored by Leiden clusters for a given mouse order.
    """
    n = len(order)
    n_cols = max(1, n_cols)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    for i, mouse in enumerate(order):
        r, c = divmod(i, n_cols)
        ax = axes[r][c]
        if mouse not in df_by_mouse or df_by_mouse[mouse] is None:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{mouse}\n(no data)", ha="center", va="center", fontsize=10)
            continue
        df = df_by_mouse[mouse]
        required = {"x_centroid", "y_centroid", "cluster_leiden"}
        if not required.issubset(df.columns):
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"{mouse}\n(missing {sorted(required - set(df.columns))})",
                ha="center",
                va="center",
                fontsize=10,
            )
            continue
        plot_df = df
        if sample_for_scatter is not None and sample_for_scatter < len(df):
            plot_df = df.sample(sample_for_scatter, random_state=random_state)
        sns.scatterplot(
            data=plot_df,
            x="x_centroid",
            y="y_centroid",
            hue="cluster_leiden",
            palette="tab20",
            s=6,
            linewidth=0,
            alpha=0.7,
            ax=ax,
            legend=False,
        )
        ax.invert_yaxis()
        ax.set_title(str(mouse))
        ax.set_xlabel("X coordinate (µm)")
        ax.set_ylabel("Y coordinate (µm)")
    # Hide any unused axes
    for j in range(n, n_rows * n_cols):
        r, c = divmod(j, n_cols)
        axes[r][c].axis("off")
    if suptitle:
        fig.suptitle(suptitle)
    fig.tight_layout()
    plt.show()


@dataclass
class JointLeidenClusteringResult:
    """
    Outputs for joint clustering across multiple mice.
    """

    df_by_mouse: dict[str, pd.DataFrame]
    adata: AnnData
    adata_by_mouse: dict[str, AnnData]
    gene_cols: list[str]


def cluster_cells_leiden_joint(
    df_by_mouse: Mapping[str, pd.DataFrame],
    *,
    meta_columns: Iterable[str] | None = None,
    leiden_resolution: float = 0.6,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    logger: logging.Logger | None = None,
) -> JointLeidenClusteringResult:
    """
    Joint Scanpy Leiden clustering across multiple mice with optional BBKNN integration.
    Cluster labels are global and comparable across mice.
    """
    log = logger or logging.getLogger(__name__)
    if not df_by_mouse:
        raise ValueError("`df_by_mouse` is empty – nothing to cluster.")

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
        "mouse",
        "mouse_id",
    }
    meta_cols = set(meta_columns) if meta_columns is not None else default_meta

    per_mouse_gene_cols: dict[str, list[str]] = {}
    per_mouse_dfs: list[pd.DataFrame] = []
    for mouse, df in df_by_mouse.items():
        if df is None or df.empty:
            log.warning("Mouse %s has empty DataFrame; skipping.", mouse)
            continue
        df_local = df.copy()
        if "cell_id" not in df_local.columns:
            raise KeyError(f"DataFrame for mouse '{mouse}' is missing a 'cell_id' column.")
        df_local["cell_id"] = df_local["cell_id"].astype(str)
        df_local["mouse"] = str(mouse)
        gene_cols_mouse = [
            c for c in df_local.columns if c not in meta_cols and is_numeric_dtype(df_local[c])
        ]
        if not gene_cols_mouse:
            raise ValueError(f"No numeric gene columns detected for mouse '{mouse}'.")
        per_mouse_gene_cols[mouse] = gene_cols_mouse
        per_mouse_dfs.append(df_local)
    if not per_mouse_dfs:
        raise ValueError("All DataFrames were empty – no cells to cluster.")

    gene_sets = [set(cols) for cols in per_mouse_gene_cols.values()]
    gene_cols = sorted(set.intersection(*gene_sets))
    if not gene_cols:
        raise ValueError(
            "No overlapping gene columns across mice. "
            "Check that all DataFrames share the same gene panel and `meta_columns`."
        )
    if any(set(cols) != gene_sets[0] for cols in gene_sets[1:]):
        log.warning(
            "Gene columns differ between mice; using intersection of %d genes.", len(gene_cols)
        )

    log.info(
        "Joint clustering across %d mice using %d gene columns.", len(per_mouse_dfs), len(gene_cols)
    )

    combined_df = pd.concat(per_mouse_dfs, axis=0, ignore_index=False)
    combined_df["cell_uid"] = (
        combined_df["mouse"].astype(str) + "__" + combined_df["cell_id"].astype(str)
    )
    combined_df = combined_df.set_index("cell_uid", drop=False)
    if combined_df.index.has_duplicates:
        raise ValueError(
            "Non-unique 'cell_uid' index after concatenation. "
            "Check that 'cell_id' is unique per mouse."
        )

    X = combined_df[gene_cols].to_numpy(dtype=np.float32)
    adata = sc.AnnData(X)
    adata.obs = combined_df[["cell_id", "mouse"]].copy()
    adata.obs_names = combined_df.index.astype(str)
    adata.var_names = pd.Index(gene_cols, name="genes")

    log.info("Global AnnData created with shape %s", adata.shape)
    log.info("=== Joint PCA → Neighbors (batch-aware) → Leiden → UMAP ===")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sce.pp.bbknn(adata, batch_key="mouse", n_pcs=n_pcs, neighbors_within_batch=n_neighbors)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden")
    sc.tl.umap(adata)

    log.info("Global Leiden cluster counts:\n%s", adata.obs["leiden"].value_counts())

    combined_df["cluster_leiden"] = adata.obs["leiden"].reindex(combined_df.index).values

    df_out: dict[str, pd.DataFrame] = {}
    adata_out: dict[str, AnnData] = {}
    for mouse in df_by_mouse.keys():
        mask = combined_df["mouse"] == str(mouse)
        if not mask.any():
            log.warning("No cells for mouse '%s' in combined data.", mouse)
            continue
        df_mouse = combined_df.loc[mask].copy()
        df_out[mouse] = df_mouse
        adata_mouse = adata[adata.obs["mouse"] == str(mouse)].copy()
        adata_out[mouse] = adata_mouse
        log.info(
            "Mouse %s: %d cells, Leiden clusters: %s",
            mouse,
            adata_mouse.n_obs,
            adata_mouse.obs["leiden"].value_counts().to_dict(),
        )

    return JointLeidenClusteringResult(
        df_by_mouse=df_out,
        adata=adata,
        adata_by_mouse=adata_out,
        gene_cols=gene_cols,
    )


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
    freq_trend_df: pd.DataFrame  # OLS coefficients for frequency trends per cluster
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

    # -----------------------------
    # 1–3) Clustering
    # -----------------------------
    clustering = cluster_cells_leiden(
        df,
        meta_columns=meta_columns,
        leiden_resolution=leiden_resolution,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        logger=log,
    )
    df = clustering.df
    adata = clustering.adata
    gene_cols = clustering.gene_cols

    # Align by cell_id to be robust to any ordering differences
    if "cell_id" in df.columns:
        df_idx = df["cell_id"].astype(str)
        if "distance_to_plaque" in df.columns:
            dist_series = pd.to_numeric(df["distance_to_plaque"], errors="coerce")
            dist_aligned = pd.Series(dist_series.values, index=df_idx).reindex(adata.obs_names)
            adata.obs["distance_to_plaque"] = dist_aligned.values
        if "dist_bin_simple" in df.columns:
            bin_series = df["dist_bin_simple"].astype(str)
            bin_aligned = pd.Series(bin_series.values, index=df_idx).reindex(adata.obs_names)
            adata.obs["dist_bin_simple"] = bin_aligned.values

    # UMAP plot of clusters
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
    # Proportion in [0,1] for OLS trends
    freq_df["prop"] = freq_df["pct"].astype(float) / 100.0

    # -----------------------------------------------
    # 6) OLS trends per cluster: prop ~ bin_mid
    # -----------------------------------------------
    ols_rows: list[dict] = []
    for cl in sorted(freq_df["cluster_leiden"].dropna().unique(), key=str):
        sub = freq_df.loc[
            freq_df["cluster_leiden"] == cl, ["prop", "bin_mid", "total_cells"]
        ].dropna()
        if len(sub) < 2:
            ols_rows.append(
                {
                    "cluster": str(cl),
                    "intercept": np.nan,
                    "slope": np.nan,
                    "stderr_slope": np.nan,
                    "pval_slope": np.nan,
                    "r2": np.nan,
                    "n_bins": int(len(sub)),
                }
            )
            continue
        try:
            model = smf.ols("prop ~ bin_mid", data=sub).fit()
            intercept = float(model.params.get("Intercept", np.nan))
            slope = float(model.params.get("bin_mid", np.nan))
            stderr_slope = float(model.bse.get("bin_mid", np.nan))
            pval_slope = float(model.pvalues.get("bin_mid", np.nan))
            r2 = float(model.rsquared)
            ols_rows.append(
                {
                    "cluster": str(cl),
                    "intercept": intercept,
                    "slope": slope,
                    "stderr_slope": stderr_slope,
                    "pval_slope": pval_slope,
                    "r2": r2,
                    "n_bins": int(model.nobs),
                }
            )
        except Exception as e:
            ols_rows.append(
                {
                    "cluster": str(cl),
                    "intercept": np.nan,
                    "slope": np.nan,
                    "stderr_slope": np.nan,
                    "pval_slope": np.nan,
                    "r2": np.nan,
                    "n_bins": int(len(sub)),
                    "error": str(e),
                }
            )
    freq_trend_df = pd.DataFrame(ols_rows)

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
            intercept = float(model.params.get("Intercept", np.nan))
            slope = float(model.params.get("distance_to_plaque", np.nan))
            se_intercept = float(model.bse.get("Intercept", np.nan))
            se_slope = float(model.bse.get("distance_to_plaque", np.nan))
            z_slope = float(model.tvalues.get("distance_to_plaque", np.nan))
            pval = float(model.pvalues.get("distance_to_plaque", np.nan))
            ci = model.conf_int()
            ci_lo = (
                float(ci.loc["distance_to_plaque", 0])
                if "distance_to_plaque" in ci.index
                else np.nan
            )
            ci_hi = (
                float(ci.loc["distance_to_plaque", 1])
                if "distance_to_plaque" in ci.index
                else np.nan
            )
        except Exception as e:
            log.warning("Logit failed for cluster %s (%s). Setting slope/pval = NaN.", cl, e)
            intercept, slope, se_intercept, se_slope, z_slope, pval, ci_lo, ci_hi = (
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            )

        results.append(
            {
                "cluster": str(cl),
                "intercept": intercept,
                "slope": slope,
                "se_intercept": se_intercept,
                "se_slope": se_slope,
                "z_slope": z_slope,
                "pval": pval,
                "ci95_lo_slope": ci_lo,
                "ci95_hi_slope": ci_hi,
            }
        )

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
        freq_trend_df=freq_trend_df,
        logit_df=logit_df,
        expr_z=expr_z,
        top_z=top_z,
    )
