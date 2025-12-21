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
from scipy.optimize import linear_sum_assignment
import seaborn as sns
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
import statsmodels.formula.api as smf
import matplotlib.patches as mpatches


def evaluate_cluster_alignment(
    clustered_df_by_mouse_joint: dict,
    available_mouse_order: list,
    instructor_csv_path: str,
    tissue_map: dict,
    type_to_expanded_type: dict | None = None,
):
    """
    Evaluate clustering results against instructor labels using overlap, Jaccard/F1 similarity,
    Hungarian assignment, and adjusted index metrics.

    Parameters
    ----------
    clustered_df_by_mouse_joint : dict
        Dictionary mapping mouse-age identifiers to DataFrames containing "cell_id" and "cluster_leiden".
    available_mouse_order : list
        List of keys indicating the order in which to process the mice.
    instructor_csv_path : str
        Path to the instructor CSV containing columns ['tissue', 'cell_id', 'predicted_label', 'cell_type'].
    tissue_map : dict
        Mapping from mouse_age to tissue name.
    type_to_expanded_type : dict, optional
        Mapping of compact cell type names to expanded descriptive names.

    Returns
    -------
    dict
        Contains merged tables, contingency matrices, mappings, metrics, and optionally expanded label types.
    """

    my_parts = []
    for mouse_age in available_mouse_order:
        df = clustered_df_by_mouse_joint[mouse_age][
            ["cell_id", "cluster_leiden"]
        ].copy()
        df["mouse_age"] = mouse_age
        df["tissue"] = tissue_map[mouse_age]
        my_parts.append(df)

    my_all = pd.concat(my_parts, ignore_index=True)
    dup_key = my_all.duplicated(["tissue", "cell_id"]).mean()
    if dup_key > 0:
        print(
            f"WARNING: {dup_key:.2%} duplicate (tissue, cell_id) in the input table. Investigate upstream."
        )

    my_all["cluster_leiden"] = pd.to_numeric(
        my_all["cluster_leiden"], errors="raise"
    ).astype(int)

    instr = pd.read_csv(instructor_csv_path)
    instr = instr[instr["tissue"].isin(set(tissue_map.values()))].copy()
    instr["predicted_label"] = pd.to_numeric(
        instr["predicted_label"], errors="raise"
    ).astype(int)

    label_to_celltype = (
        instr.groupby("predicted_label")["cell_type"]
        .agg(lambda x: x.mode().iat[0] if not x.mode().empty else x.iloc[0])
        .to_dict()
    )

    merged = my_all.merge(
        instr[["tissue", "cell_id", "predicted_label", "cell_type"]],
        on=["tissue", "cell_id"],
        how="inner",
    )

    cont = (
        pd.crosstab(
            merged["cluster_leiden"].astype(int),
            merged["predicted_label"].astype(int),
        )
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    row_sum = cont.sum(axis=1).astype(float)
    col_sum = cont.sum(axis=0).astype(float)

    pairs = cont.stack().rename("overlap").reset_index()
    pairs.columns = ["my_cluster", "instr_cluster", "overlap"]

    pairs = pairs.merge(
        row_sum.rename("my_size"), left_on="my_cluster", right_index=True
    )
    pairs = pairs.merge(
        col_sum.rename("instr_size"), left_on="instr_cluster", right_index=True
    )

    den_j = pairs["my_size"] + pairs["instr_size"] - pairs["overlap"]
    pairs["jaccard"] = pairs["overlap"] / np.where(den_j == 0, np.nan, den_j)
    den_f = pairs["my_size"] + pairs["instr_size"]
    pairs["f1"] = 2.0 * pairs["overlap"] / np.where(den_f == 0, np.nan, den_f)

    best_instr_to_my = (
        pairs.sort_values(
            ["instr_cluster", "jaccard", "overlap"], ascending=[True, False, False]
        )
        .groupby("instr_cluster", as_index=False)
        .head(1)
        .copy()
    )
    best_instr_to_my["instr_cell_type"] = best_instr_to_my["instr_cluster"].map(
        label_to_celltype
    )

    cost = -cont.to_numpy()
    r_ind, c_ind = linear_sum_assignment(cost)
    one_to_one = pd.DataFrame(
        {
            "my_cluster": cont.index.to_numpy()[r_ind],
            "instr_cluster": cont.columns.to_numpy()[c_ind],
        }
    )
    one_to_one["overlap"] = [
        cont.loc[i, j]
        for i, j in zip(
            one_to_one["my_cluster"], one_to_one["instr_cluster"], strict=False
        )
    ]
    one_to_one["my_size"] = one_to_one["my_cluster"].map(row_sum)
    one_to_one["instr_size"] = one_to_one["instr_cluster"].map(col_sum)
    one_to_one["jaccard"] = one_to_one["overlap"] / (
        one_to_one["my_size"] + one_to_one["instr_size"] - one_to_one["overlap"]
    ).replace(0, np.nan)
    one_to_one["f1"] = (
        2
        * one_to_one["overlap"]
        / (one_to_one["my_size"] + one_to_one["instr_size"]).replace(0, np.nan)
    )
    one_to_one["instr_cell_type"] = one_to_one["instr_cluster"].map(label_to_celltype)

    ari = adjusted_rand_score(merged["predicted_label"], merged["cluster_leiden"])
    ami = adjusted_mutual_info_score(
        merged["predicted_label"], merged["cluster_leiden"]
    )

    instr_to_my_map = dict(
        zip(
            best_instr_to_my["instr_cluster"],
            best_instr_to_my["my_cluster"],
            strict=False,
        )
    )
    my_to_instr_1to1_map = dict(
        zip(one_to_one["my_cluster"], one_to_one["instr_cluster"], strict=False)
    )

    my_to_instr_argmax = cont.idxmax(axis=1).to_dict()
    my_all["instructor_label_argmax_within_my_cluster"] = my_all["cluster_leiden"].map(
        my_to_instr_argmax
    )

    label_to_type = (
        instr[["predicted_label", "cell_type"]]
        .drop_duplicates(subset=["predicted_label"])
        .assign(
            _rest=lambda d: d["cell_type"].astype(str).str.split(n=1).str[1].fillna("-")
        )
        .set_index("predicted_label")["_rest"]
        .to_dict()
    )

    if type_to_expanded_type:
        my_label_to_type = {
            k: type_to_expanded_type.get(label_to_type.get(v, "-"), "-")
            for k, v in my_to_instr_1to1_map.items()
        }
    else:
        my_label_to_type = {}

    return {
        "merged": merged,
        "contingency": cont,
        "best_instr_to_my": best_instr_to_my,
        "one_to_one": one_to_one,
        "ari": ari,
        "ami": ami,
        "instr_to_my_map": instr_to_my_map,
        "my_to_instr_1to1_map": my_to_instr_1to1_map,
        "my_label_to_type": my_label_to_type,
        "my_all": my_all,
        "label_to_type": label_to_type,
        "pairs": pairs,
    }


@dataclass
class LeidenClusteringResult:
    """
    Outputs for clustering-only pipeline.
    """

    df: pd.DataFrame
    adata: AnnData
    gene_cols: list[str]


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

    gene_cols: list[str] = [
        c
        for c in df.columns
        if c not in meta_cols and pd.api.types.is_numeric_dtype(df[c])
    ]
    if len(gene_cols) == 0:
        raise ValueError(
            "No numeric gene columns detected. Check `meta_columns` or input dtypes."
        )

    log.info("Detected %d gene columns", len(gene_cols))

    adata = sc.AnnData(df[gene_cols].to_numpy())

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

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden")
    sc.tl.umap(adata)

    log.info("Leiden cluster counts:\n%s", adata.obs["leiden"].value_counts())

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
            ax.text(
                0.5, 0.5, f"{mouse}\n(no data)", ha="center", va="center", fontsize=10
            )
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
    add_legend_column: bool = True,
    legend_title: str = "Leiden cluster",
    palette_name: str = "tab20",
) -> None:
    """
    Plot a grid of spatial scatterplots colored by Leiden clusters for a given mouse order.
    If add_legend_column=True, allocate an extra column on the right with a single legend
    describing the cluster_id → color mapping (shared across all panels).
    """

    required = {"x_centroid", "y_centroid", "cluster_leiden"}
    clusters = []
    for m in order:
        df = df_by_mouse.get(m)
        if df is None or not required.issubset(df.columns):
            continue
        clusters.extend(pd.unique(df["cluster_leiden"]))
    unique_clusters = pd.unique(pd.Series(clusters))

    def _safe_sort_key(x):
        try:
            return (0, float(x))
        except Exception:
            return (1, str(x))

    hue_order = sorted(unique_clusters, key=_safe_sort_key)

    colors = sns.color_palette(palette_name, n_colors=len(hue_order))
    cluster_to_color = {cl: col for cl, col in zip(hue_order, colors)}

    plot_cols = max(1, n_cols)
    total_cols = plot_cols + (1 if add_legend_column else 0)
    n = len(order)
    n_rows = (n + plot_cols - 1) // plot_cols

    if add_legend_column:
        fig_w = figsize[0] * (total_cols / plot_cols)
        fig_h = figsize[1]
        use_figsize = (fig_w, fig_h)
    else:
        use_figsize = figsize

    fig, axes = plt.subplots(n_rows, total_cols, figsize=use_figsize, squeeze=False)

    for i, mouse in enumerate(order):
        r, c = divmod(i, plot_cols)
        ax = axes[r][c]
        df = df_by_mouse.get(mouse)
        if df is None:
            ax.axis("off")
            ax.text(
                0.5, 0.5, f"{mouse}\n(no data)", ha="center", va="center", fontsize=10
            )
            continue

        if not required.issubset(df.columns):
            missing = sorted(required - set(df.columns))
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"{mouse}\n(missing {missing})",
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
            hue_order=hue_order,
            palette=cluster_to_color,
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

    total_plot_slots = n_rows * plot_cols
    for j in range(n, total_plot_slots):
        r, c = divmod(j, plot_cols)
        axes[r][c].axis("off")

    if add_legend_column:

        legend_ax = axes[0][total_cols - 1]
        for row in range(1, n_rows):
            axes[row][total_cols - 1].axis("off")

        handles = [
            mpatches.Patch(color=cluster_to_color[cl], label=str(cl))
            for cl in hue_order
        ]
        legend_ax.legend(
            handles=handles,
            title=legend_title,
            loc="center left",
            frameon=False,
        )
        legend_ax.axis("off")

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
            raise KeyError(
                f"DataFrame for mouse '{mouse}' is missing a 'cell_id' column."
            )
        df_local["cell_id"] = df_local["cell_id"].astype(str)
        df_local["mouse"] = str(mouse)
        gene_cols_mouse = [
            c
            for c in df_local.columns
            if c not in meta_cols and is_numeric_dtype(df_local[c])
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
            "Gene columns differ between mice; using intersection of %d genes.",
            len(gene_cols),
        )

    log.info(
        "Joint clustering across %d mice using %d gene columns.",
        len(per_mouse_dfs),
        len(gene_cols),
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
    sce.pp.bbknn(
        adata, batch_key="mouse", n_pcs=n_pcs, neighbors_within_batch=n_neighbors
    )
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden")
    sc.tl.umap(adata)

    log.info("Global Leiden cluster counts:\n%s", adata.obs["leiden"].value_counts())

    combined_df["cluster_leiden"] = (
        adata.obs["leiden"].reindex(combined_df.index).values
    )

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

    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")
    if prox_thresh <= 0 or distal_thresh <= prox_thresh:
        raise ValueError("Thresholds must satisfy: 0 < prox_thresh < distal_thresh.")

    df_out = df if in_place else df.copy()

    x = df_out[column].astype(float)
    x_clean = x.dropna()

    _logger.info("=== Distance to nearest plaque (µm) summary ===")
    desc = x.describe()
    _logger.info(desc)

    if show_plots:
        plt.figure(figsize=figsize)
        sns.histplot(x_clean, bins=n_bins, kde=kde)
        plt.xlabel("Distance to nearest plaque (µm)")
        plt.ylabel("Cell count")
        plt.title("Distribution of cell distances to nearest plaque (linear scale)")
        plt.tight_layout()
        plt.show()

    q5, q25, q50, q75, q95 = np.percentile(x_clean, [5, 25, 50, 75, 95])
    quantiles: dict[str, float] = {
        "5%": float(q5),
        "25%": float(q25),
        "50% (median)": float(q50),
        "75%": float(q75),
        "95%": float(q95),
    }
    _logger.info(quantiles)

    bins = [0.0, float(prox_thresh), float(distal_thresh), np.inf]
    labels: list[str] = [
        f"proximal (<{prox_thresh:g} µm)",
        f"intermediate ({prox_thresh:g}–{distal_thresh:g} µm)",
        f"distal (>{distal_thresh:g} µm)",
    ]
    df_out[category_col] = pd.cut(
        x, bins=bins, labels=labels, include_lowest=True, right=False
    )

    _logger.info("\n=== Cell counts by proximity class ===")
    counts = df_out[category_col].value_counts().reindex(labels, fill_value=0)
    _logger.info(counts)

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

    df: pd.DataFrame
    adata: AnnData
    gene_cols: list[str]
    freq_df: pd.DataFrame
    freq_trend_df: pd.DataFrame
    logit_df: pd.DataFrame
    expr_z: pd.DataFrame
    top_z: pd.Series


def analyze_leiden_spatial(
    combined_df: pd.DataFrame,
    *,
    meta_columns: Iterable[str] | None = None,
    cluster_col: str = "cluster_leiden",
    bin_width_um: float = 25.0,
    proximal_threshold_um: float = 30.0,
    marker_genes: Mapping[str, Sequence[str]] | None = None,
    plot: bool = True,
    logger: logging.Logger | None = None,
    expanded_types: dict | None = None,
) -> LeidenAnalysisResult:
    """
    Spatial analysis using pre-existing Leiden cluster labels in `combined_df`.

    This function:
      1) Validates required columns (cell_id, distance_to_plaque, and cluster labels).
      2) Detects numeric gene-expression columns (excluding metadata).
      3) Builds an AnnData object and (optionally) computes UMAP for visualization
         via scale → PCA → neighbors → UMAP. No clustering is performed.
      4) Keeps the rest of the original downstream analysis:
         - spatial scatter (optional)
         - cluster frequency vs. distance (binned) + OLS trends
         - per-cluster logistic regression of membership vs. distance
         - z-scored marker enrichment per cluster
    """
    log = logger or logging.getLogger(__name__)
    df = combined_df.copy()

    required_cols = {"cell_id", "distance_to_plaque"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input dataframe is missing required columns: {missing}")

    if cluster_col not in df.columns:
        if "cluster_leiden" in df.columns:
            cluster_col = "cluster_leiden"
        elif "clustering_leiden" in df.columns:
            cluster_col = "clustering_leiden"
        else:
            raise ValueError(
                f"Cluster label column '{cluster_col}' not found, and neither "
                f"'cluster_leiden' nor 'clustering_leiden' is present."
            )

    if "cluster_leiden" not in df.columns or cluster_col != "cluster_leiden":
        df["cluster_leiden"] = df[cluster_col]

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
        "cluster_leiden",
        "clustering_leiden",
    }
    meta_cols = set(meta_columns) if meta_columns is not None else default_meta

    gene_cols: list[str] = [
        c
        for c in df.columns
        if c not in meta_cols and pd.api.types.is_numeric_dtype(df[c])
    ]
    if len(gene_cols) == 0:
        raise ValueError(
            "No numeric gene columns detected. Check `meta_columns` or input dtypes."
        )

    log.info("Detected %d gene columns", len(gene_cols))
    log.info("Using pre-existing cluster labels from column '%s'.", cluster_col)
    log.info("Cluster counts:\n%s", df["cluster_leiden"].value_counts(dropna=False))

    df["cell_id"] = df["cell_id"].astype(str)
    X = df[gene_cols].to_numpy(dtype=np.float32, copy=False)

    adata = sc.AnnData(X)
    adata.obs = pd.DataFrame(index=df["cell_id"].values)
    adata.obs["cell_id"] = df["cell_id"].values
    adata.obs["cluster_leiden"] = df["cluster_leiden"].values
    adata.obs["distance_to_plaque"] = pd.to_numeric(
        df["distance_to_plaque"], errors="coerce"
    ).values
    if "dist_bin_simple" in df.columns:
        adata.obs["dist_bin_simple"] = df["dist_bin_simple"].astype(str).values

    adata.obs_names = df["cell_id"].values

    adata.obs_names_make_unique()
    adata.var_names = pd.Index(gene_cols, name="genes")

    dist_num = pd.to_numeric(df["distance_to_plaque"], errors="coerce")
    max_dist = np.nanpercentile(dist_num, 99)
    bins = np.arange(0, max_dist + bin_width_um, bin_width_um)
    df["distance_bin_cont"] = pd.cut(dist_num, bins=bins)

    freq_df = (
        df.groupby(["distance_bin_cont", "cluster_leiden"])
        .size()
        .reset_index(name="n_cells")
    )
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
    freq_df["prop"] = freq_df["pct"].astype(float) / 100.0

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
            ols_rows.append(
                {
                    "cluster": str(cl),
                    "intercept": float(model.params.get("Intercept", np.nan)),
                    "slope": float(model.params.get("bin_mid", np.nan)),
                    "stderr_slope": float(model.bse.get("bin_mid", np.nan)),
                    "pval_slope": float(model.pvalues.get("bin_mid", np.nan)),
                    "r2": float(model.rsquared),
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

    df["is_proximal"] = (dist_num <= proximal_threshold_um).astype(int)

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
            log.warning(
                "Logit failed for cluster %s (%s). Setting slope/pval = NaN.", cl, e
            )
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
        m = logit_df["pval"].notna().sum()
        logit_df["adj_pval"] = np.minimum(logit_df["pval"] * max(m, 1), 1.0)
    else:
        logit_df["adj_pval"] = np.nan

    if plot and not freq_df.empty:
        freq_df["cluster_leiden_expanded"] = freq_df["cluster_leiden"].map(
            expanded_types
        )
        freq_df["cluster_leiden_int"] = freq_df["cluster_leiden"].astype(int)
        logit_df_sign = logit_df.loc[logit_df["adj_pval"] < 0.01].copy()
        freq_df["cluster_leiden_int"] = pd.to_numeric(
            freq_df["cluster_leiden_int"], errors="raise"
        ).astype("int64")
        logit_df_sign["cluster"] = pd.to_numeric(
            logit_df_sign["cluster"], errors="raise"
        ).astype("int64")

        freq_df_sign = freq_df.merge(
            logit_df_sign,
            how="inner",
            left_on="cluster_leiden_int",
            right_on="cluster",
        )
        freq_df_sign = freq_df.merge(
            logit_df_sign,
            how="inner",
            left_on="cluster_leiden_int",
            right_on="cluster",
        )
        plt.figure(figsize=(10, 6))
        sns.lineplot(
            data=freq_df_sign,
            x="bin_mid",
            y="pct",
            hue="cluster_leiden_expanded",
            marker="o",
            alpha=0.7,
        )
        plt.xlabel("Distance to plaque")
        plt.ylabel("% of cells in each bin")
        plt.legend(bbox_to_anchor=(1.05, 1), title="Cluster", ncol=1)
        plt.tight_layout()
        plt.show()

    if marker_genes is None:
        marker_genes = {
            "microglia": ["Apoe", "Hexb", "Cst3", "C4b"],
            "astrocyte": ["Gfap", "Serpina3n", "Vim", "S100a6"],
            "neuron": ["Nrep", "Snap25", "Rbfox3"],
            "oligodendrocyte": ["Mbp", "Plp1", "Cldn11"],
        }

    genes_to_check = [
        g for genes in marker_genes.values() for g in genes if g in df.columns
    ]
    expr_z = pd.DataFrame()
    top_z = pd.Series(dtype=object)

    if genes_to_check:
        df["cluster_leiden"] = df["cluster_leiden"].map(expanded_types)
        expr = df.groupby("cluster_leiden")[genes_to_check].mean(numeric_only=True)

        expr_z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0)

        if plot and not expr_z.empty:
            plt.figure(figsize=(9, 5))
            sns.heatmap(
                expr_z,
                cmap="vlag",
                center=0,
                cbar_kws={"label": "Z-score"},
            )
            plt.xlabel("")
            plt.ylabel("")
            plt.tight_layout()
            plt.show()

        top_z = expr_z.idxmax(axis=0).rename("max_in_cluster")

    else:
        log.warning(
            "No marker genes found in dataframe columns; skipping enrichment heatmap."
        )

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
