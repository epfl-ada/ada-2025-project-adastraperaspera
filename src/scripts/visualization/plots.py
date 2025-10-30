from __future__ import annotations

from collections.abc import Sequence
from typing import Union

import geopandas as gpd
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPoly
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import spearmanr
import seaborn as sns
from shapely.geometry import MultiPolygon, Polygon
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import multipletests

Number = Union[int, float, np.number]


def plot_model_performance(results_df):
    """Bar plot comparing R² scores across models."""
    plt.figure(figsize=(6, 4))
    results_melted = results_df.melt(
        id_vars="model",
        value_vars=["train_r2", "test_r2"],
        var_name="Dataset",
        value_name="R²",
    )
    sns.barplot(data=results_melted, x="model", y="R²", hue="Dataset", palette="viridis")
    plt.title("Model Performance Comparison (Train vs Test R²)")
    plt.xlabel("Model")
    plt.ylabel("R²")
    plt.legend(title="")
    plt.tight_layout()
    plt.show()


def plot_top_gene_importances(importance_df, top_n=20):
    """Heatmap of top predictive genes across models."""
    # Normalize importance per model
    normed = importance_df.groupby("model", group_keys=False).apply(
        lambda d: d.assign(norm_importance=d["importance"] / d["importance"].max())
    )

    # Take top_n per model
    top_genes = normed.groupby("model", group_keys=False).apply(
        lambda d: d.nlargest(top_n, "norm_importance")
    )

    # Pivot for heatmap
    pivot = top_genes.pivot_table(
        index="gene", columns="model", values="norm_importance", fill_value=0
    )

    # Order genes by average importance
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

    plt.figure(figsize=(10, max(6, top_n * 0.3)))
    sns.heatmap(pivot, cmap="mako", linewidths=0.5, cbar_kws={"label": "Normalized Importance"})
    plt.title(f"Top {top_n} Predictive Genes Across Models")
    plt.xlabel("Model")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.show()


def plot_spatial_overlay(df, gene, alpha=0.7, sample_size=20000):
    """
    Scatter overlay showing per-cell expression intensity in tissue coordinates.
    Colors high expression regions with a semi-transparent viridis map.
    """
    import matplotlib.pyplot as plt

    data = df.sample(min(sample_size, len(df)), random_state=42)
    plt.figure(figsize=(6, 6))
    sc = plt.scatter(
        data["x_centroid"],
        data["y_centroid"],
        c=data[gene],
        cmap="viridis",
        s=6,
        alpha=alpha,
        linewidth=0,
    )
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"{gene} spatial expression map")
    plt.colorbar(sc, label="Expression (log₁₊)")
    plt.tight_layout()
    plt.show()


def plot_predicted_proximity(df, score_col="plaque_proximity_score", sample_size=20000):
    """
    Shows predicted proximity score across tissue space.
    Darker/brighter regions = cells predicted closer to plaques.
    """
    import matplotlib.pyplot as plt

    data = df.sample(min(sample_size, len(df)), random_state=42)
    plt.figure(figsize=(6, 6))
    sc = plt.scatter(
        data["x_centroid"],
        data["y_centroid"],
        c=data[score_col],
        cmap="magma",
        s=6,
        alpha=0.8,
        linewidth=0,
    )
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title("Predicted plaque proximity score")
    plt.colorbar(sc, label="Model-predicted proximity (relative units)")
    plt.tight_layout()
    plt.show()


def plot_multi_gene_signature(df, genes, distance_col="distance_to_plaque"):
    """
    Plots mean expression of multiple genes along plaque distance as a heatmap.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns

    bins = np.linspace(0, df[distance_col].max(), 40)
    df["distance_bin"] = pd.cut(df[distance_col], bins=bins)

    grouped = df.groupby("distance_bin")[genes].mean().reset_index()
    heatmap_data = grouped.set_index("distance_bin")[genes].T

    plt.figure(figsize=(10, 4))
    sns.heatmap(heatmap_data, cmap="rocket_r", cbar_kws={"label": "Mean expression"})
    plt.xlabel("Distance bin (µm)")
    plt.ylabel("Gene")
    plt.title("Spatial gradient of multi-gene expression around plaques")
    plt.tight_layout()
    plt.show()


def plot_gene_near_plaques(df, gene, dist_thresh=30.0, sample_size=20000):
    """
    Highlights cells near plaques in grey and colors only the most predictive gene expression.
    """

    data = df.sample(min(sample_size, len(df)), random_state=42)

    # Split proximal vs distal
    near = data[data["distance_to_plaque"] <= dist_thresh]
    far = data[data["distance_to_plaque"] > dist_thresh]

    plt.figure(figsize=(6, 6))
    # Plot background (all cells)
    plt.scatter(
        far["x_centroid"],
        far["y_centroid"],
        color="lightgrey",
        s=4,
        alpha=0.3,
        linewidth=0,
        label="Distal cells",
    )
    # Overlay plaque-near colored by gene expression
    sc = plt.scatter(
        near["x_centroid"],
        near["y_centroid"],
        c=near[gene],
        cmap="inferno",
        s=8,
        alpha=0.8,
        linewidth=0,
    )
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"{gene}: Expression near plaques (≤ {dist_thresh} µm)")
    cbar = plt.colorbar(sc, label="Expression (log₁₊)")
    plt.tight_layout()
    plt.show()


def plot_spatial_with_plaques(df, plaques_poly, gene=None, score_col=None, sample_size=20000):
    """
    Visualizes either gene expression or model-predicted proximity,
    with plaque outlines overlaid from the plaque_polygons DataFrame.
    - If gene is given: color by gene expression
    - If score_col is given: color by model score
    """

    data = df.sample(min(sample_size, len(df)), random_state=42)

    # Decide which column to plot
    if gene:
        color_values = data[gene]
        title = f"{gene} expression + plaque outlines"
        cmap = "inferno"
    elif score_col:
        color_values = data[score_col]
        title = f"{score_col} + plaque outlines"
        cmap = "magma"
    else:
        raise ValueError("Provide either `gene` or `score_col`.")

    plt.figure(figsize=(7, 7))
    sc = plt.scatter(
        data["x_centroid"],
        data["y_centroid"],
        c=color_values,
        cmap=cmap,
        s=6,
        alpha=0.8,
        linewidth=0,
    )
    plt.gca().invert_yaxis()
    plt.axis("off")

    # Overlay plaque polygons
    try:
        if isinstance(plaques_poly, pd.DataFrame) and "geometry" in plaques_poly.columns:
            gdf = gpd.GeoDataFrame(plaques_poly, geometry="geometry")
        elif isinstance(plaques_poly, gpd.GeoDataFrame):
            gdf = plaques_poly
        else:
            raise ValueError("plaques_poly must contain a 'geometry' column.")

        gdf.boundary.plot(ax=plt.gca(), color="cyan", linewidth=0.7, alpha=0.8, label="Plaques")
    except Exception as e:
        print(f"Could not overlay plaques: {e}")

    plt.title(title)
    plt.colorbar(sc, label="Intensity")
    plt.tight_layout()
    plt.show()


def plot_pred_vs_true(y_true, y_pred, model_name="Model"):
    """
    Scatter plot showing predicted vs. true spatial distances.
    Includes line of identity and R² annotation.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.metrics import r2_score

    r2 = r2_score(y_true, y_pred)

    plt.figure(figsize=(5, 5))
    sns.scatterplot(x=y_true, y=y_pred, s=10, alpha=0.4, color="teal")
    plt.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], "r--", lw=1)
    plt.xlabel("Ground truth distance (µm)")
    plt.ylabel("Predicted distance (µm)")
    plt.title(f"{model_name}: Predicted vs. Actual\nR² = {r2:.3f}")
    plt.tight_layout()
    plt.show()


def plot_residual_hist(y_true, y_pred, model_name="Model"):
    """
    Histogram + KDE of residuals (true - predicted).
    Reveals bias in over/under-prediction.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    residuals = y_true - y_pred
    plt.figure(figsize=(6, 4))
    sns.histplot(residuals, kde=True, color="darkslateblue", bins=50)
    plt.axvline(0, color="r", linestyle="--", linewidth=1)
    plt.xlabel("Residual (True − Predicted, µm)")
    plt.ylabel("Cell count")
    plt.title(f"{model_name}: Residual Distribution")
    plt.tight_layout()
    plt.show()


def plot_spatial_residual_map(df, y_true, y_pred, sample_size=20000, model_name="Model"):
    """
    Plots per-cell residuals in spatial coordinates.
    Blue = over-predicted (model thinks it's farther),
    Red = under-predicted (model thinks it's closer).
    """

    y_true_vals = np.asarray(y_true).reshape(-1)
    y_pred_vals = np.asarray(y_pred).reshape(-1)
    n = min(sample_size, len(df))

    sample_idx = np.random.choice(len(df), size=n, replace=False)
    data = df.iloc[sample_idx].copy()

    data["residual"] = y_true_vals[sample_idx] - y_pred_vals[sample_idx]

    plt.figure(figsize=(7, 7))
    sc = plt.scatter(
        data["x_centroid"],
        data["y_centroid"],
        c=data["residual"],
        cmap="coolwarm",
        s=6,
        alpha=0.8,
        linewidth=0,
        vmin=-np.percentile(abs(data["residual"]), 99),
        vmax=np.percentile(abs(data["residual"]), 99),
    )
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"{model_name}: Spatial Residual Map")
    plt.colorbar(sc, label="True − Predicted distance (µm)")
    plt.tight_layout()
    plt.show()


def plot_residual_figure(
    df,
    y_true,
    y_pred,
    gene_cols,
    plaques_poly=None,
    oligo_marker="Plp1",
    sample_size=20000,
):
    """
    Composite figure for analyzing residual structure and cell-type bias.
    Panels:
      A: Predicted vs True distance, colored by oligodendrocyte marker
      B: Spatial residual map (with optional plaque outlines)
      C: Mean |residual| by major cell-type markers
      D: correlation heatmap for top residual-associated genes
    """

    y_true_vals = np.asarray(y_true)
    y_pred_vals = np.asarray(y_pred)
    residual = y_true_vals - y_pred_vals
    abs_resid = np.abs(residual)
    r2 = r2_score(y_true_vals, y_pred_vals)
    df = df.copy()
    df["residual"] = residual
    df["abs_residual"] = abs_resid

    n = min(sample_size, len(df))
    sample_idx = np.random.choice(len(df), size=n, replace=False)
    data = df.iloc[sample_idx]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    plt.subplots_adjust(wspace=0.35)
    sns.set_style("white")

    plp_expr = np.log1p(data[oligo_marker])
    sc = axes[0].scatter(
        y_true_vals[sample_idx],
        y_pred_vals[sample_idx],
        c=plp_expr,
        cmap="viridis",
        s=10,
        alpha=0.6,
        linewidth=0,
    )
    lims = [0, max(y_true_vals.max(), y_pred_vals.max())]
    axes[0].plot(lims, lims, "r--", lw=1)
    axes[0].set_xlabel("True distance (µm)")
    axes[0].set_ylabel("Predicted distance (µm)")
    axes[0].set_title(f"(A) Predicted vs True (colored by {oligo_marker})\nR² = {r2:.3f}")
    fig.colorbar(sc, ax=axes[0], label=f"{oligo_marker} expression (log₁₊)")

    res = data["residual"]
    sc2 = axes[1].scatter(
        data["x_centroid"],
        data["y_centroid"],
        c=res,
        cmap="coolwarm",
        s=6,
        alpha=0.8,
        linewidth=0,
        vmin=-np.percentile(abs(res), 99),
        vmax=np.percentile(abs(res), 99),
    )
    axes[1].invert_yaxis()
    axes[1].axis("off")
    axes[1].set_title("(B) Spatial residual map")
    fig.colorbar(sc2, ax=axes[1], label="True − Predicted (µm)")

    if plaques_poly is not None:
        try:
            import geopandas as gpd

            if "geometry" not in plaques_poly.columns:
                raise ValueError("plaques_poly must include a 'geometry' column.")
            gpd.GeoDataFrame(plaques_poly, geometry="geometry").boundary.plot(
                ax=axes[1], color="cyan", linewidth=0.7, alpha=0.7
            )
        except Exception as e:
            print(f"⚠️ Plaque overlay skipped: {e}")

    marker_genes = {
        "Astrocyte (Gfap)": "Gfap",
        "Microglia (C1qa)": "C1qa",
        "Neuron (Slc17a7)": "Slc17a7",
        "Oligodendrocyte (Plp1)": oligo_marker,
    }

    mean_resids = {}
    for label, g in marker_genes.items():
        mask = df[g] > np.percentile(df[g], 75)
        mean_resids[label] = df.loc[mask, "abs_residual"].mean()

    bars = pd.Series(mean_resids).sort_values(ascending=False)
    sns.barplot(x=bars.values, y=bars.index, palette="crest", ax=axes[2])
    axes[2].set_xlabel("Mean |Residual| (µm)")
    axes[2].set_ylabel("")
    axes[2].set_title("(C) Mean residual by cell-type marker")

    plt.suptitle("Residual Structure and Cell-type Bias", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.show()

    return {"r2": r2, "mean_resids": bars}


def plot_plaques(
    df,
    *,
    brain_geom: Polygon | None = None,
    sample_hulls: int = 20,
    seed: int = 42,
    figsize=(8, 8),
    ax: plt.Axes | None = None,
):
    """
    Plot plaque geometries with styling by convexity and sampled convex hulls.

    Parameters
    ----------
    df : pandas.DataFrame
        Must have a 'geometry' column with shapely (Multi)Polygon objects, and an
        'is_convex' boolean column. 'plaque_id' and 'area' will be created if missing.
    brain_geom : shapely Polygon, optional
        If provided and valid, its outline is drawn as the Brain ROI.
    sample_hulls : int
        Number of random plaques to overlay convex hulls for (capped by len(df)).
    seed : int
        Random seed for hull sampling.
    figsize : tuple
        Matplotlib figure size if `ax` is not provided.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on. A new figure/axes is created if None.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if "geometry" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'geometry' column.")

    # Work on a shallow copy to avoid mutating caller's frame
    P = df.copy()

    if "plaque_id" not in P.columns:
        P["plaque_id"] = np.arange(1, len(P) + 1, dtype=int)
    if "area" not in P.columns:
        P["area"] = P["geometry"].map(lambda g: getattr(g, "area", np.nan))

    # Create axes if needed
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True

    # ---- internal helpers -------------------------------------------------

    def _iter_polygons(g):
        """Yield Polygon objects from Polygon or MultiPolygon; ignore others."""
        if isinstance(g, Polygon):
            yield g
        elif isinstance(g, MultiPolygon):
            for sub in g.geoms:
                if isinstance(sub, Polygon):
                    yield sub

    def _add_geom(
        g,
        *,
        facecolor="none",
        edgecolor="black",
        alpha=0.6,
        linestyle="-",
        linewidth=0.8,
    ):
        """Add a (Multi)Polygon geometry to axes."""
        for poly in _iter_polygons(g):
            # Skip degenerate or invalid rings gracefully
            if poly.is_empty or not poly.is_valid or poly.exterior is None:
                continue
            x, y = poly.exterior.xy
            ax.add_patch(
                MplPoly(
                    list(zip(x, y, strict=False)),
                    closed=True,
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                    linestyle=linestyle,
                    alpha=alpha,
                )
            )

    # ---- draw brain ROI (optional) ----------------------------------------
    if brain_geom is not None and getattr(brain_geom, "is_valid", False):
        bx, by = brain_geom.exterior.xy
        ax.plot(bx, by, color="blue", lw=1.0, label="Brain ROI")

    # ---- draw plaques by convexity ----------------------------------------
    # Use itertuples for speed and avoid repeated attribute lookups
    for row in P.itertuples(index=False):
        g = getattr(row, "geometry", None)
        is_convex = getattr(row, "is_convex", None)
        if g is None:
            continue

        if bool(is_convex):
            _add_geom(
                g,
                facecolor="none",
                edgecolor="green",
                linestyle="-",
                linewidth=0.8,
            )
        else:
            _add_geom(
                g,
                facecolor="none",
                edgecolor="red",
                linestyle="--",
                linewidth=0.8,
            )

    # ---- sample convex hull overlays --------------------------------------
    n = min(len(P), int(sample_hulls))
    if n > 0:
        hull_sample = P.sample(n=n, random_state=seed)
        for row in hull_sample.itertuples(index=False):
            g = getattr(row, "geometry", None)
            if g is None:
                continue
            _add_geom(
                g.convex_hull,
                facecolor="none",
                edgecolor="orange",
                linestyle=":",
                linewidth=1.0,
            )

    # ---- cosmetics ---------------------------------------------------------
    ax.set_aspect("equal", "box")
    ax.set_title(
        "Plaque geometries after normalization\n"
        "Green = convex, Red dashed = non-convex, Orange dotted = convex hulls",
        fontsize=11,
    )
    ax.set_xlabel("X coordinate (µm)")
    ax.set_ylabel("Y coordinate (µm)")

    # Build a clean legend with proxy artists (avoids duplicate entries)
    legend_elems = [
        Line2D([0], [0], color="blue", lw=1.0, label="Brain ROI"),
        Line2D([0], [0], color="green", lw=0.8, linestyle="-", label="Convex"),
        Line2D([0], [0], color="red", lw=0.8, linestyle="--", label="Non-convex"),
        Line2D(
            [0],
            [0],
            color="orange",
            lw=1.0,
            linestyle=":",
            label="Convex hull (sample)",
        ),
    ]
    # Only include items that were actually drawn
    handles, labels = [], []
    if brain_geom is not None and getattr(brain_geom, "is_valid", False):
        handles.append(legend_elems[0])
        labels.append(legend_elems[0].get_label())
    handles.extend(legend_elems[1:])
    labels.extend([e.get_label() for e in legend_elems[1:]])
    ax.legend(handles, labels, loc="upper right", frameon=False)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax


def analyze_plaque_distance(
    cells_with_distances: pd.DataFrame,
    *,
    column: str = "nearest_plaque_dist",
    prox_thresh: float = 30.0,
    distal_thresh: float = 100.0,
    out_col: str = "dist_bin_proximal_distal",
    bins: int = 60,
    figsize: tuple[float, float] = (7, 3),
    ax: plt.Axes | None = None,
    logger: object | None = None,
    inplace: bool = True,
):
    """
    Compute summary stats, plot a histogram, and bin distances into
    {'proximal','intermediate','distal'} based on thresholds.

    Parameters
    ----------
    cells_with_distances : pd.DataFrame
        Input dataframe containing the distance column.
    column : str
        Name of the distance column.
    prox_thresh : float
        Distance <= prox_thresh -> 'proximal'.
    distal_thresh : float
        Distance >= distal_thresh -> 'distal'.
    out_col : str
        Name of the output bin column to create.
    bins : int
        Number of histogram bins.
    figsize : (w, h)
        Figure size if `ax` is None.
    ax : matplotlib.axes.Axes or None
        Existing axes to plot on. A new figure/axes is created if None.
    logger : object or None
        Logger with an `.info()` method. If None, prints to stdout.
    inplace : bool
        If True, add `out_col` to `cells_with_distances` in place. Otherwise return a copy.

    Returns
    -------
    df : pd.DataFrame
        The dataframe (original or a copy) with `out_col` added.
    ax : matplotlib.axes.Axes
        The axes containing the histogram.
    """
    if column not in cells_with_distances.columns:
        raise KeyError(f"'{column}' not found in dataframe columns.")

    # Work on either the original df or a copy
    df = cells_with_distances if inplace else cells_with_distances.copy()

    # Drop NaNs for stats/plot only
    d = df[column].dropna()

    # Stats
    stats = {
        "min": float(d.min()) if len(d) else np.nan,
        "median": float(d.median()) if len(d) else np.nan,
        "mean": float(d.mean()) if len(d) else np.nan,
        "max": float(d.max()) if len(d) else np.nan,
        "n": int(len(d)),
        "n_nan": int(df[column].isna().sum()),
    }
    if logger is not None:
        logger.info(stats)
    else:
        print(stats)

    # Plot
    created_fig = False
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
        created_fig = True

    ax.hist(d, bins=bins, alpha=0.8)
    ax.set_xlabel("Nearest plaque boundary distance")
    ax.set_ylabel("Count")
    ax.set_title("Cells: boundary distance to nearest plaque")
    if created_fig:
        plt.tight_layout()
        plt.show()

    # Binning
    def _bin(v: float):
        if pd.isna(v):
            return np.nan
        if v <= prox_thresh:
            return "proximal"
        if v >= distal_thresh:
            return "distal"
        return "intermediate"

    df[out_col] = df[column].map(_bin)

    # Report bin counts
    bin_counts = df[out_col].value_counts(dropna=False)
    if logger is not None:
        logger.info({"bin_counts": bin_counts.to_dict()})
    else:
        print({"bin_counts": bin_counts.to_dict()})

    return df, ax


def _apply_axis_formatting(
    ax: plt.Axes,
    *,
    xlim: tuple[Number, Number] | None = None,
    ylim: tuple[Number, Number] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    xscale: str = "linear",
    yscale: str = "linear",
    grid: bool = True,
    tight: bool = True,
    rotate_xticks: int | None = None,
) -> None:
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if xscale:
        ax.set_xscale(xscale)
    if yscale:
        ax.set_yscale(yscale)
    if grid:
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)
    if rotate_xticks:
        ax.tick_params(axis="x", rotation=rotate_xticks)
    if tight:
        plt.tight_layout()


def hist1d(
    data: Sequence[Number] | np.ndarray,
    *,
    bins: int = 60,
    density: bool = False,
    thresholds: Sequence[tuple[Number, str, str]] | None = None,
    # thresholds: list of tuples (x_value, color, linestyle)
    figsize: tuple[int, int] = (7, 4),
    alpha: float = 0.85,
    xlim: tuple[Number, Number] | None = None,
    ylim: tuple[Number, Number] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    xscale: str = "linear",
    yscale: str = "linear",
) -> plt.Axes:
    data = np.asarray(data)
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(data, bins=bins, alpha=alpha, density=density)
    if thresholds:
        for x, color, ls in thresholds:
            ax.axvline(x, color=color, linestyle=ls, linewidth=1.0)
    _apply_axis_formatting(
        ax,
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        xscale=xscale,
        yscale=yscale,
    )
    return ax


def line_with_ci(
    x: Sequence,
    y: Sequence[Number],
    *,
    yerr: Sequence[Number] | None = None,
    marker: str = "o",
    linewidth: float = 1.0,
    capsize: float = 3.0,
    figsize: tuple[int, int] = (6, 4),
    xlim: tuple[Number, Number] | None = None,
    ylim: tuple[Number, Number] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    rotate_xticks: int | None = 0,
    xscale: str = "linear",
    yscale: str = "linear",
) -> plt.Axes:
    fig, ax = plt.subplots(figsize=figsize)
    if yerr is None:
        ax.plot(x, y, marker=marker, linewidth=linewidth)
    else:
        ax.errorbar(x, y, yerr=yerr, marker=marker, linewidth=linewidth, capsize=capsize)
    _apply_axis_formatting(
        ax,
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        xscale=xscale,
        yscale=yscale,
        rotate_xticks=rotate_xticks,
    )
    return ax


def scatter2d(
    x: Sequence[Number],
    y: Sequence[Number],
    *,
    s: float = 8.0,
    alpha: float = 0.7,
    figsize: tuple[int, int] = (6, 5),
    xlim: tuple[Number, Number] | None = None,
    ylim: tuple[Number, Number] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    xscale: str = "linear",
    yscale: str = "linear",
) -> plt.Axes:
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x, y, s=s, alpha=alpha)
    _apply_axis_formatting(
        ax,
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        xscale=xscale,
        yscale=yscale,
    )
    return ax


def violin_plot(
    df: pd.DataFrame,
    col: str,
    title: str = None,
    xlabel: str = None,
    ylabel: str = None,
):
    """
    Create a violin plot for a specified column in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe containing the data.
    col : str
        The column name in the dataframe to visualize.
    title : str, optional
        The plot title. Defaults to None.
    xlabel : str, optional
        The label for the x-axis. Defaults to the column name if None.
    ylabel : str, optional
        The label for the y-axis. Defaults to "Density" if None.
    """
    # Validate column
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in DataFrame")

    # Plot style
    sns.set(style="whitegrid", palette="pastel")

    # Create the figure
    plt.figure(figsize=(8, 5))
    sns.violinplot(y=df[col], inner="box", cut=0)

    # Titles and labels
    plt.title(title or f"Distribution of {col}", fontsize=14, pad=12)
    plt.xlabel(xlabel or "", fontsize=12)
    plt.ylabel(ylabel or col, fontsize=12)

    # Clean layout
    plt.tight_layout()
    plt.show()


def cell_type_prop_by_dist(props):
    """
    Plot stacked bar chart of cell type proportions by distance bin.
    """
    pivot_props = props.pivot(
        index="distance_bin", columns="cell_type", values="proportion"
    ).fillna(0)
    pivot_props.plot(kind="bar", stacked=True, figsize=(8, 5), colormap="tab20")
    plt.ylabel("Proportion")
    plt.xlabel("Distance bin (µm)")
    plt.title("Stacked cell type proportions by plaque distance")
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.show()


def cell_type_comp_by_dist(props):
    """
    Plot bar chart of cell type composition by distance bin.
    """
    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=props,
        x="distance_bin",
        y="proportion",
        hue="cell_type",
    )
    plt.title("Cell type composition by distance to plaque")
    plt.xlabel("Distance to plaque (µm, binned)")
    plt.ylabel("Proportion of cells")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title="Cell type")
    plt.tight_layout()
    plt.show()


def mean_exp_dist(mean_melt):
    """
    Plot mean expression per distance bin for PIG genes."""
    g = sns.catplot(
        data=mean_melt,
        x="distance_bin",
        y="mean_expr",
        hue="gene",
        kind="bar",
        height=4,
        aspect=1.6,
    )
    g.set_axis_labels("Distance to plaque (µm, binned)", "Mean expression")
    g.fig.suptitle("Mean PIG expression per distance bin")
    plt.tight_layout()
    plt.show()


def pig_expr_by_dist(mean_by_bin, pig_cols):
    plt.figure(figsize=(9, 5))
    x = np.arange(len(mean_by_bin))
    for g in pig_cols:
        plt.plot(x, mean_by_bin[g], marker="o", linewidth=2, label=g)
    plt.xticks(x, mean_by_bin["distance_bin"].astype(str), rotation=30)
    plt.xlabel("Distance (bins)")
    plt.ylabel("Mean expression (log1p)")
    plt.title("PIGs : mean expression by distance bin")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, ncol=1)
    plt.tight_layout()
    plt.show()


def pig_comp_heatmap(pig_mat, prop_mat, pig_cols):
    corrs = pd.DataFrame(index=pig_cols, columns=prop_mat.columns, dtype=float)
    pvals = pd.DataFrame(index=pig_cols, columns=prop_mat.columns, dtype=float)
    for g in pig_cols:
        y = pig_mat[g].to_numpy()
        for ct in prop_mat.columns:
            x = prop_mat[ct].to_numpy()
            if len(y) >= 2:
                r, p = spearmanr(y, x, nan_policy="omit")
            else:
                r, p = (np.nan, np.nan)
            corrs.loc[g, ct] = r
            pvals.loc[g, ct] = p

    # FDR
    mask = np.isfinite(pvals.values)
    flat = pvals.values[mask]
    rej, qvals, *_ = multipletests(flat, method="fdr_bh")
    q = pvals.copy()
    q.values[mask] = qvals

    plt.figure(figsize=(min(14, 6 + 0.25 * len(corrs.columns)), 8))
    sns.heatmap(
        corrs.astype(float),
        cmap="coolwarm",
        center=0,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Spearman ρ"},
    )
    plt.title("PIG correlation ↔ cellular type proportion (per bins distance)")
    plt.xlabel("Cellular type")
    plt.ylabel("PIG Gene")
    plt.tight_layout()
    plt.show()


def add_plaques_to_plotly(fig, plaques_gdf, name="Plaques", line_color="lime", line_width=2):
    def _add_ring(ring, fig):
        xs, ys = ring.xy  # -> array('d', ...)
        fig.add_trace(
            go.Scatter(
                x=list(xs),
                y=list(ys),  # conversion
                mode="lines",
                line=dict(color=line_color, width=line_width),
                name=name,
                hoverinfo="skip",
                showlegend=False,  # avoid legend duplication
            )
        )

    for geom in plaques_gdf.geometry.dropna():
        if geom.is_empty:
            continue
        gtype = geom.geom_type
        if gtype == "Polygon":
            _add_ring(geom.exterior, fig)
            for interior in geom.interiors:  # holes
                _add_ring(interior, fig)
        elif gtype == "MultiPolygon":
            for poly in geom.geoms:
                _add_ring(poly.exterior, fig)
                for interior in poly.interiors:
                    _add_ring(interior, fig)
        elif gtype in ("LineString", "LinearRing"):
            _add_ring(geom, fig)
        else:
            # fallback: try to access exterior if possible
            if hasattr(geom, "exterior") and geom.exterior is not None:
                _add_ring(geom.exterior, fig)

    return fig


def plot_plaques_dist(cells_with_dist, plaques_gdf):
    fig = px.scatter(
        cells_with_dist,
        x="x_centroid",
        y="y_centroid",
        color="distance_to_plaque",
        color_continuous_scale="plasma",
        title="Cell-to-plaque distance map",
        width=900,
        height=800,
        render_mode="webgl",  # faster for large datasets
    )

    fig.update_traces(marker=dict(size=3), selector=dict(mode="markers"))
    fig = add_plaques_to_plotly(fig, plaques_gdf, name="Plaques", line_color="lime", line_width=2)

    fig.update_yaxes(scaleanchor="x", scaleratio=1)

    # export
    fig.write_html("src/data/figures/cell_to_plaque_map_interactive.html")
    # fig.write_image("src/data/figures/cell_to_plaque_map_interactive.png", scale=2)
    fig.show()


def _topn_args(df, metrics, vis_map, n, p_col, fdr_col, title):
    """Build the (args) for the slider step to update data for all metric traces."""
    new_data = []
    for m in metrics:
        if m not in df.columns:
            # placeholder (won't be visible anyway)
            new_data.append({})
            continue
        top = df.nlargest(int(n), m).copy().sort_values(m, ascending=True)
        hover = f"<b>%{{y}}</b><br>{m}: %{{x:.4g}}"
        if p_col:
            hover += "<br>p: %{customdata[0]:.2e}"
        if fdr_col:
            hover += "<br>FDR: %{customdata[1]:.2e}"
        custom = (
            np.stack(
                [
                    top[p_col].to_numpy() if p_col else np.full(len(top), np.nan),
                    top[fdr_col].to_numpy() if fdr_col else np.full(len(top), np.nan),
                ],
                axis=1,
            )
            if (p_col or fdr_col)
            else None
        )

        new_data.append(
            {
                "x": [top[m].to_numpy()],
                "y": [top["gene"].to_numpy()],
                "customdata": [custom] if custom is not None else [None],
                "hovertemplate": [hover],
                "marker": [
                    dict(
                        color=np.where(top[m].to_numpy() >= 0, "rgb(31,120,180)", "rgb(227,26,28)")
                    )
                ],
            }
        )
    # visibility mask stays the same; layout title & xaxis will be set by the dropdown
    return [{"data": new_data}, {}]


def draw_figures(plot_func, img_pth="std.png", *args, **kwargs):
    fig = plot_func(*args, **kwargs)
    fig.write_html(img_pth)
    fig.show()


def cellular_comp_by_dist(pivot_prop):
    ax = pivot_prop.plot(kind="bar", stacked=True, figsize=(8, 5), width=0.85, colormap="tab20")
    ax.set_xlabel("Distance to plaque (µm)")
    ax.set_ylabel("Proportion")
    ax.set_title("Cellular composition by plaque distance")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, ncol=1)
    plt.tight_layout()
    plt.show()
