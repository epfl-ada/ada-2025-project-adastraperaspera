from __future__ import annotations

from collections.abc import Sequence
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

Number = Union[int, float, np.number]

from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPoly
import numpy as np
from shapely.geometry import MultiPolygon, Polygon


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
        Line2D([0], [0], color="orange", lw=1.0, linestyle=":", label="Convex hull (sample)"),
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
