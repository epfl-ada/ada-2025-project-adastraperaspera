from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import math
from pathlib import Path
from typing import Any, Union

import geopandas as gpd
import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPoly
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from numpy.typing import NDArray
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
from shapely.geometry import MultiPolygon, Polygon
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch

from src.utils.logging_utils import logger

from ..preprocessing.partition import gaussian_kde

Number = Union[int, float, np.number]


def plot_leiden_logit_slopes(
    logit_df: pd.DataFrame,
    my_label_to_type: Mapping[int, str],
    *,
    pval_col: str = "adj_pval",
    slope_col: str = "slope",
    cluster_col: str = "cluster",
    pval_threshold: float = 0.01,
    blue: str = "#4C78A8",
    red: str = "#E45756",
    gap: int = 1,
    fig_width_min: float = 10.0,
    fig_width_per_bar: float = 0.35,
    fig_height: float = 15.0,
    ypad_scale: float = 0.001,
    rotation: int = 90,
    ax: plt.Axes | None = None,
    title: str | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """
    Filter significant entries and plot slopes as a split bar chart:
      - left half: negative slopes (ascending), colored blue, labels at bottom
      - right half: positive slopes (ascending), colored red, labels at top

    Returns (fig, ax, sig_df) where sig_df is the filtered significant DataFrame
    with an added 'cell_type' column.

    Notes:
      - Entries with slope == 0 are excluded by default (neither neg nor pos).
      - Cluster ids are mapped via my_label_to_type[int(cluster_id)].
        If a key is missing, falls back to str(cluster_id).
    """
    required = {pval_col, slope_col, cluster_col}
    missing = required - set(logit_df.columns)
    if missing:
        raise KeyError(f"logit_df is missing required columns: {sorted(missing)}")

    # 0) Sort by adjusted p-value (ascending)
    df = logit_df.sort_values(pval_col, ascending=True).reset_index(drop=True)

    # 1) Keep only significant entries
    sig = df.loc[df[pval_col] < pval_threshold].copy()

    # 3) Map cluster id -> cell type label
    def _map_cluster(x: Any) -> str:
        try:
            return my_label_to_type[int(x)]
        except Exception:
            return str(x)

    sig["cell_type"] = sig[cluster_col].apply(_map_cluster)

    # 2) Split into negative/positive slopes; sort each side ascending by slope
    neg = sig.loc[sig[slope_col] < 0].sort_values(slope_col, ascending=True)
    pos = sig.loc[sig[slope_col] > 0].sort_values(slope_col, ascending=True)

    nL, nR = len(neg), len(pos)
    xL = np.arange(nL)
    xR = np.arange(nR) + nL + gap

    if ax is None:
        fig_w = max(fig_width_min, fig_width_per_bar * (nL + nR + gap))
        fig, ax = plt.subplots(figsize=(fig_w, fig_height))
    else:
        fig = ax.figure

    ax.bar(xL, neg[slope_col].to_numpy(), color=blue)
    ax.bar(xR, pos[slope_col].to_numpy(), color=red)

    # Zero line for reference
    ax.axhline(0, linewidth=1)

    # Clean x-axis (we write cell types directly on bars)
    ax.set_xticks([])
    ax.set_ylabel(slope_col)

    if title:
        ax.set_title(title)

    # Annotation offsets
    max_abs = float(sig[slope_col].abs().max()) if len(sig) else 1.0
    ypad = ypad_scale * max_abs if max_abs > 0 else 0.1

    # Bottom of left (negative) bars
    for x, y, txt in zip(xL, neg[slope_col].to_numpy(), neg["cell_type"].to_numpy(), strict=False):
        ax.text(x, y - ypad, str(txt), ha="center", va="top", rotation=rotation)

    # Top of right (positive) bars
    for x, y, txt in zip(xR, pos[slope_col].to_numpy(), pos["cell_type"].to_numpy(), strict=False):
        ax.text(x, y + ypad, str(txt), ha="center", va="bottom", rotation=rotation)

    ax.margins(x=0.01)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax, sig


def plot_half_distance(
    df: pd.DataFrame,
    gene_col: str = "gene",
    slope_col: str = "slope",
    qval_col: str | None = "qval",
    filter_negative: bool = True,
    sort: str = "half",  # "half" | "abs_half_desc" | "qval_then_half"
    top_n: int | None = None,
    log_scale: bool = False,
    tissue_radius_um: float | None = None,
    annotate: bool = True,
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] | None = None,
):
    """
    Compute and plot distance to halve expression: d_{1/2} = ln(2) / |slope| (µm).

    Parameters
    ----------
    df : DataFrame
        Must contain `gene_col`, `slope_col`, and optionally `qval_col`.
    filter_negative : bool
        If True, keep only rows with slope < 0.
    sort : str
        "half" -> ascending d_{1/2} (shortest first).
        "abs_half_desc" -> descending |d_{1/2}| (longest first).
        "qval_then_half" -> ascending qval then ascending d_{1/2}. Requires qval_col.
    top_n : int | None
        If set, plot only the first N after sorting.
    log_scale : bool
        If True, use a log x-axis.
    tissue_radius_um : float | None
        If set, draw a vertical reference line at this radius (µm).
    annotate : bool
        If True, write the numeric value next to each bar.
    ax : matplotlib Axes | None
        If provided, draw on this axes; otherwise create a new figure/axes.
    figsize : (w, h) | None
        Figure size. If None, height is scaled by number of bars.

    Returns
    -------
    ax : matplotlib Axes
    plot_df : DataFrame
        The dataframe (indexed by gene) with a 'half_dist_um' column used for plotting.
    """
    if filter_negative:
        work = df.loc[
            df[slope_col] < 0, [gene_col, slope_col] + ([qval_col] if qval_col else [])
        ].copy()
    else:
        work = df[[gene_col, slope_col] + ([qval_col] if qval_col else [])].copy()

    # Compute half-distance (µm) and clean
    work["half_dist_um"] = np.log(2) / work[slope_col].abs()
    work.replace([np.inf, -np.inf], np.nan, inplace=True)
    work.dropna(subset=["half_dist_um"], inplace=True)

    # Sorting
    if sort == "half":
        work.sort_values("half_dist_um", inplace=True, ascending=True)
    elif sort == "abs_half_desc":
        work["abs_half"] = work["half_dist_um"].abs()
        work.sort_values("abs_half", inplace=True, ascending=False)
    elif sort == "qval_then_half":
        if qval_col is None or qval_col not in work.columns:
            raise ValueError(
                "qval_then_half sorting requires qval_col to be provided and present in df."
            )
        work.sort_values([qval_col, "half_dist_um"], inplace=True, ascending=[True, True])
    else:
        raise ValueError("sort must be one of {'half','abs_half_desc','qval_then_half'}")

    if top_n is not None:
        work = work.head(int(top_n))

    plot_df = work.set_index(gene_col)

    # Create axes
    n = len(plot_df)
    if ax is None:
        if figsize is None:
            figsize = (11, max(4, 0.28 * n))
        _, ax = plt.subplots(figsize=figsize)

    # Plot
    plot_df["half_dist_um"].plot(kind="barh", ax=ax)
    ax.set_xlabel("Distance to halve expression (µm)")
    ax.set_ylabel("Gene")
    ax.set_title(f"Half-distance d₁/₂ = ln(2)/|β₁| (n={n})")
    ax.invert_yaxis()  # shortest at top if ascending sort

    if log_scale:
        ax.set_xscale("log")
        pad_factor = 1.03
    else:
        pad = (plot_df["half_dist_um"].max() or 0) * 0.01

    # Reference line for tissue radius
    if tissue_radius_um is not None:
        ax.axvline(tissue_radius_um, linestyle="--", linewidth=1)
        ax.text(
            tissue_radius_um,
            -0.6,
            f"radius = {tissue_radius_um:,.0f} µm",
            rotation=0,
            ha="right",
            va="bottom",
        )

    # Annotations
    if annotate:
        for i, v in enumerate(plot_df["half_dist_um"].values):
            if np.isnan(v):
                continue
            if log_scale:
                ax.text(v * pad_factor, i, f"{v:,.0f} µm", va="center", ha="left", fontsize=9)
            else:
                ax.text(v + pad, i, f"{v:,.0f} µm", va="center", ha="left", fontsize=9)

    plt.tight_layout()
    return ax, plot_df


def load_png_rgb(path: Path) -> np.ndarray:
    """Load PNG as RGB array (H, W, 3)."""
    if not path.exists():
        raise FileNotFoundError(str(path))

    img = plt.imread(str(path))

    # Grayscale -> RGB
    if img.ndim == 2:
        img = np.repeat(img[..., None], 3, axis=2)

    # RGBA -> RGB
    if img.ndim == 3 and img.shape[2] == 4:
        img = img[..., :3]

    if img.ndim != 3 or img.shape[2] != 3:
        raise ValueError(f"Unexpected image shape {img.shape} for {path}")

    return img


def pad_to_max(img: np.ndarray, target_h: int, target_w: int) -> np.ndarray:
    """Center-pad with black margins up to (target_h, target_w)."""
    h, w, c = img.shape
    if c != 3:
        raise ValueError("Expected RGB (H,W,3)")

    pad_h = max(0, target_h - h)
    pad_w = max(0, target_w - w)

    top = pad_h // 2
    bottom = pad_h - top
    left = pad_w // 2
    right = pad_w - left

    return np.pad(
        img,
        pad_width=((top, bottom), (left, right), (0, 0)),
        mode="constant",
        constant_values=0,
    )


def _infer_grid_shape(key_to_pos: dict[str, tuple[int, int]]) -> tuple[int, int]:
    """Infer (n_rows, n_cols) from the maximum indices in key_to_pos."""
    max_row = max(r for r, _ in key_to_pos.values())
    max_col = max(c for _, c in key_to_pos.values())
    return max_row + 1, max_col + 1


def make_image_grid(
    *,
    files: dict[str, Path],
    key_to_pos: dict[str, tuple[int, int]],
    col_ticks: Iterable[str],
    row_ticks: Iterable[str],
    flip_h_keys: Iterable[str] = (),
    figsize: tuple[float, float] = (15, 9),
) -> None:
    """
    Create and display an image grid with unified axes.

    Notes
    -----
    - Each image is vertically cropped to remove the bottom 10% strip
      before padding and plotting.
    """
    flip_h_keys = set(flip_h_keys)

    # 1) Load images and crop bottom 10%
    imgs: dict[str, np.ndarray] = {}
    shapes = []

    for key, path in files.items():
        img = load_png_rgb(path)
        if key in flip_h_keys:
            img = np.fliplr(img)

        # --- Crop bottom 10% strip ---
        h = img.shape[0]
        keep_h = max(1, int(round(h * 0.9)))  # keep the top 90%
        img = img[:keep_h, ...]
        # --------------------------------

        imgs[key] = img
        shapes.append(img.shape[:2])

    # 2) Target size (based on cropped images)
    max_h = max(h for h, _ in shapes)
    max_w = max(w for _, w in shapes)

    # 3) Grid creation
    n_rows, n_cols = _infer_grid_shape(key_to_pos)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)

    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1 or n_cols == 1:
        axes = np.atleast_2d(axes)

    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    # 4) Plot images (after padding to a common size)
    for key, (r, c) in key_to_pos.items():
        img = pad_to_max(imgs[key], target_h=max_h, target_w=max_w)
        axes[r, c].imshow(img, aspect="equal", interpolation="nearest")

    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.98, wspace=0.02, hspace=0.02)

    # 5) Unified axis overlay
    positions = [ax.get_position() for ax in axes.ravel()]
    left = min(p.x0 for p in positions)
    right = max(p.x1 for p in positions)
    bottom = min(p.y0 for p in positions)
    top = max(p.y1 for p in positions)

    big_ax = fig.add_axes([left, bottom, right - left, top - bottom], frameon=False)
    big_ax.set_zorder(10)
    big_ax.patch.set_alpha(0)

    big_ax.set_xlim(0, n_cols)
    big_ax.set_ylim(0, n_rows)

    big_ax.set_xticks([i + 0.5 for i in range(n_cols)], labels=list(col_ticks))
    big_ax.set_yticks([n_rows - (i + 0.5) for i in range(n_rows)], labels=list(row_ticks))

    big_ax.set_xlabel("Age (months)")
    big_ax.set_ylabel("Type")
    big_ax.tick_params(axis="both", which="both", length=0)

    plt.show()


def plot_kde_with_stats(
    data: np.ndarray | pd.Series,
    title: str | None = None,
    xlabel: str | None = None,
    color: str | None = None,
    show: bool = True,
) -> dict:
    """
    Plots the KDE of a numeric Series or array with annotated vertical lines
    for min, 25%, 50%, 75%, max, and mean values (each with a unique color).

    Parameters
    ----------
    data : array-like or pd.Series
        The numeric data to plot.
    title : str, optional
        Plot title.
    xlabel : str, optional
        Label for the x-axis.
    color : str, optional
        Color for the KDE curve.
    show : bool, default=True
        Whether to display the plot immediately.

    Returns
    -------
    stats : dict
        Dictionary of computed statistics: min, q25, q50, q75, max, mean.
    """
    stats = {
        "min": data.min(),
        "q25": data.quantile(0.25),
        "q50": data.quantile(0.5),
        "q75": data.quantile(0.75),
        "max": data.max(),
        "mean": data.mean(),
    }

    # Assign distinct colors for each statistic
    stat_colors = {
        "min": "tab:red",
        "q25": "tab:orange",
        "q50": "tab:green",
        "q75": "tab:blue",
        "max": "tab:purple",
        "mean": "tab:brown",
    }

    plt.figure(figsize=(8, 5))
    sns.kdeplot(data, fill=True, color=color or "tab:gray", alpha=0.5, label="KDE")

    # Plot vertical lines with unique colors and styles
    for key, value in stats.items():
        linestyle = "--" if "q" in key else "-."
        plt.axvline(
            value,
            color=stat_colors[key],
            linestyle=linestyle,
            linewidth=1.5,
            label=f"{key}: {value:.3f}",
        )

    plt.legend()
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    plt.ylabel("Density")
    plt.grid(alpha=0.3)

    if show:
        plt.show()

    return stats


def set_pub_style():
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "figure.autolayout": False,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.linewidth": 0.4,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    sns.set_style("whitegrid")


def savefig(path: str, fig=None):
    (fig or plt.gcf()).savefig(path, bbox_inches="tight", pad_inches=0.02)


def plot_gene_distributions(
    df: pd.DataFrame,
    genes: Sequence[str],
    all_genes: Sequence[str] | None = None,
    bins: int = 50,
    cols: int = 4,
    use_log1p: bool = True,
    kde: bool = True,
    clip_quantiles: tuple[float, float] = (0.0, 0.995),
    figsize: tuple[int, int] | None = None,
    suptitle: str | None = None,
    show_zero_fraction: bool = True,
    show_iqr: bool = True,
    show_median: bool = True,
    hist_alpha: float = 0.6,
    kde_lw: float = 1.6,
    dpi: int = 120,
    save_path: str | None = None,
):
    arrays = []
    present = [g for g in genes if g in df.columns]
    if not present:
        raise ValueError("None of the requested genes are present in the dataframe.")

    if all_genes is None:
        all_genes = genes
    all_present = [g for g in all_genes if g in df.columns]

    avg_arrays = []
    for g in all_present:
        x = df[g].to_numpy(dtype=float)
        if use_log1p:
            x = np.log1p(x)
        x = x[np.isfinite(x)]
        if x.size:
            avg_arrays.append(x)

    for g in present:
        x = df[g].to_numpy(dtype=float)
        if use_log1p:
            x = np.log1p(x)
        x = x[np.isfinite(x)]
        arrays.append(x)

    all_vals = (
        np.concatenate(arrays + avg_arrays) if (len(arrays) + len(avg_arrays)) else np.array([])
    )
    if all_vals.size == 0:
        raise ValueError("No finite values to plot.")
    lo = np.quantile(all_vals, clip_quantiles[0])
    hi = np.quantile(all_vals, clip_quantiles[1])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))

    bin_edges = np.linspace(lo, hi, bins + 1)
    bw = float(bin_edges[1] - bin_edges[0]) if bins > 0 else 1.0
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    gene_hist_color = "#4C78A8"
    gene_kde_color = "#1F77B4"
    gene_iqr_color = "#B3D9FF"

    avg_hist_color = "#F4A261"
    avg_kde_color = "#F28E2B"
    avg_iqr_color = "#FFE0B2"

    iqr_alpha = 0.45

    avg_hist = None
    avg_xs = None
    avg_kde_counts = None
    avg_q25 = avg_q75 = None
    avg_mean = None
    avg_zero_frac = None

    if avg_arrays:
        hists = []
        q25s = []
        q75s = []
        means = []
        for x in avg_arrays:
            h, _ = np.histogram(x, bins=bin_edges, density=False)
            hists.append(h.astype(float))
            q25, q75 = np.quantile(x, [0.25, 0.75])
            q25s.append(float(q25))
            q75s.append(float(q75))
            means.append(float(np.mean(x)))
        avg_hist = np.mean(np.stack(hists, axis=0), axis=0) if hists else None
        avg_q25 = float(np.mean(q25s)) if q25s else None
        avg_q75 = float(np.mean(q75s)) if q75s else None
        avg_mean = float(np.mean(means)) if means else None

        if kde:
            avg_xs = np.linspace(lo, hi, 400)
            pdfs_counts = []
            for x in avg_arrays:
                if x.size > 5:
                    x_pos = x[x > 0]
                    x_kde = x_pos if x_pos.size > 5 else x
                    try:
                        pdf = gaussian_kde(x_kde)(avg_xs)
                        pdfs_counts.append(pdf * float(x.size) * bw)
                    except Exception:
                        pass
            if pdfs_counts:
                avg_kde_counts = np.mean(np.stack(pdfs_counts, axis=0), axis=0)

        if show_zero_fraction and all_present:
            zfs = []
            for g in all_present:
                raw = df[g].to_numpy(dtype=float)
                if use_log1p:
                    zfs.append(float(np.mean(raw <= 0.0)))
                else:
                    zfs.append(float(np.mean(raw == 0.0)))
            avg_zero_frac = float(np.mean(zfs)) if zfs else None

    y_max_global = 1.0
    for x in arrays:
        if x.size:
            y_max_global = max(
                y_max_global, float(np.max(np.histogram(x, bins=bin_edges, density=False)[0]))
            )
    if avg_hist is not None and np.isfinite(avg_hist).any():
        y_max_global = max(y_max_global, float(np.nanmax(avg_hist)))
    if avg_kde_counts is not None and np.isfinite(avg_kde_counts).any():
        y_max_global = max(y_max_global, float(np.nanmax(avg_kde_counts)))

    n = len(present)
    rows = (n + cols - 1) // cols
    if figsize is None:
        figsize = (4.2 * cols, 3.2 * rows)

    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.linewidth": 0.4,
            "grid.alpha": 0.35,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    fig, axes = plt.subplots(rows, cols, figsize=figsize, sharex=True, sharey=True)
    fig.set_facecolor("white")
    axes = np.atleast_1d(axes).ravel()

    bar_w = 0.44 * bw
    shift = 0.23 * bw
    bar_alpha = 1.0

    for i, (ax, gene) in enumerate(zip(axes, present, strict=False)):
        ax.set_facecolor("#FCFCFF")

        x_raw = df[gene].to_numpy(dtype=float)
        x = x_raw.copy()
        if use_log1p:
            x = np.log1p(x)
        x = x[np.isfinite(x)]
        if x.size == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", color="black")
            ax.axis("off")
            continue

        if show_iqr and avg_q25 is not None and avg_q75 is not None:
            ax.axvspan(
                avg_q25,
                avg_q75,
                facecolor=avg_iqr_color,
                edgecolor=avg_iqr_color,
                alpha=iqr_alpha,
                linewidth=0.6,
                zorder=0.6,
            )

        q25 = q75 = None
        if show_iqr:
            q25, q75 = np.quantile(x, [0.25, 0.75])
            ax.axvspan(
                float(q25),
                float(q75),
                facecolor=gene_iqr_color,
                edgecolor=gene_iqr_color,
                alpha=iqr_alpha,
                linewidth=0.6,
                zorder=0.7,
            )

        gene_counts, _ = np.histogram(x, bins=bin_edges, density=False)

        if avg_hist is not None:
            ax.bar(
                centers - shift,
                avg_hist,
                width=bar_w,
                align="center",
                alpha=bar_alpha,
                color=avg_hist_color,
                edgecolor=avg_hist_color,
                linewidth=0.35,
                zorder=1.0,
            )

        ax.bar(
            centers + shift,
            gene_counts,
            width=bar_w,
            align="center",
            alpha=bar_alpha,
            color=gene_hist_color,
            edgecolor=gene_hist_color,
            linewidth=0.35,
            zorder=1.1,
        )

        gene_pdf_counts = None
        if kde and x.size > 5:
            x_pos = x[x > 0]
            x_kde = x_pos if x_pos.size > 5 else x
            try:
                xs = np.linspace(lo, hi, 400)
                gene_pdf_counts = gaussian_kde(x_kde)(xs) * float(x.size) * bw
                ax.plot(xs, gene_pdf_counts, lw=kde_lw, color=gene_kde_color, zorder=2.2)
            except Exception:
                gene_pdf_counts = None

        if avg_kde_counts is not None and avg_xs is not None:
            ax.plot(avg_xs, avg_kde_counts, lw=kde_lw, color=avg_kde_color, zorder=2.1)

        if show_median:
            gene_mean = float(np.mean(x))
            ax.axvline(gene_mean, color=gene_kde_color, ls="--", lw=1.2, zorder=2.35)
            if avg_mean is not None:
                ax.axvline(avg_mean, color=avg_kde_color, ls="--", lw=1.2, zorder=2.3)

        if show_zero_fraction:
            gene_zero_frac = None
            if use_log1p:
                gene_zero_frac = float(np.mean(x_raw <= 0.0))
            else:
                gene_zero_frac = float(np.mean(x_raw == 0.0))

            nz_gene = np.flatnonzero(gene_counts > 0)
            if nz_gene.size:
                jg = int(nz_gene[0])
                ax.annotate(
                    f"{gene_zero_frac*100:.1f}%",
                    xy=(centers[jg] + shift, float(gene_counts[jg])),
                    xytext=(8, 3),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                    fontsize=8.5,
                    color=gene_hist_color,
                    zorder=3.0,
                )

            if avg_hist is not None and avg_zero_frac is not None:
                nz_avg = np.flatnonzero(avg_hist > 0)
                if nz_avg.size:
                    ja = int(nz_avg[0])
                    dy = 3 if (nz_gene.size and int(nz_gene[0]) != ja) else 10
                    ax.annotate(
                        f"{avg_zero_frac*100:.1f}%",
                        xy=(centers[ja] - shift, float(avg_hist[ja])),
                        xytext=(-8, dy),
                        textcoords="offset points",
                        ha="right",
                        va="bottom",
                        fontsize=8.5,
                        color=avg_hist_color,
                        zorder=3.0,
                    )

        ax.set_xlim(lo, hi)
        ax.set_title(gene)
        ax.grid(True, axis="y")
        ax.grid(False, axis="x")
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))

        r = i // cols
        c = i % cols
        if r < rows - 1:
            ax.tick_params(labelbottom=False)
        if c > 0:
            ax.tick_params(labelleft=False)

        y_local = float(np.max(gene_counts)) if gene_counts.size else 1.0
        if avg_hist is not None and np.isfinite(avg_hist).any():
            y_local = max(y_local, float(np.nanmax(avg_hist)))
        if gene_pdf_counts is not None and np.isfinite(gene_pdf_counts).any():
            y_local = max(y_local, float(np.nanmax(gene_pdf_counts)))
        if avg_kde_counts is not None and np.isfinite(avg_kde_counts).any():
            y_local = max(y_local, float(np.nanmax(avg_kde_counts)))
        y_max_global = max(y_max_global, y_local)

    for k in range(n, len(axes)):
        axes[k].axis("off")

    for ax in axes[:n]:
        ax.set_ylim(0, y_max_global * 1.10)

    handles = []
    labels = []

    if avg_hist is not None:
        handles.append(Patch(facecolor=avg_hist_color, edgecolor=avg_hist_color, alpha=1.0))
        labels.append("Hist (avg)")

    handles.append(Patch(facecolor=gene_hist_color, edgecolor=gene_hist_color, alpha=1.0))
    labels.append("Hist (gene)")

    if kde and avg_kde_counts is not None:
        handles.append(Line2D([0], [0], color=avg_kde_color, lw=kde_lw))
        labels.append("KDE (avg)")

    if kde:
        handles.append(Line2D([0], [0], color=gene_kde_color, lw=kde_lw))
        labels.append("KDE (gene)")

    if show_median and avg_mean is not None:
        handles.append(Line2D([0], [0], color=avg_kde_color, lw=1.2, ls="--"))
        labels.append("Mean (avg)")

    if show_median:
        handles.append(Line2D([0], [0], color=gene_kde_color, lw=1.2, ls="--"))
        labels.append("Mean (gene)")

    if show_iqr and avg_q25 is not None and avg_q75 is not None:
        handles.append(Patch(facecolor=avg_iqr_color, edgecolor=avg_iqr_color, alpha=iqr_alpha))
        labels.append("IQR (avg)")

    if show_iqr:
        handles.append(Patch(facecolor=gene_iqr_color, edgecolor=gene_iqr_color, alpha=iqr_alpha))
        labels.append("IQR (gene)")

    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
            ncol=4,
            frameon=True,
            fancybox=True,
            framealpha=0.92,
            edgecolor="0.85",
            fontsize=9,
            handlelength=2.2,
            handletextpad=0.6,
            columnspacing=1.2,
        )

    if suptitle:
        fig.suptitle(suptitle, y=0.998, fontsize=12)

    fig.supxlabel("log1p(counts)" if use_log1p else "counts", y=0.03, fontsize=11)
    fig.supylabel("Count", x=0.01, fontsize=11)

    plt.tight_layout(rect=[0.02, 0.05, 1, 0.94])

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    plt.show()
    return fig


def plot_weird_gene_panels(
    df: pd.DataFrame,
    genes: Sequence[str],
    *,
    bins: int = 50,
    cols: int = 4,
    # aesthetics / behavior aligned with plot_gene_distributions:
    clip_quantiles: tuple[float, float] = (
        0.0,
        1.0,
    ),  # (0, 1) = full range; set to (0, 0.995) to trim tails
    kde: bool = True,
    show_zero_fraction: bool = True,
    show_iqr: bool = True,
    show_median: bool = True,
    hist_alpha: float = 0.5,
    kde_lw: float = 1.6,
    dpi: int = 120,
    linear_title: str = "Weirdest genes — linear scale (counts)",
    log_title: str = "Weirdest genes — log1p scale",
    figsize: tuple[int, int] | None = None,
    show: bool = True,
    save_linear_path: str | None = None,
    save_log_path: str | None = None,
):
    """
    Publication-ready panels for 'weird' genes using log1p(expression)
    (histogram + KDE + median + IQR + zero %).
    """

    present = [g for g in genes if g in df.columns]
    if not present:
        raise ValueError("None of the requested genes are present in the dataframe.")

    # ---- gather arrays for global axis limits (linear and log1p separately)
    X_lin = []
    X_log = []
    for g in present:
        x = df[g].to_numpy(dtype=float)
        x = x[np.isfinite(x) & (x >= 0)]
        if x.size:
            X_lin.append(x)
            X_log.append(np.log1p(x))
    if not X_lin:
        raise ValueError("No finite, non-negative values to plot.")

    all_lin = np.concatenate(X_lin)
    all_log = np.concatenate(X_log)

    def _limits(arr, q=(0.0, 1.0)):
        lo = np.quantile(arr, q[0])
        hi = np.quantile(arr, q[1])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
        # ensure non-degenerate limits
        if hi == lo:
            hi = lo + 1.0
        return lo, hi

    lo_lin, hi_lin = _limits(all_lin, clip_quantiles)
    lo_log, hi_log = _limits(all_log, clip_quantiles)

    # pre-compute shared bin edges
    bins_lin = np.linspace(lo_lin, hi_lin, bins + 1)
    bins_log = np.linspace(lo_log, hi_log, bins + 1)

    # ---- shared rcParams for both figures
    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.linewidth": 0.4,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    n = len(present)
    rows = (n + cols - 1) // cols
    if figsize is None:
        figsize = (4.2 * cols, 3.2 * rows)

    # ---------------- Linear (counts) ----------------
    fig_lin, axes_lin = plt.subplots(rows, cols, figsize=figsize)
    axes_lin = np.atleast_1d(axes_lin).ravel()

    for ax, g in zip(axes_lin, present, strict=False):
        x = df[g].to_numpy(dtype=float)
        x = x[np.isfinite(x) & (x >= 0)]
        if x.size == 0:
            ax.axis("off")
            continue

        # hist
        counts, _, _ = ax.hist(
            x, bins=bins_lin, density=True, alpha=hist_alpha, color="#4C78A8", edgecolor="none"
        )

        # kde (positive support preferred)
        if kde and x.size > 5:
            x_pos = x[x > 0]
            x_kde = x_pos if x_pos.size > 5 else x
            try:
                xs = np.linspace(lo_lin, hi_lin, 400)
                pdf = gaussian_kde(x_kde)(xs)
                ax.plot(xs, pdf, lw=kde_lw, color="#1F77B4")
            except Exception:
                pass

        # median & IQR
        if show_median or show_iqr:
            q25, q50, q75 = np.quantile(x, [0.25, 0.50, 0.75])
            if show_iqr:
                ax.axvspan(q25, q75, color="0.85", alpha=0.6, zorder=0, label="IQR")
            if show_median:
                ax.axvline(q50, color="#E45756", ls="--", lw=1.2, label="Median")

        # zeros
        if show_zero_fraction:
            zf = (x == 0).mean() * 100.0
            ax.text(
                0.98,
                0.95,
                f"zeros: {zf:.1f}%",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8,
                color="0.3",
            )

        ax.set_xlim(lo_lin, hi_lin)
        ax.set_title(g)
        ax.set_xlabel("Expression (counts)")
        ax.set_ylabel("Density")

    for k in range(n, len(axes_lin)):
        axes_lin[k].axis("off")

    fig_lin.suptitle(linear_title, y=0.995, fontsize=12)
    fig_lin.tight_layout(rect=[0, 0, 1, 0.98])
    if save_linear_path:
        fig_lin.savefig(save_linear_path, bbox_inches="tight")

    # ---------------- Log1p ----------------
    fig_log, axes_log = plt.subplots(rows, cols, figsize=figsize)
    axes_log = np.atleast_1d(axes_log).ravel()

    for ax, g in zip(axes_log, present, strict=False):
        x = df[g].to_numpy(dtype=float)
        x = x[np.isfinite(x) & (x >= 0)]
        if x.size == 0:
            ax.axis("off")
            continue
        xlog = np.log1p(x)

        counts, _, _ = ax.hist(
            xlog, bins=bins_log, density=True, alpha=hist_alpha, color="#4C78A8", edgecolor="none"
        )

        if kde and xlog.size > 5:
            x_pos = xlog[x > 0]  # keep positives based on raw
            x_kde = x_pos if x_pos.size > 5 else xlog
            try:
                xs = np.linspace(lo_log, hi_log, 400)
                pdf = gaussian_kde(x_kde)(xs)
                ax.plot(xs, pdf, lw=kde_lw, color="#1F77B4")
            except Exception:
                pass

        if show_median or show_iqr:
            q25, q50, q75 = np.quantile(xlog, [0.25, 0.50, 0.75])
            if show_iqr:
                ax.axvspan(q25, q75, color="0.85", alpha=0.6, zorder=0)
            if show_median:
                ax.axvline(q50, color="#E45756", ls="--", lw=1.2)

        if show_zero_fraction:
            zf = (x == 0).mean() * 100.0
            ax.text(
                0.98,
                0.95,
                f"zeros: {zf:.1f}%",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8,
                color="0.3",
            )

        ax.set_xlim(lo_log, hi_log)
        ax.set_title(g)
        ax.set_xlabel("Expression (log1p counts)")
        ax.set_ylabel("Density")

    for k in range(n, len(axes_log)):
        axes_log[k].axis("off")

    fig_log.suptitle(log_title, y=0.995, fontsize=12)
    fig_log.tight_layout(rect=[0, 0, 1, 0.98])
    if save_log_path:
        fig_log.savefig(save_log_path, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig_lin)
        plt.close(fig_log)

    return fig_lin, fig_log


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
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.show()


def plot_spatial_overlay(df, gene, alpha=0.7, sample_size=20000):
    """
    Scatter overlay showing per-cell expression intensity in tissue coordinates.
    Colors high expression regions with a semi-transparent viridis map.
    """
    import matplotlib.pyplot as plt

    data = df.sample(min(sample_size, len(df)), random_state=42)
    # Plot raw expression
    # Undo the log1p transformation
    data[gene] = np.expm1(data[gene])
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
    plt.colorbar(sc, label="Expression")
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
    cbar = plt.colorbar(sc, label="Expression (log1p)")
    plt.tight_layout()
    plt.show()


def plot_spatial_with_plaques(
    df, plaques_poly, gene=None, score_col=None, sample_size=20000, logger=None
):
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
        if logger is not None:
            logger.warning(f"Could not overlay plaques: {e}")
        else:
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
    logger=None,
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
    fig.colorbar(sc, ax=axes[0], label=f"{oligo_marker} expression (log1p)")

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
            if logger is not None:
                logger.warning(f"Plaque overlay skipped: {e}")
            else:
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
    # --- new rotation controls ---
    rotate_180: bool = True,
    rotation_origin: tuple[float, float] | str = "auto",
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
    rotate_180 : bool, optional
        If True, rotate all geometries by 180° clockwise about `rotation_origin`.
        (Note: 180° CW == 180° CCW.)
    rotation_origin : tuple(x, y) | "auto", optional
        Rotation pivot. If "auto", uses the center of the combined bounds of all
        geometries (plaques and brain ROI). If a (x, y) tuple, rotates around it.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    from shapely.affinity import rotate as shp_rotate

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

    # ---- rotation prep ----------------------------------------------------
    def _compute_auto_origin() -> tuple[float, float]:
        """Center of combined bounds of plaques and brain geometry if present."""
        bounds = []
        for g in P["geometry"]:
            if g is not None and hasattr(g, "bounds"):
                try:
                    b = g.bounds  # (minx, miny, maxx, maxy)
                    if b and np.isfinite(b).all():
                        bounds.append(b)
                except Exception:
                    pass
        if brain_geom is not None and getattr(brain_geom, "is_valid", False):
            try:
                b = brain_geom.bounds
                if b and np.isfinite(b).all():
                    bounds.append(b)
            except Exception:
                pass
        if not bounds:
            return (0.0, 0.0)
        arr = np.array(bounds, dtype=float)
        minx, miny = arr[:, 0].min(), arr[:, 1].min()
        maxx, maxy = arr[:, 2].max(), arr[:, 3].max()
        return ((minx + maxx) / 2.0, (miny + maxy) / 2.0)

    if isinstance(rotation_origin, str):
        if rotation_origin != "auto":
            raise ValueError("rotation_origin must be a (x, y) tuple or 'auto'.")
        origin_point = _compute_auto_origin()
    else:
        origin_point = rotation_origin

    def _rot(g):
        """Rotate geometry by 180° about origin if enabled."""
        if not rotate_180 or g is None:
            return g
        # 180° CW == 180° CCW; shapely rotates CCW by default.
        return shp_rotate(g, 180.0, origin=origin_point, use_radians=False)

    # ---- draw brain ROI (optional) ----------------------------------------
    if brain_geom is not None and getattr(brain_geom, "is_valid", False):
        g_b = _rot(brain_geom)
        bx, by = g_b.exterior.xy
        ax.plot(bx, by, color="blue", lw=1.0, label="Brain ROI")

    # ---- draw plaques by convexity ----------------------------------------
    for row in P.itertuples(index=False):
        g = getattr(row, "geometry", None)
        is_convex = getattr(row, "is_convex", None)
        if g is None:
            continue
        g_draw = _rot(g)

        if bool(is_convex):
            _add_geom(
                g_draw,
                facecolor="none",
                edgecolor="green",
                linestyle="-",
                linewidth=0.8,
            )
        else:
            _add_geom(
                g_draw,
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
            hull = g.convex_hull
            _add_geom(
                _rot(hull),
                facecolor="none",
                edgecolor="orange",
                linestyle=":",
                linewidth=1.0,
            )

    # ---- cosmetics ---------------------------------------------------------
    ax.set_aspect("equal", "box")
    title = "Plaque geometries after normalization"
    if rotate_180:
        title += " (rotated 180°)"
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("X coordinate (µm)")
    ax.set_ylabel("Y coordinate (µm)")

    legend_elems = [
        Line2D([0], [0], color="blue", lw=1.0, label="Brain ROI"),
        Line2D([0], [0], color="green", lw=0.8, linestyle="-", label="Convex"),
        Line2D([0], [0], color="red", lw=0.8, linestyle="--", label="Non-convex"),
        Line2D([0], [0], color="orange", lw=1.0, linestyle=":", label="Convex hull (sample)"),
    ]
    handles, labels = [], []
    if brain_geom is not None and getattr(brain_geom, "is_valid", False):
        handles.append(legend_elems[0])
        labels.append(legend_elems[0].get_label())
    handles.extend(legend_elems[1:])
    labels.extend([e.get_label() for e in legend_elems[1:]])
    ax.legend(handles, labels, loc="lower left", frameon=False)

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


def plot_cell_type_comp_by_dist(props):
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


def plot_mean_exp_dist(mean_melt):
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


def plot_pig_expr_by_dist(mean_by_bin, pig_cols):
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


def plot_pig_comp_heatmap(pig_mat, prop_mat, pig_cols):
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

    # Black out non‑significant cells at 0.01 FDR by setting them to NaN
    sig = q <= 0.01
    corrs_masked = corrs.astype(float).where(sig)

    # The colormap then renders NaNs as black
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    try:
        cmap.set_bad("black")
    except Exception:
        pass

    plt.figure(figsize=(min(14, 6 + 0.25 * len(corrs.columns)), 8))
    sns.heatmap(
        corrs_masked,
        cmap=cmap,
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


def plot_cellular_comp_by_dist(pivot_prop):
    ax = pivot_prop.plot(kind="bar", stacked=True, figsize=(8, 5), width=0.85, colormap="tab20")
    ax.set_xlabel("Distance to plaque (µm)")
    ax.set_ylabel("Proportion")
    ax.set_title("Cellular composition by plaque distance")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, ncol=1)
    plt.tight_layout()
    plt.show()


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
    # ---- Basic validation (unchanged) ----
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

    # --- compute a consistent order + rounded labels for the distance bins
    order, labels = _sorted_bins_and_labels(combined_df_normalized[distance_col])

    # ---- Layout (unchanged) ----
    n_rows = max(1, math.ceil(len(top_genes) / n_cols))
    fig_height = n_rows * float(fig_height_per_row)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(float(fig_width), fig_height), sharey=sharey)

    if isinstance(axes, np.ndarray):
        axes_flat = axes.ravel().tolist()
    else:
        axes_flat = [axes]

    # ---- Plot each gene ----
    for i, g in enumerate(top_genes):
        ax = axes_flat[i]

        # Base boxplot with ordered bins
        sns.boxplot(
            data=combined_df_normalized,
            x=distance_col,
            y=g,
            ax=ax,
            order=order,  # NEW
            showfliers=False,
            palette=palette,
        )

        # Overlay nonzero points (sampled), with same order
        nonzero = combined_df_normalized.loc[combined_df_normalized[g] > 0, [distance_col, g]]
        if len(nonzero) > 0 and max_points > 0:
            sample_n = min(int(max_points), len(nonzero))
            sample = nonzero.sample(sample_n, random_state=random_state)
            sns.stripplot(
                data=sample,
                x=distance_col,
                y=g,
                ax=ax,
                order=order,  # NEW
                color="black",
                size=jitter_size,
                alpha=jitter_alpha,
                jitter=jitter,
            )

        ax.set_title(str(g))
        ax.set_xlabel("")
        ax.set_ylabel(y_label)
        ax.set_xticklabels(labels, rotation=45, ha="right")  # NEW

    # Hide any unused axes
    for j in range(len(top_genes), len(axes_flat)):
        axes_flat[j].axis("off")

    fig.suptitle(suptitle, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if show:
        plt.show()

    return fig, axes_flat, top_genes


def _sorted_bins_and_labels(series):
    """
    Return (order, labels) for a distance-bin column that may contain
    Interval categories. Falls back to string labels if not Intervals.
    """
    vals = series.dropna()
    # Keep category order if it's categorical, otherwise use unique values
    cats = (
        list(vals.cat.categories)
        if pd.api.types.is_categorical_dtype(vals)
        else list(pd.unique(vals))
    )

    def _mid(v):
        # Use midpoint for sorting if it looks like an Interval
        try:
            return float(v.mid)
        except Exception:
            return float("inf")

    def _label(v):
        # Pretty rounded bounds if Interval, otherwise str()
        try:
            return f"{int(round(v.left))}-{int(round(v.right))}"
        except Exception:
            return str(v)

    order = sorted(cats, key=_mid)
    labels = [_label(v) for v in order]
    return order, labels


def plot_expression_heatmap(
    summary_df: pd.DataFrame,
    *,
    index_col: str = "gene",
    column_col: str = "distance_bin",
    value_col: str = "mean_expr",
    title: str = "PIG Expression by Distance to Plaque (TgCRND8 17.9m)",
    figsize: tuple[float, float] = (8, 6),
    cmap: str = "mako_r",
    cbar_label: str = "Mean log1+ expression",
    xlabel: str | None = "Distance to plaque (µm, binned)",
    ylabel: str | None = "Gene",
    ax: Axes | None = None,
    show: bool = True,
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
    # Pivot to genes x distance bins
    heatmap_df = summary_df.pivot(index=index_col, columns=column_col, values=value_col)

    # --- Build a robust order + label set for ANY column type ---
    cols_list = list(heatmap_df.columns)

    def _mid(x):
        # Midpoint key for sorting
        try:
            # pd.Interval has .mid; if categorical w/ interval categories, items are Intervals too
            return float(x.mid)
        except Exception:
            return np.inf  # non-intervals (e.g., 'NA') go to the end

    def _label(x):
        # Rounded label for ticks
        try:
            return f"{int(round(x.left))}-{int(round(x.right))}"
        except Exception:
            return str(x)

    order_idx = np.argsort([_mid(c) for c in cols_list])
    cols_sorted = [cols_list[i] for i in order_idx]
    labels = [_label(c) for c in cols_sorted]

    # Reorder columns by midpoint
    heatmap_df = heatmap_df[cols_sorted]

    # --- Plot ---
    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        heatmap_df,
        cmap=cmap,
        cbar_kws={"label": cbar_label},
        # linewidths=0.2,
        # linecolor="white",
    )

    # Apply the rounded bin labels on the x-axis
    ax.set_xticklabels(labels, rotation=35, ha="right")

    # Titles & axes
    ax.set_title(title, pad=8)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

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
    y_label: str = "Mean log-expression (± 95% CI ≈ 1.96×SEM)",
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
    Plot mean gene expression (± 95% CI ≈ 1.96×SEM) against distance-to-plaque bins for a set of genes.

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
        Global y-axis label. Defaults to "Mean log-expression (± 95% CI ≈ 1.96×SEM)".
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
            yerr=(1.96 * df_g[sem_col]),
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


def plot_top_spatial_genes(
    stats_df: pd.DataFrame,
    top_n: int = 20,
    metric: str = "spearman_r",
    figsize: tuple = (6, 6),
    pig_genes: list[str] | None = None,
) -> None:
    """Bar plot of genes most associated with plaque distance.
    Adds '*' to PIGs and ensures consistent matching.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if metric not in stats_df.columns:
        raise ValueError(f"Metric '{metric}' not found in stats_df columns.")

    # Normalize PIG list for safe matching
    if pig_genes is None:
        pig_genes = []
    pig_set = {g.strip().upper() for g in pig_genes}

    # Pick top genes
    top = stats_df.nlargest(top_n, metric).copy()

    # Normalize df gene names for matching
    top["_gene_norm"] = top["gene"].astype(str).str.strip().str.upper()

    # Mark PIGs
    top["is_pig"] = top["_gene_norm"].isin(pig_set)
    top["gene_label"] = top.apply(lambda r: f"{r['gene']} *" if r["is_pig"] else r["gene"], axis=1)

    # Warn if no PIGs marked
    if not top["is_pig"].any() and len(pig_genes) > 0:
        print(
            "None of the top genes match provided PIG list. "
            "Check naming (e.g., symbol vs Ensembl)."
        )

    # Plotting
    order = top["gene_label"]

    plt.figure(figsize=figsize)
    sns.set_style("whitegrid")
    ax = sns.barplot(
        data=top,
        x=metric,
        y="gene_label",
        order=order,
        palette="vlag" if "spearman" in metric.lower() else "crest",
    )

    # Axis labels with units
    if "spearman" in metric.lower():
        xlabel = "Spearman ρ (unitless)"
    elif "slope" in metric.lower():
        xlabel = "OLS slope (Δ log1+ expression per µm)"
    else:
        xlabel = metric

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Gene")
    ax.set_title(f"Top {len(top)} genes by {metric}")

    # Zero reference
    if top[metric].min() < 0 < top[metric].max():
        ax.axvline(0, color="0.4", lw=1, ls="--")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.show()


def plot_cell_to_plaque_map_visible(
    cells_df,
    plaques_gdf,
    x_col="x_centroid",
    y_col="y_centroid",
    dist_col="distance_to_plaque",
    cmap="plasma",
    figsize=(7, 6),
    clip_quantiles=(0.01, 0.99),
    vmin=None,
    vmax=None,
    max_points=None,
    point_size=4,
    point_alpha=0.85,
    plaque_edgecolor="cyan",
    plaque_linewidth=1.0,
    invert_y=False,  # <- changed default (set True if image needs flipping)
    show_scalebar=True,
    scalebar_um=100,
    title="Cell–plaque distance map (µm)",
):
    """
    Spatial map showing distance to plaque.
    """
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import numpy as np

    C = cells_df
    if max_points and len(C) > max_points:
        C = C.sample(n=max_points, random_state=0)

    dvals = C[dist_col].astype(float)
    lo = np.quantile(dvals, clip_quantiles[0]) if vmin is None else vmin
    hi = np.quantile(dvals, clip_quantiles[1]) if vmax is None else vmax
    norm = Normalize(vmin=lo, vmax=hi)

    fig, ax = plt.subplots(figsize=figsize)

    sc = ax.scatter(
        C[x_col],
        C[y_col],
        c=C[dist_col],
        cmap=cmap,
        norm=norm,
        s=point_size,
        alpha=point_alpha,
        edgecolors="none",
        rasterized=True,
    )

    if plaques_gdf is not None and len(plaques_gdf):
        plaques_gdf.plot(
            ax=ax,
            facecolor="none",
            edgecolor=plaque_edgecolor,
            linewidth=plaque_linewidth,
            zorder=10,
        )

    ax.set_aspect("equal")
    if invert_y:
        ax.invert_yaxis()

    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title(title)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Distance to plaque (µm)")
    cbar.ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))

    # Legend outside to avoid covering tissue
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            markersize=6,
            linestyle="None",
            markerfacecolor="gray",
            alpha=0.8,
            label="Cells",
        ),
        Line2D([0], [0], color=plaque_edgecolor, lw=plaque_linewidth, label="Plaque boundary"),
    ]
    ax.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),  # inside border
        frameon=True,
        facecolor="white",
        edgecolor="0.85",
        borderpad=0.3,
        handlelength=1.2,
        handletextpad=0.4,
        fontsize=9,
    )

    # Scalebar
    if show_scalebar:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xb = x0 + 0.05 * (x1 - x0)
        yb = y1 - 0.05 * (y1 - y0)
        ax.plot([xb, xb + scalebar_um], [yb, yb], color="k", lw=2)
        ax.text(
            xb + scalebar_um / 2, yb, f"{int(scalebar_um)} µm", ha="center", va="bottom", fontsize=9
        )

    plt.tight_layout()
    plt.show()
    return fig, ax


def plot_gene_trends(
    mean_expr: pd.DataFrame,
    genes: list[str],
    *,
    sem_expr: pd.DataFrame | None = None,  # <- pass SEM per bin here (optional)
    ci: str = "95ci",  # "95ci" (1.96*SEM) or "sem"
    xlabel: str = "Distance to plaque (µm, binned)",
    ylabel: str = "Mean expression (log1+)",
    title: str = "Spatial gene expression gradients (mean ± CI)",
    figsize: tuple = (8.5, 5),
    legend_loc: str = "best",
    min_visible_err: float | None = None,  # e.g., 0.003 to avoid invisible caps (optional)
) -> None:
    """
    Plot mean expression across distance bins for selected genes,
    with optional SEM-based uncertainty ribbons.

    mean_expr: index = distance_bin (Interval), columns = genes (log1p means)
    sem_expr:  same index/columns as mean_expr (SEM on log1p scale)
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if mean_expr.empty:
        raise ValueError("mean_expr is empty; check your inputs")

    # x-coordinates & tick labels from IntervalIndex
    bins = list(mean_expr.index)
    x = np.arange(len(bins))
    xticklabels = [f"{b.left:.0f}-{b.right:.0f}" for b in bins]

    # scale for CI
    if sem_expr is not None:
        if ci.lower() == "95ci":
            scale = 1.96
        elif ci.lower() == "sem":
            scale = 1.0
        else:
            raise ValueError("ci must be '95ci' or 'sem'")

    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.linewidth": 0.4,
        }
    )

    plt.figure(figsize=figsize)
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    any_plotted = False
    for i, g in enumerate(genes):
        if g not in mean_expr.columns:
            logger.warning(f"[WARN] Gene '{g}' not in mean_expr; skipping.")
            continue

        y = mean_expr[g].to_numpy(dtype=float)
        c = color_cycle[i % len(color_cycle)]

        # shaded ribbon + error bars if SEM provided
        if sem_expr is not None and g in sem_expr.columns:
            e = sem_expr[g].to_numpy(dtype=float) * scale
            if min_visible_err is not None:
                e = np.maximum(e, float(min_visible_err))

            # ribbon
            plt.fill_between(x, y - e, y + e, color=c, alpha=0.18, linewidth=0)
            # line with caps
            plt.errorbar(x, y, yerr=e, color=c, marker="o", lw=1.6, capsize=4, label=g)
        else:
            plt.plot(x, y, color=c, marker="o", lw=1.8, label=g)

        any_plotted = True

    if not any_plotted:
        plt.text(0.5, 0.5, "No valid genes found", ha="center", va="center")
        plt.axis("off")
        plt.show()
        return

    plt.xticks(x, xticklabels, rotation=35, ha="right")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, pad=6)
    plt.legend(title="Genes", frameon=False, loc=legend_loc, ncol=1)
    plt.tight_layout()
    plt.show()


def plot_mean_heatmap(
    mean_expr: pd.DataFrame,
    top_n: int = 25,
    *,
    zscore: bool = True,
    pig_genes: list[str] | None = None,
    title: str | None = None,
    figsize: tuple = (9, 6),
) -> None:
    """
    Visualize top N genes with strongest spatial variation across plaque distance.
    Marks PIG genes with *.
    """
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
    M = mean_expr[gene_cols].copy()

    if isinstance(M.index, pd.IntervalIndex):
        M = M.sort_index(key=lambda x: x.map(lambda i: i.mid))

    grad = M.diff().abs().sum().sort_values(ascending=False)
    top_genes = grad.head(top_n).index
    sub_df = M[top_genes]

    if zscore:
        sub_df = sub_df.apply(lambda s: (s - s.mean()) / (s.std(ddof=1) + 1e-12))

    # clean rounded distance labels
    if isinstance(sub_df.index, pd.IntervalIndex):
        xticklabels = [f"{int(round(b.left))}-{int(round(b.right))}" for b in sub_df.index]
    else:
        xticklabels = [str(x) for x in sub_df.index]

    pigs = set(pig_genes or [])
    row_labels = [f"{g} *" if g in pigs else g for g in sub_df.columns]

    plt.figure(figsize=figsize)
    xticklabels = _format_bin_labels(sub_df.index)

    ax = sns.heatmap(
        sub_df.T,
        cmap="vlag" if zscore else "magma",
        center=0 if zscore else None,
        cbar_kws={"label": "Z-scored mean expression" if zscore else "Mean expression (log1+)"},
        linewidths=0.2,
        linecolor="white",
    )

    ax.set_xlabel("Distance to plaque (µm, binned)")
    ax.set_ylabel("Gene")
    ax.set_xticklabels(xticklabels, rotation=35, ha="right")
    ax.set_yticklabels(row_labels, rotation=0)

    for tk in ax.get_yticklabels():
        if tk.get_text().endswith("*"):
            tk.set_fontweight("bold")

    if title is None:
        title = f"Top {top_n} genes varying with plaque distance (* = PIG)"
    ax.set_title(title, pad=8)

    plt.tight_layout()
    plt.show()


def _format_bin_labels(index):
    """Return pretty 'low–high' labels for an index whose values may be Intervals."""
    vals = list(index)

    def _fmt(v):
        # Works for Interval and prints fallback for plain values
        if hasattr(v, "left") and hasattr(v, "right"):
            return f"{int(round(v.left))}-{int(round(v.right))}"
        return str(v)

    return [_fmt(v) for v in vals]


def plot_resid_distance_bivariate_pro(
    cells_df: pd.DataFrame,
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
    dist_col: str = "nearest_plaque_center_dist",
    n_resid_bins: int = 3,
    n_dist_bins: int = 3,
    figsize: tuple[float, float] = (6.5, 6),
) -> None:
    """
    Bivariate spatial plot:
      - residual bins: negative / near-zero / positive (quantile-based)
      - distance bins: near / mid / far (quantile-based)
    Colored with a 3×3 bivariate palette and a small legend grid.
    """
    # --- 1) Build working frame -------------------------------------------------
    df = cells_df.copy()
    df["y_true"] = np.asarray(y_true)
    df["y_pred"] = np.asarray(y_pred)
    df["resid"] = df["y_true"] - df["y_pred"]

    # Drop rows with missing essentials
    df = df.dropna(subset=[x_col, y_col, dist_col, "resid"])

    # --- 2) Bin residuals and distances (1D only!) ------------------------------
    # residual bins: quantiles
    resid_q = df["resid"].quantile([0, 1 / 3, 2 / 3, 1]).to_numpy()
    # ensure strictly increasing
    resid_edges = np.unique(resid_q)
    if resid_edges.size < 4:
        # fallback: simple linspace if too many ties
        resid_edges = np.linspace(df["resid"].min(), df["resid"].max(), 4)

    df["resid_bin"] = pd.cut(
        df["resid"],
        bins=resid_edges,
        labels=["neg", "mid", "pos"][: (len(resid_edges) - 1)],
        include_lowest=True,
    )

    # distance bins: quantiles
    dist_q = df[dist_col].quantile([0, 1 / 3, 2 / 3, 1]).to_numpy()
    dist_edges = np.unique(dist_q)
    if dist_edges.size < 4:
        dist_edges = np.linspace(df[dist_col].min(), df[dist_col].max(), 4)

    df["dist_bin"] = pd.cut(
        df[dist_col],
        bins=dist_edges,
        labels=["near", "mid", "far"][: (len(dist_edges) - 1)],
        include_lowest=True,
    )

    # Drop any rows that didn’t get binned (edge cases)
    df = df.dropna(subset=["resid_bin", "dist_bin"])

    # --- 3) Build a simple 3×3 bivariate color grid -----------------------------
    # rows = residual bins (neg → pos), cols = distance bins (near → far)
    # we’ll manually define a 3×3 palette: cooler for near, warmer for far
    bivariate_grid = np.array(
        [
            ["#3b4cc0", "#5470d8", "#6f92f3"],  # negative residual
            ["#79b5b6", "#8fcf9b", "#a6e675"],  # mid residual
            ["#e78b5d", "#f3a151", "#f9c74f"],  # positive residual
        ]
    )

    resid_levels = ["neg", "mid", "pos"]
    dist_levels = ["near", "mid", "far"]

    color_map = {}
    for i, rb in enumerate(resid_levels):
        for j, db in enumerate(dist_levels):
            color_map[(rb, db)] = bivariate_grid[i, j]

    # Map each row to a color
    colors = df.apply(lambda r: color_map.get((r["resid_bin"], r["dist_bin"]), "#cccccc"), axis=1)

    # --- 4) Plot ----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(
        df[x_col],
        df[y_col],
        c=colors,
        s=4,
        alpha=0.9,
        edgecolors="none",
    )
    ax.invert_yaxis()
    ax.set_aspect("equal", "box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Bivariate map: residual × distance to plaque", fontsize=12)

    # --- 5) Mini 3×3 legend grid ------------------------------------------------
    legend_ax = fig.add_axes([0.78, 0.58, 0.16, 0.3])  # [left, bottom, width, height]
    legend_ax.set_xticks(range(3))
    legend_ax.set_yticks(range(3))
    legend_ax.set_xticklabels(dist_levels, rotation=30, ha="right", fontsize=8)
    legend_ax.set_yticklabels(["neg", "mid", "pos"], fontsize=8)
    for i in range(3):
        for j in range(3):
            legend_ax.add_patch(
                plt.Rectangle(
                    (j - 0.5, i - 0.5),
                    width=1,
                    height=1,
                    facecolor=bivariate_grid[i, j],
                    edgecolor="0.7",
                )
            )
    legend_ax.set_xlim(-0.5, 2.5)
    legend_ax.set_ylim(-0.5, 2.5)
    legend_ax.set_title("dist\nnear → far", fontsize=8)
    legend_ax.set_ylabel("resid\nneg → pos", fontsize=8)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    legend_ax.axis("off")

    plt.tight_layout()
    plt.show()


def plot_spatial_bubble_resid_distance_pro(
    cells_df: pd.DataFrame,
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
    dist_col: str = "nearest_plaque_center_dist",
    figsize: tuple[float, float] = (6.5, 6),
) -> None:
    """
    Bubble plot:
      - x,y = spatial coordinates
      - marker size encodes true distance to plaque
      - marker color encodes residual (true - pred), blue→red.
    """
    df = cells_df.copy()
    df["y_true"] = np.asarray(y_true)
    df["y_pred"] = np.asarray(y_pred)
    df["resid"] = df["y_true"] - df["y_pred"]

    df = df.dropna(subset=[x_col, y_col, dist_col, "resid"])

    # Scale marker size from distance (sqrt for nicer dynamic range)
    dist = df[dist_col].to_numpy()
    dist_norm = (dist - dist.min()) / (dist.max() - dist.min() + 1e-9)
    sizes = 10 + 40 * np.sqrt(dist_norm)  # between ~10 and 50

    # Color by residual
    resid = df["resid"].to_numpy()
    vmax = np.percentile(np.abs(resid), 99)
    norm = Normalize(vmin=-vmax, vmax=vmax)

    fig, ax = plt.subplots(figsize=figsize)
    sc = ax.scatter(
        df[x_col],
        df[y_col],
        c=resid,
        s=sizes,
        cmap="coolwarm",
        norm=norm,
        linewidth=0,
        alpha=0.8,
    )
    ax.invert_yaxis()
    ax.set_aspect("equal", "box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Bubble map of distance (size) and residuals (color)", fontsize=12)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Residual (true − predicted distance, µm)")

    # Legend for size (distance)
    for frac, label in zip([0.1, 0.5, 0.9], ["near", "mid", "far"], strict=False):
        d_val = np.quantile(dist, frac)
        size_val = 10 + 40 * np.sqrt((d_val - dist.min()) / (dist.max() - dist.min() + 1e-9))
        ax.scatter([], [], s=size_val, color="gray", alpha=0.7, label=f"{label} distance")
    ax.legend(
        loc="upper right",
        title="Distance encoding",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
    )

    plt.tight_layout()
    plt.show()


def plot_kde_true_pred_distance_pro(
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    model_name: str = "xgb",
    figsize: tuple[float, float] = (6, 4),
) -> None:
    """
    KDE plot of true vs predicted plaque distance on same axes.
    """
    y_true_arr = np.asarray(y_true, dtype=float)
    y_pred_arr = np.asarray(y_pred, dtype=float)

    mask = np.isfinite(y_true_arr) & np.isfinite(y_pred_arr)
    y_true_arr = y_true_arr[mask]
    y_pred_arr = y_pred_arr[mask]

    fig, ax = plt.subplots(figsize=figsize)
    sns.kdeplot(y_true_arr, ax=ax, label="True distance", linewidth=2)
    sns.kdeplot(
        y_pred_arr, ax=ax, label=f"Predicted distance ({model_name})", linewidth=2, linestyle="--"
    )

    ax.set_xlabel("Distance to plaque (µm)")
    ax.set_ylabel("Density")
    ax.set_title(f"KDE of true vs predicted plaque distance ({model_name})", fontsize=12)
    ax.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def plot_resid_vs_distance_scatter_pro(
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    model_name: str = "xgb",
    figsize: tuple[float, float] = (6, 4),
) -> None:
    """
    Scatter plot: residual (true - predicted) vs true distance.
    Pearson r is reported in the legend.
    """
    y_true_arr = np.asarray(y_true, dtype=float)
    y_pred_arr = np.asarray(y_pred, dtype=float)
    resid = y_true_arr - y_pred_arr

    mask = np.isfinite(y_true_arr) & np.isfinite(resid)
    x = y_true_arr[mask]
    r = resid[mask]

    if x.size < 2:
        print("Not enough points for correlation.")
        return

    pearson_r, pval = pearsonr(x, r)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x, r, s=6, alpha=0.4, label="Cells")

    ax.axhline(0, color="0.5", linestyle="--", linewidth=1)
    ax.set_xlabel("True distance to plaque (µm)")
    ax.set_ylabel("Residual (true − predicted, µm)")
    ax.set_title(f"Residuals vs true distance ({model_name})", fontsize=12)

    # Fake handle for legend entry with r
    corr_label = f"Pearson r = {pearson_r:.2f} (p = {pval:.1e})"
    ax.scatter([], [], c="none", edgecolor="none", label=corr_label)

    ax.legend(frameon=True, fontsize=9)
    plt.tight_layout()
    plt.show()


def _format_interval(iv) -> str:
    """Format an interval as rounded bounds."""
    left = float(iv.left)
    right = float(iv.right)
    return f"{left:.1f} to {right:.1f}"


def _relative_luminance(rgb: tuple[float, float, float]) -> float:
    """WCAG-ish relative luminance for deciding white vs black text on a color."""
    r, g, b = rgb

    def _lin(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    r_lin, g_lin, b_lin = _lin(r), _lin(g), _lin(b)
    return 0.2126 * r_lin + 0.7152 * g_lin + 0.0722 * b_lin


def _bivariate_palette_rowwise(
    *,
    n_resid: int,
    n_dist: int,
    neg_cmap: str = "Blues",
    mid_cmap: str = "Greys",
    pos_cmap: str = "Oranges",
    dist_lo: float = 0.25,
    dist_hi: float = 0.95,
) -> "np.ndarray":
    """
    Build a bivariate palette where:
      - residual bin selects the row hue family (neg=blue, mid=grey, pos=orange)
      - distance bin selects lightness/saturation within that family (near=light, far=dark)

    This intentionally creates strong, monotonic lightness steps so bins remain
    distinguishable at very small marker sizes. :contentReference[oaicite:3]{index=3}
    """
    import numpy as np
    import matplotlib as mpl

    if n_resid != 3:
        raise ValueError("rowwise palette is designed for n_bins=3 (neg/mid/pos residual bins).")

    def sample_rgb(cmap_name: str) -> list[tuple[float, float, float]]:
        cmap = mpl.cm.get_cmap(cmap_name)
        vals = np.linspace(dist_lo, dist_hi, n_dist)
        return [tuple(cmap(v)[:3]) for v in vals]

    neg_row = sample_rgb(neg_cmap)
    mid_row = sample_rgb(mid_cmap)
    pos_row = sample_rgb(pos_cmap)

    # rows: residual bins low->mid->high, cols: distance bins near->far
    return np.array([neg_row, mid_row, pos_row], dtype=object)


def plot_bivariate_resid_distance_spatial(
    cells_df: "pd.DataFrame",
    y_true: "pd.Series",
    y_pred: "np.ndarray",
    *,
    model_name: str = "model",
    dist_col: str = "nearest_plaque_center_dist",
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
    n_bins: int = 3,
    point_size: float = 3.0,
    point_alpha: float = 0.90,
    legend_show_counts: bool = True,
    legend_count_fmt: str = "{:,}",
    legend_count_fontsize: int = 7,
    legend_count_min: int | None = None,
    palette_dist_lo: float = 0.25,
    palette_dist_hi: float = 0.95,
    ax: "plt.Axes | None" = None,
    # ---  plaque overlay ---
    plaque_x_col: str = "plaque_x",
    plaque_y_col: str = "plaque_y",
    show_plaques: bool = True,
    plaque_marker: str = "*",
    plaque_color: str = "red",
    plaque_edgecolor: str = "red",
    plaque_size: float = 5,
    plaque_alpha: float = 0.95,
    plaque_linewidth: float = 0.8,
    plaque_zorder: float = 5,
    plaque_round_decimals: int | None = None,
    plaque_label: str = "Plaque",
    show_plaque_legend: bool = False,
    show_marker_legend: bool = True,
    marker_legend_loc: str = "lower right",
    marker_legend_fontsize: int = 8,
    marker_legend_framealpha: float = 0.90,
    marker_legend_cell_markersize: float = 4.0,
    marker_legend_plaque_markersize: float = 8.0,
    marker_legend_pad: float = 0.25,
) -> "plt.Axes":
    """
    Bivariate choropleth-style plot of residual bin × distance bin in spatial coordinates.

    Adds:
      - optional overlay of unique plaque locations from (plaque_x_col, plaque_y_col)
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    if n_bins != 3:
        raise ValueError("This implementation expects n_bins=3 to map residual bins to neg/mid/pos.")

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    df = cells_df.copy()

    # Robust alignment: assume y_true corresponds row-wise to cells_df
    resid = np.asarray(y_true, dtype=float) - np.asarray(y_pred, dtype=float)
    df["_resid_"] = resid
    df["_dist_"] = pd.to_numeric(df[dist_col], errors="coerce")
    df["_x_"] = pd.to_numeric(df[x_col], errors="coerce")
    df["_y_"] = pd.to_numeric(df[y_col], errors="coerce")

    # Drop non-finite (cells)
    mask = np.isfinite(df["_resid_"].to_numpy()) & np.isfinite(df["_dist_"].to_numpy())
    mask &= np.isfinite(df["_x_"].to_numpy()) & np.isfinite(df["_y_"].to_numpy())
    df = df.loc[mask].copy()

    # Bin residuals and distances
    df["_resid_bin_"] = pd.qcut(df["_resid_"], q=n_bins, duplicates="drop")
    df["_dist_bin_"] = pd.qcut(df["_dist_"], q=n_bins, duplicates="drop")

    # Enforce 3 bins (if qcut drops bins due to ties, fail loudly)
    if df["_resid_bin_"].cat.categories.size != 3 or df["_dist_bin_"].cat.categories.size != 3:
        raise ValueError(
            "qcut produced fewer than 3 bins (ties). Consider adding small jitter "
            "or switching to pd.cut with fixed edges."
        )

    resid_cats = df["_resid_bin_"].cat.categories  # low->high
    dist_cats = df["_dist_bin_"].cat.categories    # near->far (low->high)

    df["_resid_idx_"] = df["_resid_bin_"].cat.codes  # 0,1,2
    df["_dist_idx_"] = df["_dist_bin_"].cat.codes    # 0,1,2

    # Palette
    colors = _bivariate_palette_rowwise(
        n_resid=3,
        n_dist=3,
        dist_lo=palette_dist_lo,
        dist_hi=palette_dist_hi,
    )

    # Point colors
    ridx = df["_resid_idx_"].to_numpy()
    didx = df["_dist_idx_"].to_numpy()
    point_colors = [colors[r, d] for r, d in zip(ridx, didx, strict=False)]

    # Cells
    ax.scatter(
        df["_x_"],
        df["_y_"],
        c=point_colors,
        s=point_size,
        alpha=point_alpha,
        edgecolors="none",
        rasterized=True,
        zorder=2,
    )

    # --- Plaques (unique coordinates) ---
    if show_plaques and (plaque_x_col in cells_df.columns) and (plaque_y_col in cells_df.columns):
        plaques = cells_df[[plaque_x_col, plaque_y_col]].copy()
        plaques["_px_"] = pd.to_numeric(plaques[plaque_x_col], errors="coerce")
        plaques["_py_"] = pd.to_numeric(plaques[plaque_y_col], errors="coerce")
        plaques = plaques.loc[np.isfinite(plaques["_px_"].to_numpy()) & np.isfinite(plaques["_py_"].to_numpy())]

        # Optional rounding helps if plaque coords are floats with tiny noise
        if plaque_round_decimals is not None:
            plaques["_px_"] = plaques["_px_"].round(plaque_round_decimals)
            plaques["_py_"] = plaques["_py_"].round(plaque_round_decimals)

        plaques = plaques.drop_duplicates(subset=["_px_", "_py_"])

        ax.scatter(
            plaques["_px_"],
            plaques["_py_"],
            marker=plaque_marker,
            s=plaque_size,
            c=plaque_color,
            alpha=plaque_alpha,
            edgecolors=plaque_edgecolor,
            linewidths=plaque_linewidth,
            rasterized=False,   # keep plaque markers crisp in vector exports
            zorder=plaque_zorder,
            label=plaque_label,
        )

    ax.set_aspect("equal", "box")
    ax.set_xlabel("X (µm)", fontsize=10)
    ax.set_ylabel("Y (µm)", fontsize=10)

    import matplotlib.ticker as mticker
    ax.xaxis.set_major_locator(mticker.MultipleLocator(1000))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(1000))
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax.tick_params(axis="both", which="major", direction="out", length=4, width=1, labelsize=9)
    ax.tick_params(axis="both", which="minor", direction="out", length=2, width=0.8)
    ax.invert_yaxis()
    ax.set_title("Residual vs distance bivariate map", fontsize=10)

    # Counts for legend squares
    counts = (
        df.groupby(["_resid_idx_", "_dist_idx_"], observed=True)
        .size()
        .reindex(pd.MultiIndex.from_product([range(3), range(3)]), fill_value=0)
        .to_numpy()
        .reshape(3, 3)
    )

    # Legend inset (placed outside left)
    inset = ax.inset_axes([-0.44, 0.25, 0.30, 0.50])

    for r in range(3):
        for d in range(3):
            y0 = 2 - r
            x0 = d
            rgb = colors[r, d]
            inset.add_patch(plt.Rectangle((x0, y0), 1, 1, color=rgb, transform=inset.transData))

            if legend_show_counts:
                n = int(counts[r, d])
                if (legend_count_min is None) or (n >= legend_count_min):
                    txt_color = "black" if _relative_luminance(rgb) > 0.55 else "white"
                    inset.text(
                        x0 + 0.5, y0 + 0.5,
                        legend_count_fmt.format(n),
                        ha="center", va="center",
                        fontsize=legend_count_fontsize,
                        color=txt_color,
                    )

    inset.set_xlim(0, 3)
    inset.set_ylim(0, 3)

    inset.set_xticks(np.arange(3) + 0.5)
    inset.set_yticks(np.arange(3) + 0.5)

    inset.set_xticklabels(
        [f"{dist_cats[j].left:.0f}-{dist_cats[j].right:.0f} µm" for j in range(3)],
        rotation=45,
        ha="right",
        fontsize=7,
    )
    inset.set_yticklabels([_format_interval(c) for c in resid_cats[::-1]], fontsize=7)

    inset.set_xlabel("Distance bin", fontsize=7)
    inset.set_ylabel("Residual bin", fontsize=7)
    inset.tick_params(length=0)
    for spine in inset.spines.values():
        spine.set_visible(False)
    
    # --- Marker legend (bottom-right): star = plaque centroid; circles = cells (all 9 colors) ---
    if show_marker_legend:
        from matplotlib.lines import Line2D
        from matplotlib.legend_handler import HandlerTuple

        handles = []
        labels = []

        # Plaque centroid handle (only if plaques are actually shown)
        if show_plaques:
            plaque_handle = Line2D(
                [0], [0],
                marker=plaque_marker,
                linestyle="None",
                color="none",
                markerfacecolor=plaque_color,
                markeredgecolor=plaque_edgecolor,
                markeredgewidth=plaque_linewidth,
                markersize=marker_legend_plaque_markersize,
            )
            handles.append(plaque_handle)
            labels.append("Plaque centroid")

        # Cells: show all 9 colors, but compactly as 3 rows (residual bins) × 3 cols (distance bins)
        # Each legend row is a tuple of 3 colored circles (near->far left->right)
        cell_row_labels = ["Cell (low residual)", "Cell (mid residual)", "Cell (high residual)"]
        for r in range(3):  # 0,1,2 matches your palette rows
            row = tuple(
                Line2D(
                    [0], [0],
                    marker="o",
                    linestyle="None",
                    color="none",
                    markerfacecolor=colors[r, d],
                    markeredgecolor="none",
                    markersize=marker_legend_cell_markersize,
                )
                for d in range(3)
            )
            handles.append(row)
            labels.append(cell_row_labels[r])

        ax.legend(
            handles,
            labels,
            handler_map={tuple: HandlerTuple(ndivide=None, pad=marker_legend_pad)},
            loc=marker_legend_loc,
            frameon=True,
            framealpha=marker_legend_framealpha,
            fontsize=marker_legend_fontsize,
            borderpad=0.35,
            handletextpad=0.6,
            labelspacing=0.5,
        )
    return ax


def plot_radius_color_spatial(
    cells_df: "pd.DataFrame",
    y_true: "pd.Series",
    y_pred: "np.ndarray",
    *,
    model_name: str = "model",
    dist_col: str = "nearest_plaque_center_dist",
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
    # NEW defaults for dense data:
    render: str = "binned_rgba",  # {"binned_rgba", "scatter"}
    gridsize: int = 300,
    dist_cmap: str = "cividis",   # sequential; distance as color :contentReference[oaicite:6]{index=6}
    min_alpha: float = 0.05,
    max_alpha: float = 1.00,
    # scatter-mode controls (only used if render="scatter")
    max_points: int = 250_000,
    size_min: float = 4.0,
    size_max: float = 36.0,
    size_gamma: float = 0.6,
    scatter_alpha: float = 0.35,
    ax: "plt.Axes | None" = None,
) -> "plt.Axes":
    """
    Residuals in space (improved for heavy overplotting).

    Default render="binned_rgba":
      - distance to plaque -> color (sequential colormap)
      - mean |residual|    -> opacity (alpha)
      This avoids the “all points overlap so size doesn't register” problem and is a
      standard approach for large scatter/dense point clouds. :contentReference[oaicite:7]{index=7}

    Optional render="scatter":
      - distance -> color (sequential)
      - |residual| -> size
      - sign(residual) -> marker shape (o vs v)
      Automatically subsamples if too many points.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.patches import Patch

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    df = cells_df.copy()
    resid = np.asarray(y_true, dtype=float) - np.asarray(y_pred, dtype=float)

    df["_resid_"] = resid
    df["_abs_resid_"] = np.abs(resid)
    df["_dist_"] = pd.to_numeric(df[dist_col], errors="coerce")
    df["_x_"] = pd.to_numeric(df[x_col], errors="coerce")
    df["_y_"] = pd.to_numeric(df[y_col], errors="coerce")

    # Drop non-finite
    m = np.isfinite(df["_x_"].to_numpy()) & np.isfinite(df["_y_"].to_numpy())
    m &= np.isfinite(df["_dist_"].to_numpy()) & np.isfinite(df["_resid_"].to_numpy())
    df = df.loc[m].copy()

    x = df["_x_"].to_numpy()
    y = df["_y_"].to_numpy()
    dist = df["_dist_"].to_numpy()
    absr = df["_abs_resid_"].to_numpy()
    r = df["_resid_"].to_numpy()

    ax.set_aspect("equal", "box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("x coordinate (µm)", fontsize=10)
    ax.set_ylabel("y coordinate (µm)", fontsize=10)

    if render == "binned_rgba":
        # 2D binning in (y, x) so imshow aligns naturally
        H, yedges, xedges = np.histogram2d(y, x, bins=gridsize)
        sum_dist, _, _ = np.histogram2d(y, x, bins=[yedges, xedges], weights=dist)
        sum_absr, _, _ = np.histogram2d(y, x, bins=[yedges, xedges], weights=absr)

        with np.errstate(invalid="ignore", divide="ignore"):
            mean_dist = np.divide(sum_dist, H, out=np.full_like(sum_dist, np.nan), where=H > 0)
            mean_absr = np.divide(sum_absr, H, out=np.full_like(sum_absr, np.nan), where=H > 0)

        # Robust scaling
        if np.isfinite(mean_dist).any():
            d_lo, d_hi = np.nanpercentile(mean_dist, [1, 99])
        else:
            d_lo, d_hi = float(np.nanmin(dist)), float(np.nanmax(dist))

        if np.isfinite(mean_absr).any():
            a_lo, a_hi = np.nanpercentile(mean_absr, [50, 99])
        else:
            a_lo, a_hi = float(np.nanmin(absr)), float(np.nanmax(absr))

        if not np.isfinite(d_lo) or not np.isfinite(d_hi) or d_lo == d_hi:
            d_lo, d_hi = float(np.nanmin(dist)), float(np.nanmax(dist))
        if not np.isfinite(a_lo) or not np.isfinite(a_hi) or a_lo == a_hi:
            a_lo, a_hi = float(np.nanmin(absr)), float(np.nanmax(absr) + 1e-9)

        dist_norm = mpl.colors.Normalize(vmin=d_lo, vmax=d_hi, clip=True)
        abs_norm = mpl.colors.Normalize(vmin=a_lo, vmax=a_hi, clip=True)

        cmap = mpl.cm.get_cmap(dist_cmap)

        # Convert mean_dist -> RGBA via colormap
        mean_dist_filled = np.nan_to_num(mean_dist, nan=d_lo)
        rgba = cmap(dist_norm(mean_dist_filled))  # (H,W,4)

        # Opacity from mean |residual|
        mean_absr_filled = np.nan_to_num(mean_absr, nan=a_lo)
        alpha = min_alpha + (max_alpha - min_alpha) * abs_norm(mean_absr_filled)
        alpha = np.clip(alpha, 0.0, 1.0)
        alpha[H == 0] = 0.0
        rgba[..., 3] = alpha

        ax.imshow(
            rgba,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            origin="lower",
            interpolation="nearest",
            rasterized=True,
        )
        ax.invert_yaxis()

        ax.set_title(
            "Residuals in space\ncolor = mean distance, opacity = mean |residual|",
            fontsize=10,
        )

        # Distance colorbar
        sm = mpl.cm.ScalarMappable(norm=dist_norm, cmap=cmap)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("Mean distance to plaque center (µm)")

        # Opacity legend (categorical, but derived from continuous mapping)
        handles = [
            Patch(facecolor="black", alpha=min_alpha + 0.25 * (max_alpha - min_alpha), label="small |residual|"),
            Patch(facecolor="black", alpha=min_alpha + 0.60 * (max_alpha - min_alpha), label="medium |residual|"),
            Patch(facecolor="black", alpha=min_alpha + 0.95 * (max_alpha - min_alpha), label="large |residual|"),
        ]
        ax.legend(
            handles=handles,
            title="Mean |residual|\n(opacity)",
            loc="upper right",
            frameon=True,
            framealpha=0.85,
        )

        return ax

    if render == "scatter":
        # Subsample to keep plotting interactive and reduce overplotting
        n = x.shape[0]
        if n > max_points:
            rng = np.random.default_rng(0)
            idx = rng.choice(n, size=max_points, replace=False)
            x, y, dist, absr, r = x[idx], y[idx], dist[idx], absr[idx], r[idx]

        # Robust scaling for size
        vmax = np.nanpercentile(absr, 99) if np.isfinite(absr).any() else float(np.nanmax(absr))
        vmax = max(vmax, 1e-9)
        t = np.clip(absr / vmax, 0.0, 1.0) ** size_gamma
        sizes = size_min + (size_max - size_min) * t

        # Robust scaling for distance color
        d_lo, d_hi = np.nanpercentile(dist, [1, 99]) if np.isfinite(dist).any() else (np.nanmin(dist), np.nanmax(dist))
        dist_norm = mpl.colors.Normalize(vmin=d_lo, vmax=d_hi, clip=True)

        # Split by sign (marker encodes sign)
        pos = r >= 0
        neg = ~pos

        sc1 = ax.scatter(
            x[pos],
            y[pos],
            c=dist[pos],
            s=sizes[pos],
            cmap=dist_cmap,
            norm=dist_norm,
            alpha=scatter_alpha,
            marker="o",
            edgecolors="none",
            rasterized=True,
            label="residual ≥ 0",
        )
        ax.scatter(
            x[neg],
            y[neg],
            c=dist[neg],
            s=sizes[neg],
            cmap=dist_cmap,
            norm=dist_norm,
            alpha=scatter_alpha,
            marker="v",
            edgecolors="none",
            rasterized=True,
            label="residual < 0",
        )

        ax.invert_yaxis()
        ax.set_title(
            "Residuals in space\ncolor = distance, size = |residual|, marker = sign",
            fontsize=10,
        )

        cbar = plt.colorbar(sc1, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("Distance to plaque center (µm)")

        ax.legend(loc="upper right", frameon=True, framealpha=0.85)
        return ax

    raise ValueError(f"Unknown render='{render}'. Use 'binned_rgba' or 'scatter'.")


def plot_true_pred_kde(
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    model_name: str = "Model",
    figsize: tuple[float, float] = (10, 4),
    show_qq: bool = True,
    kde_points: int = 400,
    qq_points: int | None = 200,   # None => use min(len(y_true), len(y_pred))
    qq_ref_line: bool = True,
    qq_log: bool = False,          # helpful for heavy right tails (requires all values > 0)
) -> None:
    """
    KDE curves of true vs predicted distance, optionally with a two-sample Q-Q plot.

    Left panel: KDE(y_true) vs KDE(y_pred)
    Right panel (optional): Q-Q plot comparing empirical quantiles of y_true vs y_pred
    """
    # Clean + coerce numeric + drop NaN/inf
    y_true_s = pd.Series(y_true)
    y_pred_s = pd.Series(y_pred)

    y_true_arr = pd.to_numeric(y_true_s, errors="coerce").to_numpy()
    y_pred_arr = pd.to_numeric(y_pred_s, errors="coerce").to_numpy()

    y_true_arr = y_true_arr[np.isfinite(y_true_arr)]
    y_pred_arr = y_pred_arr[np.isfinite(y_pred_arr)]

    if y_true_arr.size == 0 or y_pred_arr.size == 0:
        raise ValueError("y_true and y_pred must each contain at least one finite numeric value.")

    # Figure layout
    if show_qq:
        fig, (ax_kde, ax_qq) = plt.subplots(ncols=2, figsize=figsize)
    else:
        fig, ax_kde = plt.subplots(figsize=figsize)
        ax_qq = None

    # -------- KDE panel --------
    common_min = float(min(y_true_arr.min(), y_pred_arr.min()))
    common_max = float(max(y_true_arr.max(), y_pred_arr.max()))
    xs = np.linspace(common_min, common_max, kde_points)

    kde_true = gaussian_kde(y_true_arr)
    ax_kde.plot(xs, kde_true(xs), label="True distance", lw=2)

    kde_pred = gaussian_kde(y_pred_arr)
    ax_kde.plot(xs, kde_pred(xs), label="Predicted distance", lw=2, linestyle="--")

    ax_kde.set_xlabel("Distance to plaque (µm)")
    ax_kde.set_ylabel("Density")
    ax_kde.set_title(f"{model_name}: KDE")
    ax_kde.legend(frameon=False)
    ax_kde.grid(alpha=0.25, linestyle="--", linewidth=0.5)

    # -------- Two-sample Q-Q panel --------
    if show_qq and ax_qq is not None:
        n = min(y_true_arr.size, y_pred_arr.size)
        if qq_points is not None:
            n = min(n, int(qq_points))

        if n < 2:
            raise ValueError("Need at least 2 points in each array to draw a Q-Q plot.")

        # Plotting positions (avoid exact 0 and 1)
        ps = (np.arange(1, n + 1) - 0.5) / n

        # Empirical quantiles (works for non-normal, skewed, heavy-tailed data)
        q_true = np.quantile(y_true_arr, ps, method="linear")
        q_pred = np.quantile(y_pred_arr, ps, method="linear")

        ax_qq.scatter(q_true, q_pred, s=18, alpha=0.7)

        if qq_ref_line:
            lo = float(min(q_true.min(), q_pred.min()))
            hi = float(max(q_true.max(), q_pred.max()))
            ax_qq.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)

        ax_qq.set_xlabel("True quantiles")
        ax_qq.set_ylabel("Predicted quantiles")
        ax_qq.set_title(f"{model_name}: Q-Q (two-sample)")
        ax_qq.grid(alpha=0.25, linestyle="--", linewidth=0.5)

        if qq_log:
            if (q_true <= 0).any() or (q_pred <= 0).any():
                raise ValueError("qq_log=True requires all values to be strictly positive.")
            ax_qq.set_xscale("log")
            ax_qq.set_yscale("log")

    plt.tight_layout()
    plt.show()


def plot_residual_vs_distance(
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    model_name: str = "Model",
    figsize: tuple[float, float] = (6, 4),
) -> None:
    """
    Scatter plot of residual (y) vs true distance (x), annotated
    with Pearson correlation in the legend.
    """
    y_true = pd.Series(y_true).astype(float)
    y_pred = pd.Series(y_pred).astype(float)
    residuals = y_true - y_pred

    mask = ~(y_true.isna() | residuals.isna())
    y_true = y_true[mask]
    residuals = residuals[mask]

    if len(y_true) == 0:
        raise ValueError("No finite values for residual vs distance scatter.")

    r = np.corrcoef(y_true.to_numpy(), residuals.to_numpy())[0, 1]

    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(
        x=y_true,
        y=residuals,
        s=10,
        alpha=0.5,
        edgecolor=None,
        ax=ax,
        # label=f"Cells (n={len(y_true)})",
    )

    # Optional LOWESS-like smoother via rolling median (cheap)
    order = np.argsort(y_true.to_numpy())
    x_sorted = y_true.to_numpy()[order]
    y_sorted = residuals.to_numpy()[order]
    window = max(20, len(x_sorted) // 50)
    if window > 5:
        y_smooth = pd.Series(y_sorted).rolling(window, center=True).median()
        ax.plot(
            x_sorted,
            y_smooth,
            color="black",
            lw=2,
            label="Rolling median residual",
        )

    ax.axhline(0, color="red", linestyle="--", linewidth=1)

    ax.set_xlabel("True distance to plaque, µm")
    ax.set_ylabel("True - predicted, µm")
    # ax.set_title("Residuals vs distance")

    # Add Pearson r as a separate legend entry
    # dummy = plt.Line2D([], [], color="none", label=f"Pearson r = {r:.2f}")
    print(f"Pearson r = {r:.2f}")
    handles, labels = ax.get_legend_handles_labels()
    # handles.append(dummy)
    # labels.append(dummy.get_label())
    ax.legend(handles, labels, frameon=False, loc="upper right")

    ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.show()
