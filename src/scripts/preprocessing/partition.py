from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, kurtosis, skew

# Default metadata columns for your Xenium dataset
DEFAULT_META_COLS = {
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
}


def get_gene_columns(df: pd.DataFrame, meta_cols: Iterable[str] = DEFAULT_META_COLS) -> list[str]:
    """Return numeric columns not in meta as gene feature columns."""
    num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    meta_set = set(meta_cols)
    return [c for c in num_cols if c not in meta_set]


def split_meta_vs_genes(
    df: pd.DataFrame, meta_cols: Iterable[str] = DEFAULT_META_COLS
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Split DataFrame into meta and gene feature DataFrames."""
    meta_set = set(meta_cols)
    meta_df = df.loc[:, [c for c in df.columns if c in meta_set]].copy()
    gene_cols = get_gene_columns(df, meta_cols)
    genes_df = df.loc[:, gene_cols].copy()
    info = {
        "n_numeric": df.select_dtypes(include=[np.number]).shape[1],
        "n_genes": len(gene_cols),
        "n_meta_found": meta_df.shape[1],
    }
    return meta_df, genes_df, info


@dataclass(frozen=True)
class QCBounds:
    lo_genes: float
    hi_genes: float


def _iqr_bounds(s: pd.Series, mult: float = 3.0, floor: float = 1.0) -> tuple[float, float]:
    q1, q3 = s.quantile([0.25, 0.75])
    iqr = q3 - q1
    lo = max(floor, float(q1 - mult * iqr))
    hi = float(q3 + mult * iqr)
    return lo, hi


def add_qc_metrics(df: pd.DataFrame, gene_cols: Iterable[str]) -> pd.DataFrame:
    """
    Return a copy with n_genes (nonzero gene features) and n_counts (sum across genes).
    """
    out = df.copy()
    gc = list(gene_cols)
    # Counts the number of genes with above zero expression level in each cell
    # Dimension: cells x 1
    out["n_genes"] = (out[gc] > 0).sum(axis=1).astype("int64")
    return out


def filter_cells_iqr(
    df: pd.DataFrame,
    gene_cols: Iterable[str],
    area_col: str = "cell_area",
    nucleus_col: str = "nucleus_area",
    iqr_mult: float = 3.0,
    min_floor_counts: float = 1.0,
    min_floor_genes: float = 1.0,
) -> tuple[pd.DataFrame, pd.Series, QCBounds, dict[str, int], pd.DataFrame]:
    # TODO: update docstring to be more granular
    """
    Add QC metrics, compute IQR-based bounds for n_counts / n_genes, build a mask,
    and return filtered cells.

    Returns:
        df_clean: filtered dataframe (copy)
        mask: boolean mask of kept rows (aligned with input df index)
        bounds: QCBounds with thresholds used
        summary_numbers: dict with kept/removed counts
        describe_tbl: transposed describe() of key QC columns in df_clean
    """
    df_qc = add_qc_metrics(df, gene_cols)

    lo_genes, hi_genes = _iqr_bounds(df_qc["n_genes"], mult=iqr_mult, floor=min_floor_genes)

    area_ok = df_qc.get(area_col, pd.Series(1, index=df_qc.index, dtype="int64")) > 0
    nuc_ok = df_qc.get(nucleus_col, pd.Series(0, index=df_qc.index, dtype="int64")) > 0
    genes_ok = df_qc["n_genes"].between(lo_genes, hi_genes)

    mask = area_ok & nuc_ok & genes_ok
    df_clean = df_qc.loc[mask].copy()

    kept = int(mask.sum())
    total = int(len(df_qc))
    removed = total - kept

    bounds = QCBounds(
        lo_genes=lo_genes,
        hi_genes=hi_genes,
    )
    summary_numbers = {"kept": kept, "total": total, "removed": removed}

    cols_to_show = [
        c for c in ["n_counts", "n_genes", area_col, nucleus_col] if c in df_clean.columns
    ]
    describe_tbl = df_clean[cols_to_show].describe().T

    return df_clean, mask, bounds, summary_numbers, describe_tbl


def summarize_genes(gene_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-gene summary statistics:
    mean, std, median, MAD, skew, kurtosis, and nonzero fraction.

    Args:
        gene_df (pd.DataFrame): DataFrame of shape (cells x genes), only gene columns.

    Returns:
        pd.DataFrame: Summary indexed by gene name, sorted alphabetically.
    """
    g = gene_df

    summary = (
        pd.DataFrame(
            {
                "gene": g.columns,
                "mean": g.mean(axis=0).values,
                "std": g.std(axis=0, ddof=1).values,
                "median": g.median(axis=0).values,
                "mad": g.subtract(g.median()).abs().median(axis=0).values,
                "skew": g.apply(lambda s: skew(s, bias=False, nan_policy="omit")).values,
                "kurtosis": g.apply(
                    lambda s: kurtosis(s, fisher=True, bias=False, nan_policy="omit")
                ).values,
                "nonzero_frac": (g > 0).mean(axis=0).values,
            }
        )
        .set_index("gene")
        .sort_index()
    )

    return summary


def compute_pigs_zscores(
    df: pd.DataFrame,
    gene_cols: list[str],
    pigs: list[str] | None = None,
    clip: float = 10.0,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """
    Compute z-scored expression matrix for a predefined set of PIG genes.

    Args:
        df (pd.DataFrame): Clean expression data (cells x all genes).
        gene_cols (List[str]): List of gene feature column names in df.
        pigs (List[str], optional): List of target PIG genes. If None, uses default.
        clip (float): Clip z-scores to +/- this value.

    Returns:
        pigs_z (pd.DataFrame): Z-scored expression matrix for present PIG genes.
        info (dict): {'present': [...], 'missing': [...]}
    """
    if pigs is None:
        pigs = [
            "Hexb",
            "Cst3",
            "CD63",
            "C4b",
            "Ctsd",
            "B2m",
            "H2-K1",
            "Apoe",
            "Gfap",
            "Nrep",
            "Serpina3n",
            "Cd74",
            "Cxcl10",
            "Vim",
            "S100a6",
            "Ifit3",
        ]

    # Map lowercased gene names to actual column names
    lower_map = {c.lower(): c for c in gene_cols}
    pigs_present = [lower_map[p.lower()] for p in pigs if p.lower() in lower_map]
    pigs_missing = [p for p in pigs if p.lower() not in lower_map]

    if not pigs_present:
        return pd.DataFrame(index=df.index), {"present": [], "missing": pigs_missing}

    # Extract expression matrix for present genes
    Xp = df[pigs_present].to_numpy(dtype=float)
    mu = np.nanmean(Xp, axis=0, keepdims=True)
    sd = np.nanstd(Xp, axis=0, keepdims=True)
    sd[sd == 0] = 1.0  # avoid division by zero
    Xp_z = (Xp - mu) / sd

    if clip is not None:
        Xp_z = np.clip(Xp_z, -clip, clip)

    pigs_z = pd.DataFrame(Xp_z, index=df.index, columns=pigs_present)

    return pigs_z, {"present": pigs_present, "missing": pigs_missing}


def plot_gene_distributions(
    df: pd.DataFrame,
    genes: Sequence[str],
    bins: int = 50,
    cols: int = 4,
    title: str = "Per-gene distributions (histogram + KDE)",
    figsize: tuple[int, int] | None = None,
    show: bool = True,
    save_path: str | None = None,
):
    """
    Plot histogram + KDE for a set of genes.

    Args:
        df (pd.DataFrame): Expression DataFrame (cells x genes).
        genes (Sequence[str]): List of gene column names to plot.
        bins (int): Number of histogram bins.
        cols (int): Number of subplot columns.
        title (str): Figure title.
        figsize (tuple, optional): (width, height) in inches.
            Defaults to (4*cols, 3.2*rows).
        show (bool): Whether to display the plot.
        save_path (str, optional): If provided, save the figure to this path.

    Returns:
        matplotlib.figure.Figure: The created figure.
    """
    n = len(genes)
    rows = (n + cols - 1) // cols
    if figsize is None:
        figsize = (4 * cols, 3.2 * rows)

    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    axes = np.array(axes).ravel()

    for ax, gene in zip(axes, genes, strict=False):
        x = df[gene].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        ax.hist(x, bins=bins, density=True, alpha=0.45)
        try:
            kde = gaussian_kde(x[x > 0] if np.count_nonzero(x) > 5 else x)
            xs = np.linspace(x.min(), x.max(), 400)
            ax.plot(xs, kde(xs), linewidth=1.75)
        except Exception:
            pass
        ax.set_title(gene)
        ax.set_xlabel("expression")
        ax.set_ylabel("density")

    # Hide unused axes
    for k in range(n, len(axes)):
        axes[k].axis("off")

    fig.suptitle(title, y=1.02, fontsize=14)
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def compute_weird_gene_scores(
    gene_df: pd.DataFrame,
    clip_z: float = 3.0,
    ddof_std: int = 0,
    eps: float = 1e-9,
    PIGs: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute per-gene metrics and a composite 'weird_score' capturing
    zero-inflation, dispersion, shape (skew/kurtosis), and tail ratio.

    Args:
        gene_df: DataFrame of shape (cells x genes) with numeric gene columns.
        clip_z: Clip value for the Z-transformed components (±clip_z).
        ddof_std: ddof for standard deviation (0 matches your original).
        eps: small constant to avoid division by zero.

    Returns:
        S: DataFrame with metrics and 'weird_score', sorted desc by weird_score.
        Z: DataFrame with standardized components used to build weird_score.
    """
    G = gene_df

    # Basic counts / fractions
    nz = (G > 0).sum(axis=0).astype(float)
    n = float(len(G))
    zero_frac = 1.0 - nz / max(n, eps)

    # Moments / dispersion
    mean_ = G.mean(axis=0)
    std_ = G.std(axis=0, ddof=ddof_std)
    cv = std_ / (mean_.replace(0, np.nan))

    # Shape: skew/kurt (omit NaNs)
    sk = G.apply(lambda s: skew(s, bias=False, nan_policy="omit"))
    ku = G.apply(lambda s: kurtosis(s, fisher=True, bias=False, nan_policy="omit"))

    # Tails
    q10 = G.quantile(0.10, axis=0)
    q90 = G.quantile(0.90, axis=0)
    tail_ratio = (q90 + eps) / (q10 + eps)

    # Aggregate table
    S = pd.DataFrame(
        {
            "zero_frac": zero_frac,
            "nonzero_frac": 1.0 - zero_frac,
            "mean": mean_,
            "std": std_,
            "cv": cv,
            "skew": sk,
            "kurtosis": ku,
            "q10": q10,
            "q90": q90,
            "tail_ratio": tail_ratio,
            "isPIG": G.columns.isin(PIGs),
        }
    ).replace([np.inf, -np.inf], np.nan)

    # Z-standardize selected components
    def _zcol(col: pd.Series) -> pd.Series:
        return (col - col.mean()) / (col.std(ddof=ddof_std) + eps)

    Z = pd.DataFrame(
        {
            "zero_frac": _zcol(S["zero_frac"]),
            "cv": _zcol(S["cv"]),
            "skew": _zcol(S["skew"]).abs(),
            "kurtosis": _zcol(S["kurtosis"]),
            "tail_ratio": _zcol(S["tail_ratio"]),
        }
    )

    # Clip and sum to get the composite score
    Zc = Z.clip(lower=-clip_z, upper=clip_z)
    weird_score = Zc.sum(axis=1)

    S = S.copy()
    S["weird_score"] = weird_score

    # Sort by composite score (desc = "weirder" first)
    S = S.sort_values("weird_score", ascending=False)

    return S, Z


def _panel(values: pd.Series, ax: plt.Axes, title: str, bins: int, use_log1p: bool):
    x = values.to_numpy(dtype=float)
    x = np.log1p(x) if use_log1p else x
    x = x[np.isfinite(x)]
    ax.hist(x, bins=bins, density=True, alpha=0.45)
    try:
        xs = np.linspace(x.min(), x.max(), 400)
        if use_log1p:
            kde = gaussian_kde(x)  # already log1p-transformed
        else:
            kde = gaussian_kde(x[x > 0] if np.count_nonzero(x) > 5 else x)
        ax.plot(xs, kde(xs), linewidth=1.75)
    except Exception:
        pass
    ax.set_title(f"{title}{' (log1p)' if use_log1p else ''}", fontsize=10)
    ax.set_xlabel("log1p(expression)" if use_log1p else "expression")
    ax.set_ylabel("density")


def plot_weird_gene_panels(
    df: pd.DataFrame,
    genes: Sequence[str],
    *,
    bins: int = 50,
    cols: int = 4,
    linear_title: str = "Weirdest genes - linear scale",
    log_title: str = "Weirdest genes - log1p scale",
    show: bool = True,
    save_linear_path: str | None = None,
    save_log_path: str | None = None,
) -> tuple[plt.Figure, plt.Figure]:
    """
    Plot histogram + KDE panels for a list of genes, both linear and log1p scales.

    Returns:
        (fig_linear, fig_log): matplotlib Figure objects.
    """
    n = len(genes)
    rows = (n + cols - 1) // cols
    figsize = (4 * cols, 3.2 * rows)

    # Linear panels
    fig_lin, axes_lin = plt.subplots(rows, cols, figsize=figsize)
    axes_lin = np.array(axes_lin).ravel()
    for ax, g in zip(axes_lin, genes, strict=False):
        _panel(df[g], ax, g, bins=bins, use_log1p=False)
    for k in range(n, len(axes_lin)):
        axes_lin[k].axis("off")
    fig_lin.suptitle(linear_title, y=1.02, fontsize=14)
    plt.tight_layout()
    if save_linear_path:
        fig_lin.savefig(save_linear_path, dpi=150, bbox_inches="tight")

    # Log1p panels
    fig_log, axes_log = plt.subplots(rows, cols, figsize=figsize)
    axes_log = np.array(axes_log).ravel()
    for ax, g in zip(axes_log, genes, strict=False):
        _panel(df[g], ax, g, bins=bins, use_log1p=True)
    for k in range(n, len(axes_log)):
        axes_log[k].axis("off")
    fig_log.suptitle(log_title, y=1.02, fontsize=14)
    plt.tight_layout()
    if save_log_path:
        fig_log.savefig(save_log_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig_lin)
        plt.close(fig_log)

    return fig_lin, fig_log
