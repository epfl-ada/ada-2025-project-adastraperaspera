from __future__ import annotations

import base64
from collections.abc import Sequence
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde
from shapely.affinity import rotate as shp_rotate
from shapely.geometry import MultiPolygon, Polygon
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


# Get current directory
vis_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.dirname(vis_dir)
src_dir = os.path.dirname(scripts_dir)
figures_dir = os.path.join(src_dir, "data", "figures")


def interactive_comp_pig_regression_grid(agg, PIGS, bin_order,figures_dir):
    # Create a 4x4 grid of subplots for up to 16 genes
    n = len(PIGS)
    rows = 4
    cols = 4

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=PIGS,
        shared_xaxes=False,
        shared_yaxes=False,
        horizontal_spacing=0.06,
        vertical_spacing=0.08,
    )

    btypes = agg["broad_type"].unique().tolist()

    for i, gene in enumerate(PIGS):
        row = i // cols + 1
        col = i % cols + 1

        sub = agg[agg["gene"] == gene]
        for bt in btypes:
            dsub = sub[sub["broad_type"] == bt]
            if dsub.empty:
                continue

            xcats = dsub["distance_bin"].astype(str)

            # Add error band
            fig.add_trace(
                go.Scatter(
                    x=pd.concat([xcats, xcats[::-1]]),
                    y=pd.concat(
                        [
                            dsub["mean"] + 1.96 * dsub["sem"],
                            (dsub["mean"] - 1.96 * dsub["sem"])[::-1],
                        ]
                    ),
                    mode="lines",
                    fill="toself",
                    line=dict(width=0),
                    name=f"{bt} ± 1.96×SEM",
                    legendgroup=bt,
                    showlegend=False,
                    hoverinfo="skip",
                ),
                row=row,
                col=col,
            )

            # Add mean line
            fig.add_trace(
                go.Scatter(
                    x=xcats,
                    y=dsub["mean"],
                    mode="lines+markers",
                    name=bt,
                    legendgroup=bt,
                    showlegend=(i == 0),
                    hovertemplate=(
                        f"Gene: {gene}<br>Type: {bt}<br>Bin: %{{x}}<br>"
                        "Mean expr: %{y:.3f}<extra></extra>"
                    ),
                ),
                row=row,
                col=col,
            )

        fig.update_xaxes(
            title_text="Distance bin",
            row=row,
            col=col,
            categoryorder="array",
            categoryarray=bin_order,
        )
        fig.update_yaxes(title_text="Mean expression", row=row, col=col)

    fig.update_layout(
        height=1500,
        width=3000,
        title_text="PIG Expression by Distance and Cell Type (4×4 grid)",
        showlegend=False,
        margin=dict(t=120, b=60, l=60, r=60),

        # transparent background (ONLY change)
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.show()

    out_path = os.path.join(figures_dir, "interactive_PIG_by_broad_type.html")
    fig.write_html(
        out_path,
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False},
        default_width="100%",
        default_height="100%",
    )

    logging.info(f"Saved to {out_path}")



def plot_gene_trends_interactive(
    mean_expr: pd.DataFrame,
    genes: list[str],
    ylabel: str = "Mean expression (log1p normalized)",
    xlabel: str = "Distance to plaque (µm, binned)",
    title: str = "Spatial gene expression gradients",
    line_width: int = 2,
    height: int = 480,
    width: int = 820,
    use_webgl: bool = True,
    *,
    sem_expr: pd.DataFrame | None = None,  # optional SEM matrix (same shape as mean_expr)
) -> go.Figure:
    """
    Interactive version of 'plot_gene_trends' with optional SEM error bars.

    mean_expr: wide matrix (rows=bins, cols=genes) of means.
    sem_expr:  wide matrix (rows=bins, cols=genes) of SEM (optional).
    genes: list of genes to display (subset of columns in mean_expr/sem_expr).
    """
    if mean_expr.empty:
        raise ValueError("mean_expr is empty; check your inputs")

    # ---- Order bins & rounded labels ----
    idx = mean_expr.index
    if isinstance(idx, pd.IntervalIndex):
        # sort by bin midpoints, then rebuild mean/sem accordingly
        order = sorted(idx, key=lambda iv: iv.mid)
        M = mean_expr.loc[order]
        S = sem_expr.loc[order] if sem_expr is not None else None
        bin_labels = [f"{int(round(iv.left))}-{int(round(iv.right))}" for iv in order]
    else:
        # keep existing order; still provide nice labels
        M = mean_expr.copy()
        S = sem_expr.copy() if sem_expr is not None else None
        bin_labels = [str(x) for x in M.index]

    # ---- Filter genes present ----
    genes_present = [g for g in genes if g in M.columns]
    if not genes_present:
        raise ValueError("None of the requested genes were found in 'mean_expr' columns.")
    if S is not None:
        # ensure SEM has the same genes; drop missing gracefully
        genes_present = [g for g in genes_present if g in S.columns]
        if not genes_present:
            raise ValueError("Requested genes not present in both mean_expr and sem_expr.")

    # ---- Build long frames for Plotly ----
    M2 = M[genes_present].copy()
    M2["__bin__"] = bin_labels
    long_mean = M2.melt(id_vars="__bin__", var_name="gene", value_name="mean_expr")

    if S is not None:
        S2 = S[genes_present].copy()
        S2["__bin__"] = bin_labels
        long_sem = S2.melt(id_vars="__bin__", var_name="gene", value_name="sem_expr")
        long = pd.merge(long_mean, long_sem, on=["__bin__", "gene"], how="left")
    else:
        long = long_mean
        long["sem_expr"] = None

    # ---- Choose trace type (error bars not supported in Scattergl) ----
    use_gl = bool(use_webgl and S is None)

    fig = go.Figure()
    for g in genes_present:
        sub = long[long["gene"] == g]
        # error bars if SEM exists
        err = None
        if S is not None:
            # Plotly expects 'array' (absolute size of +/- error)
            err = dict(
                type="data",
                array=sub["sem_expr"].to_numpy(),
                visible=True,
                thickness=1.2,
                width=3,
            )

        trace_cls = go.Scattergl if use_gl else go.Scatter
        fig.add_trace(
            trace_cls(
                x=sub["__bin__"],
                y=sub["mean_expr"],
                mode="lines+markers",
                name=g,
                line=dict(width=line_width),
                error_y=err,  # ignored by Scattergl; shown by Scatter
                hovertemplate=(
                    "Bin: %{x}<br>"
                    f"Gene: {g}<br>"
                    "Mean: %{y:.3f}"
                    + ("<br>SEM: %{customdata:.3f}" if S is not None else "")
                    + "<extra></extra>"
                ),
                customdata=sub["sem_expr"] if S is not None else None,
            )
        )

    # ---- Layout cosmetics ----
    fig.update_layout(
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        template="simple_white",
        height=height,
        width=width,
        legend_title="Gene",
        margin=dict(l=60, r=20, t=60, b=60),
    )
    # Keep the bin order as given
    fig.update_xaxes(tickangle=45, categoryorder="array", categoryarray=bin_labels)

    return fig


def plot_mean_heatmap_interactive(
    mean_expr: pd.DataFrame,
    top_n: int = 25,
    zscore: bool = True,
    title: str = "Top genes varying with plaque distance",
    height: int = 600,
    width: int = 900,
):
    """
    Interactive heatmap of top N genes with highest spatial variation.
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
    if not gene_cols:
        raise ValueError("No gene columns detected in 'mean_expr'.")

    expr = mean_expr[gene_cols].copy()

    # Gene selection by gradient
    grad = expr.diff().abs().sum().sort_values(ascending=False)
    top_genes = grad.head(top_n).index
    sub_df = expr[top_genes]

    # Z-score (optional)
    if zscore:
        # convertit sparse -> dense si besoin, puis zscore colonne par colonne
        sub_df = sub_df.apply(lambda x: x.sparse.to_dense() if pd.api.types.is_sparse(x) else x)
        sub_df = (sub_df - sub_df.mean()) / (sub_df.std(ddof=0).replace(0, np.nan))

    # Labels d'axe X
    if isinstance(mean_expr.index, pd.IntervalIndex):
        x_labels = [f"{b.left:.0f}-{b.right:.0f}" for b in mean_expr.index]
    else:
        x_labels = [str(x) for x in mean_expr.index]

    # 2D array for plotly
    fig = px.imshow(
        sub_df.T.values,
        x=x_labels,
        y=sub_df.columns,
        color_continuous_scale="RdBu" if zscore else "Magma",
        origin="upper",
        aspect="auto",
        labels=dict(color="Z-score" if zscore else "Mean log1p"),
        title=title,
        height=height,
        width=width,
    )
    # Center color scale at 0 for z score
    if zscore:
        vmax = float(np.nanmax(np.abs(sub_df.values)))
        fig.update_coloraxes(cmid=0.0, cmax=vmax, cmin=-vmax)

    fig.update_layout(
        template="simple_white",
        xaxis_title="Distance bin",
        yaxis_title="Gene",
        margin=dict(l=60, r=20, t=60, b=60),
    )
    return fig


def gene_distribution_selector_interactive(
    df: pd.DataFrame,
    genes: Sequence[str],
    bins: int = 50,
    title: str = "Per-gene distributions (histogram + KDE)",
    width: int = 900,
    height: int = 600,
    kde_points: int = 400,
):
    """
    Create an interactive HTML figure with dropdowns to choose:
      - the gene
      - the scale (raw vs log1p)

    df: DataFrame (cells x genes) with numeric columns for gene expression.
    genes: list of gene column names in df.
    """
    # --- Precompute per-gene traces (raw + log1p) ---
    traces = []  # list of go.Scatter / go.Histogram
    vis_map = {}  # (gene, scale) -> list of trace indices to set visible=True
    x_ranges = {}  # scale -> (global_min, global_max) for consistent axes

    # Prepare global x-limits for stability across genes
    def _clean(x):
        x = np.asarray(x, dtype=float)
        return x[np.isfinite(x)]

    raw_vals = _clean(pd.concat([df[g] for g in genes], axis=0, ignore_index=True))
    log_vals = _clean(np.log1p(raw_vals))
    x_ranges["raw"] = (float(np.nanmin(raw_vals)), float(np.nanmax(raw_vals)))
    x_ranges["log1p"] = (float(np.nanmin(log_vals)), float(np.nanmax(log_vals)))

    for g in genes:
        # Raw
        x_raw = _clean(df[g].to_numpy())
        # Log1p
        x_log = _clean(np.log1p(df[g].to_numpy()))

        # Histogram (raw)
        h_raw = go.Histogram(
            x=x_raw,
            nbinsx=bins,
            histnorm="probability density",
            name=f"{g} - hist (raw)",
            opacity=0.45,
            showlegend=False,
        )
        # KDE (raw)
        raw_kde_trace = None
        if x_raw.size > 5:
            xr = np.linspace(
                max(x_ranges["raw"][0], np.min(x_raw)),
                min(x_ranges["raw"][1], np.max(x_raw)),
                kde_points,
            )
            try:
                kde = gaussian_kde(x_raw[x_raw > 0] if (x_raw > 0).sum() > 5 else x_raw)
                yr = kde(xr)
                raw_kde_trace = go.Scatter(
                    x=xr,
                    y=yr,
                    mode="lines",
                    name=f"{g} - kde (raw)",
                    line=dict(width=2),
                    showlegend=False,
                )
            except Exception:
                pass

        # Histogram (log1p)
        h_log = go.Histogram(
            x=x_log,
            nbinsx=bins,
            histnorm="probability density",
            name=f"{g} - hist (log1p)",
            opacity=0.45,
            showlegend=False,
        )
        # KDE (log1p)
        log_kde_trace = None
        if x_log.size > 5:
            xl = np.linspace(
                max(x_ranges["log1p"][0], np.min(x_log)),
                min(x_ranges["log1p"][1], np.max(x_log)),
                kde_points,
            )
            try:
                kde_l = gaussian_kde(x_log)  # déjà > 0
                yl = kde_l(xl)
                log_kde_trace = go.Scatter(
                    x=xl,
                    y=yl,
                    mode="lines",
                    name=f"{g} - kde (log1p)",
                    line=dict(width=2),
                    showlegend=False,
                )
            except Exception:
                pass

        # Store indices for visibility toggling
        start_idx = len(traces)
        g_raw_idxs = [start_idx]  # raw hist
        traces.append(h_raw)
        if raw_kde_trace is not None:
            g_raw_idxs.append(len(traces))
            traces.append(raw_kde_trace)

        g_log_idxs = [len(traces)]  # log hist
        traces.append(h_log)
        if log_kde_trace is not None:
            g_log_idxs.append(len(traces))
            traces.append(log_kde_trace)

        vis_map[(g, "raw")] = g_raw_idxs
        vis_map[(g, "log1p")] = g_log_idxs

    # --- Build the figure with all traces (initially hide everything) ---
    fig = go.Figure(data=traces)
    for t in fig.data:
        t.visible = False

    # Initial state
    init_gene = genes[0]
    init_scale = "log1p"
    for idx in vis_map[(init_gene, init_scale)]:
        fig.data[idx].visible = True

    # --- Dropdowns ---
    # Helper to build visibility masks
    def visibility_for(g, scale):
        vis = [False] * len(traces)
        for idx in vis_map[(g, scale)]:
            vis[idx] = True
        return vis

    # Buttons for genes
    gene_buttons = []
    for g in genes:
        gene_buttons.append(
            dict(
                label=g,
                method="update",
                args=[
                    {"visible": visibility_for(g, init_scale)},
                    {
                        "title": f"{title} - {g} ({init_scale})",
                        "xaxis": (
                            {"title": "log1p(expression)"}
                            if init_scale == "log1p"
                            else {"title": "expression"}
                        ),
                    },
                ],
            )
        )

    # Buttons for scale
    scale_buttons = []
    for sc in ["raw", "log1p"]:
        scale_buttons.append(
            dict(
                label=sc,
                method="update",
                args=[
                    {"visible": visibility_for(init_gene, sc)},
                    {
                        "title": f"{title} - {init_gene} ({sc})",
                        "xaxis": (
                            {"title": "log1p(expression)"}
                            if sc == "log1p"
                            else {"title": "expression"}
                        ),
                    },
                ],
            )
        )

    fig.update_layout(
        width=width,
        height=height,
        template="simple_white",
        title=f"{title} - {init_gene} ({init_scale})",
        xaxis_title="log1p(expression)",
        yaxis_title="density",
        barmode="overlay",
        legend_title=None,
        updatemenus=[
            dict(
                buttons=gene_buttons,
                direction="down",
                showactive=True,
                x=0.02,
                xanchor="left",
                y=1.15,
                yanchor="top",
                bgcolor="white",
                bordercolor="#ccc",
            ),
            dict(
                buttons=scale_buttons,
                direction="down",
                showactive=True,
                x=0.30,
                xanchor="left",
                y=1.15,
                yanchor="top",
                bgcolor="white",
                bordercolor="#ccc",
            ),
        ],
        margin=dict(l=60, r=20, t=90, b=60),
    )

    # Consistent x ranges per scale (switch via relayout on button click)
    # We’ll attach ranges to layout meta for clarity (optional)
    fig.layout.meta = dict(xrange_raw=x_ranges["raw"], xrange_log=x_ranges["log1p"])
    fig.update_layout(
        title={
            "text": f"{title} - {init_gene} ({init_scale})",
            "x": 0.5,  # center horizontally
            "xanchor": "center",
            "y": 0.97,  # slightly below top edge
            "yanchor": "top",
        },
        updatemenus=[
            dict(
                buttons=gene_buttons,
                direction="down",
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.12,
                yanchor="top",  # a bit below title
                bgcolor="white",
                bordercolor="#ccc",
            ),
            dict(
                buttons=scale_buttons,
                direction="down",
                showactive=True,
                x=0.25,
                xanchor="left",
                y=1.12,
                yanchor="top",
                bgcolor="white",
                bordercolor="#ccc",
            ),
        ],
        margin=dict(l=60, r=20, t=100, b=60),
    )

    return fig


def plot_top_spatial_genes_interactive(
    stats_df: pd.DataFrame,
    top_n: int = 20,
    metrics: list[str] | None = None,  # ex: ["spearman_r","slope"]
    gene_col: str = "gene",
    p_col_candidates=("p_value", "pval", "p"),
    fdr_col_candidates=("fdr", "q_value", "adj_p", "qval"),
    title: str = "Top 20 genes by spatial metric",
    width: int = 820,
    height: int = 650,
) -> go.Figure:
    """
    Interactive horizontal bar plots of the genes most spatially associated.
    - No slider: top_n is fixed.
    - A single dropdown menu to select the metric.
    - Centered title, with the menu positioned just below the title (no overlap).
    """
    df = stats_df.copy()
    if gene_col not in df.columns:
        raise ValueError(f"'{gene_col}' absent de stats_df.")

    # Auto-detect metrics if not provided
    if metrics is None:
        non_metric = {gene_col, *p_col_candidates, *fdr_col_candidates}
        metrics = [c for c in df.select_dtypes(include=[np.number]).columns if c not in non_metric]
    if not metrics:
        raise ValueError("Aucune métrique numérique détectée. Fournis `metrics=[...]`.")

    # p/FDR columns for hover (if present)
    p_col = next((c for c in p_col_candidates if c in df.columns), None)
    fdr_col = next((c for c in fdr_col_candidates if c in df.columns), None)

    # Build one trace (horizontal bar) per metric; visibility will be toggled via dropdown
    traces = []
    vis_map = {}
    top_n = int(min(top_n, len(df)))

    for m in metrics:
        if m not in df.columns:
            continue
        top = (
            df.nlargest(top_n, m).copy().sort_values(m, ascending=True)
        )  # pour empiler vers le haut

        hover = f"<b>%{{y}}</b><br>{m}: %{{x:.4g}}"
        custom = None
        if p_col or fdr_col:
            hover += f"<br>{p_col or 'p'}: %{{customdata[0]:.2e}}" if p_col else ""
            hover += f"<br>{fdr_col or 'FDR'}: %{{customdata[1]:.2e}}" if fdr_col else ""
            custom = np.stack(
                [
                    top[p_col].to_numpy() if p_col else np.full(len(top), np.nan),
                    top[fdr_col].to_numpy() if fdr_col else np.full(len(top), np.nan),
                ],
                axis=1,
            )

        trace = go.Bar(
            x=top[m].to_numpy(),
            y=top[gene_col].to_numpy(),
            orientation="h",
            name=m,
            marker=dict(color=np.where(top[m] >= 0, "rgb(31,120,180)", "rgb(227,26,28)")),
            hovertemplate=hover,
            customdata=custom,
            visible=False,
        )
        traces.append(trace)
        vis_map[m] = len(traces) - 1

    if not traces:
        raise ValueError("Aucune trace créée (vérifie les noms de métriques).")

    init_metric = metrics[0]
    traces[vis_map[init_metric]].visible = True

    fig = go.Figure(data=traces)

    # Menu dropdown
    buttons = []
    for m in metrics:
        vis = [i == vis_map[m] for i in range(len(traces))]
        buttons.append(
            dict(
                label=m,
                method="update",
                args=[
                    {"visible": vis},
                    {"title": f"{title} - {m}", "xaxis": {"title": m}},
                ],
            )
        )

    fig.update_layout(
        width=width,
        height=height,
        template="simple_white",
        title=dict(
            text=f"{title} - {init_metric}",
            x=0.5,
            xanchor="center",
            y=0.96,
            yanchor="top",
        ),
        xaxis_title=init_metric,
        yaxis_title="Gene",
        margin=dict(l=140, r=30, t=110, b=50),  # top ↑ for menu
        showlegend=False,
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=0.02,
                xanchor="left",
                y=1.10,
                yanchor="top",
                bgcolor="white",
                bordercolor="#ccc",
                pad={"r": 6, "t": 6},
            )
        ],
    )
    return fig


def interactive_comp_pig_regression(agg, PIGS, bin_order):
    # --- mean and SEM per (gene, broad_type, bin)

    btypes = agg["broad_type"].unique().tolist()
    fig = go.Figure()

    # We’ll add 2 traces per broad type per gene: line (mean) + band (± 1.96×SEM).
    # We'll store indices of the "line" traces per gene to control legend visibility per dropdown.
    traces_per_gene = 2 * len(btypes)  # band + line for each type (upper bound)
    line_idxs_per_gene: list[list[int]] = []

    for gi, gene in enumerate(PIGS):
        sub = agg[agg["gene"] == gene]
        line_idxs_for_gene: list[int] = []
        for bi, bt in enumerate(btypes):
            dsub = sub[sub["broad_type"] == bt]
            if dsub.empty:
                # Skip empty groups (avoid blank traces)
                continue

            xcats = dsub["distance_bin"].astype(str)

            # Error band (invisible in legend) — use ~95% CI via 1.96×SEM
            fig.add_trace(
                go.Scatter(
                    x=pd.concat([xcats, xcats[::-1]]),
                    y=pd.concat(
                        [
                            dsub["mean"] + 1.96 * dsub["sem"],
                            (dsub["mean"] - 1.96 * dsub["sem"])[::-1],
                        ]
                    ),
                    mode="lines",
                    fill="toself",
                    line=dict(width=0),
                    name=f"{bt} ± 1.96×SEM",
                    legendgroup=bt,
                    showlegend=False,  # <-- avoid legend clutter
                    visible=(gi == 0),
                    hoverinfo="skip",
                )
            )

            # Mean line
            fig.add_trace(
                go.Scatter(
                    x=xcats,
                    y=dsub["mean"],
                    mode="lines+markers",
                    name=bt,
                    legendgroup=bt,
                    showlegend=(gi == 0),  # <-- legend only shown for first gene
                    visible=(gi == 0),
                    hovertemplate=(
                        f"Gene: {gene}<br>Type: {bt}<br>Bin: %{{x}}<br>"
                        "Mean expr: %{y:.3f}<extra></extra>"
                    ),
                )
            )
            # Record index of the just-added line trace for legend control
            line_idxs_for_gene.append(len(fig.data) - 1)

        line_idxs_per_gene.append(line_idxs_for_gene)

    # Dropdown buttons (toggle visibility blocks per gene)
    buttons = []
    total_traces = len(fig.data)
    for gi, gene in enumerate(PIGS):
        vis = [False] * total_traces
        showlegend = [False] * total_traces
        # Turn on visibility for traces that belong to this gene block
        # We added traces sequentially per gene; however, some may be missing.
        start = gi * traces_per_gene
        for idx in range(traces_per_gene):
            k = start + idx
            if k < total_traces:
                vis[k] = True
        # Ensure legend entries for the visible gene's line traces
        for k in line_idxs_per_gene[gi]:
            showlegend[k] = True
        buttons.append(
            dict(
                label=gene,
                method="update",
                args=[
                    {"visible": vis, "showlegend": showlegend},
                    {"title": f"{gene} expression by distance and cell type"},
                ],
            )
        )

    fig.update_layout(
        title={
            "text": "PIG expression by distance and cell type",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title="Distance to plaque (µm, binned)",
        yaxis_title="Mean expression (log1p)",
        legend_title="Cell family",
        xaxis=dict(categoryorder="array", categoryarray=bin_order),
        updatemenus=[
            dict(
                type="dropdown",
                buttons=buttons,
                x=1.1,  # --> move dropdown to the right
                xanchor="right",
                y=1.22,  # --> raise dropdown a bit higher above the legend
                yanchor="top",
                pad=dict(l=2, r=2, t=2, b=2),
                direction="down",  # menu expands downward
                showactive=True,
            )
        ],
        margin=dict(l=60, r=20, t=120, b=60),  # extra top space to avoid overlap
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,  # --> move legend slightly up
            xanchor="center",
            x=0.5,
        ),
    )

    # Horizontal legend above plot
    fig.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))

    fig.show()
    fig.write_html(os.path.join(figures_dir, "pig_by_distance_interactive.html"))
    logging.info(f"Saved to {os.path.join(figures_dir, 'pig_by_distance_interactive.html')}")


def plot_gene_expression_by_distance_interactive(
    summary_df: pd.DataFrame,
    pig_genes: list[str],
    distance_col: str = "distance_bin",
    gene_col: str = "gene",
    mean_col: str = "mean_expr",
    sem_col: str = "sem_expr",
    title: str = "Plaque-Induced Gene Expression vs Distance (mean ± 95% CI ≈ 1.96×SEM, log1p)",
    x_label: str = "Distance to Plaque (µm, binned)",
    y_label: str = "Mean log1p Expression",
    use_ci95: bool = True,  # multiply SEM by 1.96
):
    """
    Interactive line plot with ALL genes overlaid (different colors) and CI bands.
    Saves to frontend/public/plots/PIG_expression_vs_distance.html with transparent background.
    """

    df = summary_df.copy()

    # Keep only genes found in data
    pig_genes = [g for g in pig_genes if g in df[gene_col].unique()]
    if not pig_genes:
        raise ValueError("None of the requested genes are present in summary_df.")

    # Convert interval bins to clean labels
    def _bin_label(b):
        if isinstance(b, pd.Interval):
            return f"{int(round(b.left))}-{int(round(b.right))}"
        return str(b)

    df["bin_label"] = df[distance_col].apply(_bin_label)

    # ~95% CI if requested
    scale = 1.96 if use_ci95 else 1.0

    fig = go.Figure()

    for gene in pig_genes:
        sub = df[df[gene_col] == gene].sort_values("bin_label")

        m = sub[mean_col].to_numpy(float)
        s = sub[sem_col].to_numpy(float) * scale
        x = sub["bin_label"].tolist()

        # mean line (each gene gets a different default Plotly color)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=m,
                mode="lines+markers",
                name=f"{gene}",
                line=dict(width=2),
                marker=dict(size=7),
            )
        )

        # SEM/CI (shaded band) — same trace color family will be used automatically
        fig.add_trace(
            go.Scatter(
                x=x + x[::-1],
                y=(m + s).tolist() + (m - s)[::-1].tolist(),
                fill="toself",
                line=dict(width=0),
                hoverinfo="skip",
                name=f"{gene} CI",
                opacity=0.18,
            )
        )

    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        autosize=True,
        template="plotly_white",
        margin=dict(t=120, l=60, r=20, b=60),

        # transparent background
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    out_path = Path("frontend/public/plots") / "PIG_expression_vs_distance.html"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False, "scrollZoom": True},
        default_width="100%",
        default_height="100%",
    )

    fig.show()
    return fig

def _infer_grid_shape(key_to_pos: dict[str, tuple[int, int]]) -> tuple[int, int]:
    rs = [r for r, _ in key_to_pos.values()]
    cs = [c for _, c in key_to_pos.values()]
    return (max(rs) + 1, max(cs) + 1)


def _load_png_rgb(path: Path) -> np.ndarray:
    img = Image.open(path).convert("RGB")
    return np.asarray(img)


def _crop_bottom(img: np.ndarray, keep_ratio: float = 0.9) -> np.ndarray:
    h = img.shape[0]
    keep_h = max(1, int(round(h * keep_ratio)))
    return img[:keep_h, ...]


def _pad_to_max(img: np.ndarray, target_h: int, target_w: int) -> np.ndarray:
    h, w = img.shape[:2]
    pad_h = max(0, target_h - h)
    pad_w = max(0, target_w - w)
    # center pad
    top = pad_h // 2
    bottom = pad_h - top
    left = pad_w // 2
    right = pad_w - left
    return np.pad(img, ((top, bottom), (left, right), (0, 0)), mode="constant", constant_values=255)


def _to_data_uri(img: np.ndarray) -> str:
    pil = Image.fromarray(img.astype(np.uint8))
    from io import BytesIO

    buf = BytesIO()
    pil.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def make_wt_tg_age_grid_scatter_from_csv(
    *,
    csv_paths: dict[str, str | Path],                 # {"wt2": "...csv", "tg2": "...csv", ...}
    age_map: dict[str, tuple[str, str]],              # age_label -> (wt_key, tg_key)
    title: str | None = "WT vs TG by age (interactive)",
    filename: str = "wt_tg_age_grid_scatter.html",
    out_dir: str = "frontend/public/plots",
    max_points: int | None = 150_000,
    x_candidates=("x_centroid", "x"),
    y_candidates=("y_centroid", "y"),
    reverse_y: bool = True,
    swap_xy: bool = False,

    # --- NEW: make all panels same color ---
    uniform_color: str = "gold",
    marker_size: float = 1.8,
    marker_opacity: float = 0.65,

    # --- NEW: flip horizontally only for specific samples/keys ---
    flip_h_keys: Iterable[str] = ("wt5", "tg17"),

    row_label_wt: str = "Wild type",
    row_label_tg: str = "Transgenic",
    row_label_font_size: int = 18,
) -> go.Figure:
    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)

    flip_h_keys = set(flip_h_keys)

    # ---- helpers ----
    def pick_xy(df: pd.DataFrame):
        xcol = next((c for c in x_candidates if c in df.columns), None)
        ycol = next((c for c in y_candidates if c in df.columns), None)
        if xcol is None or ycol is None:
            raise ValueError(
                f"Missing coordinates. Looked for {x_candidates} and {y_candidates}. "
                f"Columns start: {list(df.columns)[:30]}"
            )
        return xcol, ycol

    def prep_df(key: str) -> tuple[pd.Series, pd.Series]:
        df = pd.read_csv(csv_paths[key])

        if max_points is not None and len(df) > max_points:
            df = df.sample(n=max_points, random_state=0)

        xcol, ycol = pick_xy(df)
        x = df[xcol].astype(float)
        y = df[ycol].astype(float)

        if swap_xy:
            x, y = y, x

        # horizontal flip only for selected keys (mirror within that sample's bounds)
        if key in flip_h_keys:
            xmin, xmax = float(x.min()), float(x.max())
            x = (xmin + xmax) - x

        return x, y

    # fixed age order based on insertion order in age_map
    age_labels = list(age_map.keys())
    if len(age_labels) != 3:
        raise ValueError("age_map should contain exactly 3 ages for a 2×3 grid.")

    # ---- create subplot grid ----
    fig = make_subplots(
        rows=2,
        cols=3,
        column_titles=[f"{a} months" for a in age_labels],
        horizontal_spacing=0.01,
        vertical_spacing=0.05,
    )

    # Track global ranges so all panels match
    xmins, xmaxs, ymins, ymaxs = [], [], [], []

    # single marker style for everything
    marker_common = dict(
        size=marker_size,
        opacity=marker_opacity,
        color=uniform_color,
    )

    # ---- add traces ----
    for j, age in enumerate(age_labels, start=1):
        wt_key, tg_key = age_map[age]

        # WT (row 1)
        x_wt, y_wt = prep_df(wt_key)
        xmins.append(float(x_wt.min())); xmaxs.append(float(x_wt.max()))
        ymins.append(float(y_wt.min())); ymaxs.append(float(y_wt.max()))

        fig.add_trace(
            go.Scattergl(
                x=x_wt,
                y=y_wt,
                mode="markers",
                marker=marker_common,
                showlegend=False,
                hovertemplate=f"{row_label_wt}<br>Age: {age}<br>x=%{{x:.2f}}<br>y=%{{y:.2f}}<extra></extra>",
            ),
            row=1, col=j
        )

        # TG (row 2)
        x_tg, y_tg = prep_df(tg_key)
        xmins.append(float(x_tg.min())); xmaxs.append(float(x_tg.max()))
        ymins.append(float(y_tg.min())); ymaxs.append(float(y_tg.max()))

        fig.add_trace(
            go.Scattergl(
                x=x_tg,
                y=y_tg,
                mode="markers",
                marker=marker_common,
                showlegend=False,
                hovertemplate=f"{row_label_tg}<br>Age: {age}<br>x=%{{x:.2f}}<br>y=%{{y:.2f}}<extra></extra>",
            ),
            row=2, col=j
        )

    # ---- unify axes ----
    xr = [min(xmins), max(xmaxs)]
    yr = [min(ymins), max(ymaxs)]

    for r in (1, 2):
        for c in (1, 2, 3):
            fig.update_xaxes(
                row=r, col=c,
                range=xr,
                showgrid=False,
                zeroline=False,
                visible=False,
            )
            fig.update_yaxes(
                row=r, col=c,
                range=yr,
                showgrid=False,
                zeroline=False,
                visible=False,
                scaleanchor=f"x{'' if (r == 1 and c == 1) else ((r-1)*3 + c)}",
                scaleratio=1,
                autorange="reversed" if reverse_y else True,
            )

    # ---- row labels (bigger + bold) ----
    fig.update_layout(
        annotations=list(fig.layout.annotations) + [
            dict(
                text=f"<b>{row_label_wt}</b>",
                x=0.01, y=0.97,
                xref="paper", yref="paper",
                xanchor="left", yanchor="middle",
                showarrow=False,
                font=dict(size=row_label_font_size),
            ),
            dict(
                text=f"<b>{row_label_tg}</b>",
                x=0.01, y=0.50,
                xref="paper", yref="paper",
                xanchor="left", yanchor="middle",
                showarrow=False,
                font=dict(size=row_label_font_size),
            ),
        ]
    )

    # ---- overall layout + transparency ----
    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=30, r=20, t=80 if title else 40, b=30),
        dragmode="pan",
        template="plotly_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.write_html(str(out_path), include_plotlyjs="cdn")
    return fig



def make_alignment_overlay_plot(
    *,
    df_ref: pd.DataFrame,
    df_aligned: pd.DataFrame,
    ref_label: str = "TG17 (reference)",
    aligned_label: str = "TG5 aligned → TG17",
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
    swap_xy: bool = True,  # swap x and y (requested)
    reverse_y: bool = True,  # image-like orientation (optional but usually correct)
    max_points: int | None = 200_000,
    marker_size: float = 1.8,
    opacity: float = 0.6,
    title: str = "Attempted alignment: TG5 → TG17 (interactive)",
    filename: str = "tg5_to_tg17_alignment.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    def get_xy(df: pd.DataFrame):
        x = df[x_col].astype(float)
        y = df[y_col].astype(float)
        if swap_xy:
            x, y = y, x
        return x, y

    # Downsample for speed
    if max_points is not None:
        if len(df_ref) > max_points:
            df_ref = df_ref.sample(max_points, random_state=0)
        if len(df_aligned) > max_points:
            df_aligned = df_aligned.sample(max_points, random_state=0)

    x_ref, y_ref = get_xy(df_ref)
    x_aln, y_aln = get_xy(df_aligned)

    xmin = float(min(x_ref.min(), x_aln.min()))
    xmax = float(max(x_ref.max(), x_aln.max()))
    ymin = float(min(y_ref.min(), y_aln.min()))
    ymax = float(max(y_ref.max(), y_aln.max()))

    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x=x_ref,
            y=y_ref,
            mode="markers",
            name=ref_label,
            marker=dict(size=marker_size, opacity=opacity),
            hovertemplate=f"{ref_label}<br>x=%{{x:.1f}}<br>y=%{{y:.1f}}<extra></extra>",
        )
    )

    fig.add_trace(
        go.Scattergl(
            x=x_aln,
            y=y_aln,
            mode="markers",
            name=aligned_label,
            marker=dict(size=marker_size, opacity=opacity),
            hovertemplate=f"{aligned_label}<br>x=%{{x:.1f}}<br>y=%{{y:.1f}}<extra></extra>",
        )
    )

    fig.update_layout(
        title=title,
        # height=700,
        # width=None,  # responsive inside iframe
        autosize=True,
        margin=dict(l=30, r=20, t=80, b=30),
        dragmode="pan",
        template="plotly_white",
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        # transparent background
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.update_xaxes(range=[xmin, xmax], visible=False, showgrid=False, zeroline=False)
    fig.update_yaxes(
        range=[ymin, ymax],
        visible=False,
        showgrid=False,
        zeroline=False,
        scaleanchor="x",
        scaleratio=1,
        autorange="reversed" if reverse_y else True,
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    return fig


def make_plaques_detected_plotly(
    *,
    df: pd.DataFrame,
    brain_geom: Polygon | None = None,
    sample_hulls: int = 20,
    seed: int = 42,
    rotate_180: bool = True,
    rotation_origin: tuple[float, float] | str = "auto",
    title: str | None = None,
    figsize_px: tuple[int, int] = (820, 820),
    filename: str = "plaques_detected.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    """
    Interactive Plotly visualization of detected Aβ plaques.

    Color encoding
    --------------
    - Blue   : Brain ROI
    - Green  : Convex plaques
    - Red    : Non-convex plaques
    - Orange : Sampled convex hulls

    All geometries are optionally rotated by 180° for orientation consistency.
    """

    # ------------------ sanity checks ------------------
    if "geometry" not in df.columns:
        raise ValueError("df must contain a 'geometry' column with shapely objects.")
    if "is_convex" not in df.columns:
        raise ValueError("df must contain an 'is_convex' boolean column.")

    P = df.copy()

    if "plaque_id" not in P.columns:
        P["plaque_id"] = np.arange(1, len(P) + 1, dtype=int)
    if "area" not in P.columns:
        P["area"] = P["geometry"].map(lambda g: getattr(g, "area", np.nan))

    # ------------------ helpers ------------------
    def iter_polygons(g):
        if isinstance(g, Polygon):
            yield g
        elif isinstance(g, MultiPolygon):
            for sub in g.geoms:
                if isinstance(sub, Polygon):
                    yield sub

    def compute_auto_origin():
        bounds = []
        for g in P["geometry"]:
            if g is not None and hasattr(g, "bounds"):
                bounds.append(g.bounds)
        if brain_geom is not None:
            bounds.append(brain_geom.bounds)

        arr = np.array(bounds, dtype=float)
        minx, miny = arr[:, 0].min(), arr[:, 1].min()
        maxx, maxy = arr[:, 2].max(), arr[:, 3].max()
        return (minx + maxx) / 2, (miny + maxy) / 2

    if rotation_origin == "auto":
        origin = compute_auto_origin()
    else:
        origin = rotation_origin

    def rot(g):
        if rotate_180 and g is not None:
            return shp_rotate(g, 180.0, origin=origin, use_radians=False)
        return g

    fig = go.Figure()

    # ------------------ brain ROI ------------------
    if brain_geom is not None and brain_geom.is_valid:
        g = rot(brain_geom)
        x, y = map(list, g.exterior.xy)  # IMPORTANT FIX
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                name="Brain ROI",
                line=dict(color="blue", width=1.4),
                hoverinfo="skip",
            )
        )

    # ------------------ plaques ------------------
    for row in P.itertuples(index=False):
        g = getattr(row, "geometry", None)
        is_convex = bool(getattr(row, "is_convex", False))
        pid = getattr(row, "plaque_id", None)
        area = getattr(row, "area", np.nan)

        if g is None:
            continue

        g = rot(g)

        for poly in iter_polygons(g):
            if poly.is_empty or not poly.is_valid:
                continue

            x, y = map(list, poly.exterior.xy)  # IMPORTANT FIX

            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line=dict(
                        color="green" if is_convex else "red",
                        width=0.9,
                        dash="solid" if is_convex else "dash",
                    ),
                    hovertemplate=(f"Plaque {pid}<br>" f"Area: {area:,.0f} µm²<extra></extra>"),
                    showlegend=False,
                )
            )

    # ------------------ convex hull overlays ------------------
    rng = np.random.default_rng(seed)
    n = min(len(P), int(sample_hulls))
    if n > 0:
        sample = P.sample(n=n, random_state=seed)
        for row in sample.itertuples(index=False):
            g = getattr(row, "geometry", None)
            if g is None:
                continue

            hull = rot(g.convex_hull)
            x, y = map(list, hull.exterior.xy)  # IMPORTANT FIX

            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line=dict(color="orange", width=1.2, dash="dot"),
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

    # ------------------ layout ------------------
    fig.update_layout(
        title=title or "Detected Aβ plaques after normalization",
        width=figsize_px[0],
        height=figsize_px[1],
        template="simple_white",
        margin=dict(l=40, r=40, t=60, b=40),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        ),
    )

    fig.update_xaxes(
        title="X coordinate (µm)",
        showgrid=False,
        zeroline=False,
        scaleanchor="y",
        scaleratio=1,
    )
    fig.update_yaxes(
        title="Y coordinate (µm)",
        showgrid=False,
        zeroline=False,
    )

    # ------------------ save ------------------
    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_path, include_plotlyjs="cdn")

    return fig


def make_cell_to_plaque_distance_distribution_plotly(
    df: pd.DataFrame,
    *,
    column: str = "distance_to_plaque",
    prox_thresh: float = 30.0,
    distal_thresh: float = 100.0,
    n_bins: int = 60,
    clip_quantiles: tuple[float, float] | None = (0.0, 0.99),
    show_kde: bool = True,
    kde_points: int = 400,
    bandwidth: float | None = None,
    title: str = "Distribution of cell distances to nearest plaque (linear scale)",
    filename: str = "cell_to_plaque_distance_distribution.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    """
    Interactive Plotly version of the *linear-scale* distance distribution.

    - Histogram on linear scale
    - Optional smooth density curve (KDE-like Gaussian smoothing)
    - Threshold lines at prox_thresh and distal_thresh
    - Quantile lines at 5/25/50/75/95%
    - Transparent background for embedding
    """

    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")
    if prox_thresh <= 0 or distal_thresh <= prox_thresh:
        raise ValueError("Thresholds must satisfy: 0 < prox_thresh < distal_thresh.")

    x = df[column].astype(float).dropna().to_numpy()
    if x.size == 0:
        raise ValueError(f"No valid numeric values found in '{column}'.")

    # Optional clipping (useful to avoid the long tail dominating the view)
    if clip_quantiles is not None:
        loq, hiq = clip_quantiles
        lo = np.quantile(x, loq)
        hi = np.quantile(x, hiq)
        x_plot = x[(x >= lo) & (x <= hi)]
    else:
        x_plot = x

    # Quantiles (computed on *full* distribution, not clipped)
    q5, q25, q50, q75, q95 = np.percentile(x, [5, 25, 50, 75, 95])

    # Histogram (counts)
    counts, edges = np.histogram(x_plot, bins=n_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=centers,
            y=counts,
            width=(edges[1:] - edges[:-1]),
            name="Cells",
            hovertemplate="Distance: %{x:.1f} µm<br>Count: %{y}<extra></extra>",
        )
    )

    # Optional KDE-like smoothing on histogram counts (fast, no SciPy dependency)
    if show_kde and counts.sum() > 0:
        # Heuristic bandwidth in "bin units"
        if bandwidth is None:
            # Slight smoothing proportional to bins
            bandwidth = max(1.0, n_bins / 30.0)

        # Gaussian kernel in bin space
        kx = np.arange(-int(4 * bandwidth), int(4 * bandwidth) + 1)
        kernel = np.exp(-(kx**2) / (2 * bandwidth**2))
        kernel = kernel / kernel.sum()

        smooth = np.convolve(counts.astype(float), kernel, mode="same")

        fig.add_trace(
            go.Scatter(
                x=centers,
                y=smooth,
                mode="lines",
                name="Smoothed density (KDE-like)",
                hovertemplate="Distance: %{x:.1f} µm<br>Smoothed count: %{y:.1f}<extra></extra>",
            )
        )

    # Threshold lines
    for val, lab in [
        (prox_thresh, f"proximal threshold ({prox_thresh:g} µm)"),
        (distal_thresh, f"distal threshold ({distal_thresh:g} µm)"),
    ]:
        fig.add_vline(
            x=val,
            line_width=2,
            line_dash="dash",
            annotation_text=lab,
            annotation_position="top",
        )

    # Quantile lines
    for val, lab in [
        (q5, "q5"),
        (q25, "q25"),
        (q50, "median"),
        (q75, "q75"),
        (q95, "q95"),
    ]:
        fig.add_vline(
            x=float(val),
            line_width=1,
            line_dash="dot",
            annotation_text=lab,
            annotation_position="bottom",
        )

    fig.update_layout(
        title=title,
        xaxis_title="Distance to nearest plaque (µm)",
        yaxis_title="Cell count",
        bargap=0.02,
        hovermode="x unified",
        template="simple_white",
        margin=dict(l=60, r=20, t=70, b=55),
        autosize=True,
        # transparent backgrounds for embedding
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    # Save HTML responsive
    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False},
        default_width="100%",
        default_height="100%",
    )

    return fig


def make_expression_distribution_selector_plotly(
    df: pd.DataFrame,
    genes: Sequence[str],
    *,
    all_genes: Sequence[str] | None = None,
    bins: int = 50,
    use_log1p: bool = True,
    kde: bool = True,
    kde_bandwidth: float | None = None,
    clip_quantiles: tuple[float, float] = (0.0, 0.995),
    show_iqr: bool = True,
    show_mean: bool = True,
    show_zero_fraction: bool = True,
    title: str = "Gene expression distribution",
    filename: str = "expression_distribution.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    genes = [g for g in genes if g in df.columns]
    if not genes:
        raise ValueError("None of the requested genes are present in the dataframe.")

    if all_genes is None:
        all_genes = genes
    all_genes = [g for g in all_genes if g in df.columns]

    def transform(x):
        x = x.astype(float)
        if use_log1p:
            x = np.log1p(x)
        return x[np.isfinite(x)]

    # ---------- global x-range ----------
    all_vals = []
    for g in set(genes) | set(all_genes):
        all_vals.append(transform(df[g].to_numpy()))
    all_vals = np.concatenate(all_vals)

    lo = float(np.quantile(all_vals, clip_quantiles[0]))
    hi = float(np.quantile(all_vals, clip_quantiles[1]))
    if hi <= lo:
        lo, hi = float(all_vals.min()), float(all_vals.max())

    bin_edges = np.linspace(lo, hi, bins + 1)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bw = bin_edges[1] - bin_edges[0]

    # ---------- avg reference ----------
    avg_hist = None
    avg_q25 = avg_q75 = avg_mean = avg_zero_frac = None

    if all_genes:
        hists, q25s, q75s, means, zfs = [], [], [], [], []
        for g in all_genes:
            raw = df[g].to_numpy(dtype=float)
            x = transform(raw)
            if x.size == 0:
                continue
            h, _ = np.histogram(x, bins=bin_edges)
            hists.append(h)
            q25s.append(np.quantile(x, 0.25))
            q75s.append(np.quantile(x, 0.75))
            means.append(x.mean())
            zfs.append(np.mean(raw <= 0.0) if use_log1p else np.mean(raw == 0.0))

        if hists:
            avg_hist = np.mean(hists, axis=0)
            avg_q25 = float(np.mean(q25s))
            avg_q75 = float(np.mean(q75s))
            avg_mean = float(np.mean(means))
            avg_zero_frac = float(np.mean(zfs))

    # ---------- KDE-like smoothing ----------
    def smooth(counts):
        if not kde:
            return None
        bw_bins = kde_bandwidth or max(1.0, bins / 30)
        k = int(4 * bw_bins)
        xk = np.arange(-k, k + 1)
        kernel = np.exp(-(xk**2) / (2 * bw_bins**2))
        kernel /= kernel.sum()
        return np.convolve(counts, kernel, mode="same")

    avg_smooth = smooth(avg_hist) if avg_hist is not None else None

    # ---------- traces per gene ----------
    traces = []
    visibility = []

    for gi, gene in enumerate(genes):
        raw = df[gene].to_numpy(dtype=float)
        x = transform(raw)
        gene_hist, _ = np.histogram(x, bins=bin_edges)
        gene_smooth = smooth(gene_hist)

        q25, q75 = np.quantile(x, [0.25, 0.75])
        mean = x.mean()
        zero_frac = np.mean(raw <= 0.0) if use_log1p else np.mean(raw == 0.0)

        is_visible = gi == 0

        # avg IQR
        if show_iqr and avg_q25 is not None:
            traces.append(
                go.Scatter(
                    x=[avg_q25, avg_q75, avg_q75, avg_q25],
                    y=[0, 0, max(gene_hist) * 1.1, max(gene_hist) * 1.1],
                    fill="toself",
                    fillcolor="rgba(255,224,178,0.45)",
                    line=dict(width=0),
                    showlegend=False,
                    visible=is_visible,
                )
            )
            visibility.append(is_visible)

        # gene IQR
        if show_iqr:
            traces.append(
                go.Scatter(
                    x=[q25, q75, q75, q25],
                    y=[0, 0, max(gene_hist) * 1.1, max(gene_hist) * 1.1],
                    fill="toself",
                    fillcolor="rgba(179,217,255,0.45)",
                    line=dict(width=0),
                    showlegend=False,
                    visible=is_visible,
                )
            )
            visibility.append(is_visible)

        # avg histogram
        if avg_hist is not None:
            traces.append(
                go.Bar(
                    x=centers,
                    y=avg_hist,
                    name="Hist (avg)",
                    marker_color="#F4A261",
                    opacity=1.0,
                    visible=is_visible,
                )
            )
            visibility.append(is_visible)

        # gene histogram
        traces.append(
            go.Bar(
                x=centers,
                y=gene_hist,
                name="Hist (gene)",
                marker_color="#4C78A8",
                opacity=1.0,
                visible=is_visible,
            )
        )
        visibility.append(is_visible)

        # avg KDE
        if avg_smooth is not None:
            traces.append(
                go.Scatter(
                    x=centers,
                    y=avg_smooth,
                    mode="lines",
                    name="KDE (avg)",
                    line=dict(color="#F28E2B", width=2),
                    visible=is_visible,
                )
            )
            visibility.append(is_visible)

        # gene KDE
        if gene_smooth is not None:
            traces.append(
                go.Scatter(
                    x=centers,
                    y=gene_smooth,
                    mode="lines",
                    name="KDE (gene)",
                    line=dict(color="#1F77B4", width=2),
                    visible=is_visible,
                )
            )
            visibility.append(is_visible)

        # mean lines
        if show_mean:
            traces.append(
                go.Scatter(
                    x=[mean, mean],
                    y=[0, max(gene_hist) * 1.1],
                    mode="lines",
                    line=dict(color="#1F77B4", dash="dash"),
                    name="Mean (gene)",
                    visible=is_visible,
                )
            )
            visibility.append(is_visible)

            if avg_mean is not None:
                traces.append(
                    go.Scatter(
                        x=[avg_mean, avg_mean],
                        y=[0, max(gene_hist) * 1.1],
                        mode="lines",
                        line=dict(color="#F28E2B", dash="dash"),
                        name="Mean (avg)",
                        visible=is_visible,
                    )
                )
                visibility.append(is_visible)

    # ---------- dropdown ----------
    n_traces_per_gene = len(visibility) // len(genes)

    buttons = []
    for i, g in enumerate(genes):
        vis = [False] * len(traces)
        start = i * n_traces_per_gene
        for j in range(n_traces_per_gene):
            vis[start + j] = True
        buttons.append(
            dict(
                label=g,
                method="update",
                args=[{"visible": vis}, {"title": f"{title}: {g}"}],
            )
        )

    fig = go.Figure(data=traces)

    fig.update_layout(
        title=f"{title}: {genes[0]}",
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                x=0.6,
                y=1.20,
                showactive=True,
            )
        ],
        barmode="overlay",
        autosize=True,
        margin=dict(l=60, r=20, t=90, b=55),
        template="simple_white",
        xaxis=dict(title="log1p(counts)" if use_log1p else "counts", range=[lo, hi]),
        yaxis=dict(title="Count"),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False},
        default_width="100%",
        default_height="100%",
    )

    return fig

def make_joint_clustering_umap_grid_plotly(
    adata_by_mouse: Mapping[str, "AnnData"],
    order: Sequence[str],
    *,
    n_cols: int = 3,
    title: str | None = "Leiden clusters across mice",
    filename: str = "joint_clustering_umap.html",
    out_dir: str = "frontend/public/plots",
    marker_size: float = 2.0,
    marker_opacity: float = 0.75,
    max_points: int | None = 250_000,
    reverse_y: bool = False,
) -> go.Figure:
    """
    Interactive Plotly UMAP grid colored by Leiden clusters (NO LEGEND).

    - One subplot per mouse
    - Consistent cluster colors across panels
    - Transparent background, zoom/pan enabled
    """

    n = len(order)
    n_cols = max(1, int(n_cols))
    n_rows = (n + n_cols - 1) // n_cols

    # Collect global Leiden categories
    all_leiden = []
    for mouse in order:
        ad = adata_by_mouse.get(mouse)
        if ad is None:
            continue
        if "leiden" not in ad.obs:
            raise ValueError(f"{mouse}: missing 'leiden' in adata.obs")
        all_leiden.append(ad.obs["leiden"].astype(str).to_numpy())

    if not all_leiden:
        raise ValueError("No valid mice found.")

    leiden_levels = list(pd.Categorical(np.concatenate(all_leiden)).categories)

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        subplot_titles=[str(m) for m in order],
        horizontal_spacing=0.04,
        vertical_spacing=0.08,
    )

    for i, mouse in enumerate(order):
        r, c = divmod(i, n_cols)
        row, col = r + 1, c + 1

        ad = adata_by_mouse.get(mouse)
        if ad is None:
            fig.add_annotation(
                text=f"{mouse}<br>(no data)",
                x=0.5,
                y=0.5,
                xref=f"x{'' if i == 0 else i+1} domain",
                yref=f"y{'' if i == 0 else i+1} domain",
                showarrow=False,
            )
            fig.update_xaxes(visible=False, row=row, col=col)
            fig.update_yaxes(visible=False, row=row, col=col)
            continue

        umap = ad.obsm["X_umap"]
        x = umap[:, 0].astype(float)
        y = umap[:, 1].astype(float)
        labels = ad.obs["leiden"].astype(str).to_numpy()

        if max_points is not None and len(x) > max_points:
            rng = np.random.default_rng(0)
            idx = rng.choice(len(x), size=max_points, replace=False)
            x, y, labels = x[idx], y[idx], labels[idx]

        for k in leiden_levels:
            mask = labels == k
            if not np.any(mask):
                continue

            fig.add_trace(
                go.Scattergl(
                    x=x[mask],
                    y=y[mask],
                    mode="markers",
                    marker=dict(size=marker_size, opacity=marker_opacity),
                    showlegend=False,   # NO LEGEND
                    hovertemplate=(
                        f"Mouse: {mouse}<br>"
                        f"Leiden: {k}<br>"
                        "UMAP1=%{x:.3f}<br>"
                        "UMAP2=%{y:.3f}<extra></extra>"
                    ),
                ),
                row=row,
                col=col,
            )

        fig.update_xaxes(visible=False, showgrid=False, zeroline=False, row=row, col=col)
        fig.update_yaxes(
            visible=False,
            showgrid=False,
            zeroline=False,
            autorange="reversed" if reverse_y else True,
            row=row,
            col=col,
        )

    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=20, r=20, t=70 if title else 30, b=20),
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        dragmode="pan",
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False, "scrollZoom": True},
        default_width="100%",
        default_height="100%",
    )

    return fig

def make_leiden_spatial_grid_plotly(
    df_by_mouse: Mapping[str, pd.DataFrame],
    order: Sequence[str],
    *,
    n_cols: int = 3,
    sample_for_scatter: int | None = 20_000,
    random_state: int = 0,
    title: str | None = "Spatial map of Leiden clusters across mice",
    filename: str = "joint_clustering_overlayed.html",
    out_dir: str = "frontend/public/plots",
    marker_size: float = 2.5,
    marker_opacity: float = 0.7,
    reverse_y: bool = True,  # matches ax.invert_yaxis()
) -> go.Figure:
    """
    Interactive Plotly spatial scatter grid colored by Leiden clusters (no legend).
    Expects columns: x_centroid, y_centroid, cluster_leiden
    """

    n = len(order)
    n_cols = max(1, int(n_cols))
    n_rows = (n + n_cols - 1) // n_cols

    required = {"x_centroid", "y_centroid", "cluster_leiden"}

    # global cluster levels -> consistent colors across panels
    all_clusters = []
    for mouse in order:
        df = df_by_mouse.get(mouse)
        if df is None or not required.issubset(df.columns):
            continue
        all_clusters.append(df["cluster_leiden"].astype(str).to_numpy())
    if not all_clusters:
        raise ValueError("No valid mice with required columns found.")
    cluster_levels = list(pd.Categorical(np.concatenate(all_clusters)).categories)

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        subplot_titles=[str(m) for m in order],
        horizontal_spacing=0.04,
        vertical_spacing=0.08,
    )

    rng = np.random.default_rng(random_state)

    for i, mouse in enumerate(order):
        r, c = divmod(i, n_cols)
        row, col = r + 1, c + 1

        df = df_by_mouse.get(mouse)
        if df is None:
            fig.add_annotation(
                text=f"{mouse}<br>(no data)",
                x=0.5, y=0.5,
                xref=f"x{'' if i == 0 else i+1} domain",
                yref=f"y{'' if i == 0 else i+1} domain",
                showarrow=False,
                font=dict(size=12),
            )
            fig.update_xaxes(visible=False, row=row, col=col)
            fig.update_yaxes(visible=False, row=row, col=col)
            continue

        missing = sorted(required - set(df.columns))
        if missing:
            fig.add_annotation(
                text=f"{mouse}<br>(missing {missing})",
                x=0.5, y=0.5,
                xref=f"x{'' if i == 0 else i+1} domain",
                yref=f"y{'' if i == 0 else i+1} domain",
                showarrow=False,
                font=dict(size=12),
            )
            fig.update_xaxes(visible=False, row=row, col=col)
            fig.update_yaxes(visible=False, row=row, col=col)
            continue

        plot_df = df
        if sample_for_scatter is not None and sample_for_scatter < len(df):
            idx = rng.choice(len(df), size=int(sample_for_scatter), replace=False)
            plot_df = df.iloc[idx]

        x = plot_df["x_centroid"].astype(float).to_numpy()
        y = plot_df["y_centroid"].astype(float).to_numpy()
        cl = plot_df["cluster_leiden"].astype(str).to_numpy()

        # add one trace per cluster for consistent coloring (legend disabled)
        for k in cluster_levels:
            msk = (cl == k)
            if not np.any(msk):
                continue
            fig.add_trace(
                go.Scattergl(
                    x=x[msk],
                    y=y[msk],
                    mode="markers",
                    marker=dict(size=marker_size, opacity=marker_opacity),
                    showlegend=False,
                    hovertemplate=(
                        f"Mouse: {mouse}<br>"
                        f"Leiden: {k}<br>"
                        "x=%{x:.1f} µm<br>"
                        "y=%{y:.1f} µm<extra></extra>"
                    ),
                ),
                row=row, col=col,
            )

        # axes style (like your seaborn version)
        fig.update_xaxes(
            title_text="X (µm)" if row == n_rows else "",
            showgrid=False,
            zeroline=False,
            visible=False,              # hide ticks like your other plots
            row=row, col=col,
        )
        fig.update_yaxes(
            title_text="Y (µm)" if col == 1 else "",
            showgrid=False,
            zeroline=False,
            visible=False,
            autorange="reversed" if reverse_y else True,
            scaleanchor=f"x{'' if (row == 1 and col == 1) else (i+1)}",
            scaleratio=1,
            row=row, col=col,
        )

    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=20, r=20, t=70 if title else 30, b=20),
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        dragmode="pan",
        showlegend=False,
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False, "scrollZoom": True},
        default_width="100%",
        default_height="100%",
    )

    return fig

def make_cluster_frequency_vs_distance_plotly(
    freq_df: pd.DataFrame,
    logit_df: pd.DataFrame,
    *,
    p_adj_thresh: float = 0.01,
    title: str = "Cluster frequency (%) vs. distance to nearest plaque",
    filename: str = "cluster_frequency_distance_to_plaque.html",
    out_dir: str = "frontend/public/plots",
    marker_size: int = 7,
    line_width: int = 2,
    opacity: float = 0.75,
) -> go.Figure:
    """
    Interactive Plotly version of:
      sns.lineplot(x='bin_mid', y='pct', hue='cluster_leiden')

    Uses ONLY outputs already produced by analyze_leiden_spatial().
    """

    if freq_df.empty:
        raise ValueError("freq_df is empty.")

    required = {"bin_mid", "pct", "cluster_leiden"}
    missing = required - set(freq_df.columns)
    if missing:
        raise ValueError(f"freq_df missing required columns: {sorted(missing)}")

    dfp = freq_df.copy()
    dfp["bin_mid"] = pd.to_numeric(dfp["bin_mid"], errors="coerce")
    dfp["pct"] = pd.to_numeric(dfp["pct"], errors="coerce")
    dfp["cluster_leiden"] = dfp["cluster_leiden"].astype(str)
    dfp = dfp.dropna(subset=["bin_mid", "pct", "cluster_leiden"])

    # ---------------------------------------
    # Filter to significant clusters (same as seaborn)
    # ---------------------------------------
    if logit_df is not None and not logit_df.empty and "adj_pval" in logit_df.columns:
        sig_clusters = (
            logit_df.loc[logit_df["adj_pval"] < p_adj_thresh, "cluster"]
            .astype(str)
            .unique()
        )
        dfp = dfp[dfp["cluster_leiden"].isin(sig_clusters)]

    if dfp.empty:
        raise ValueError("No data left after significance filtering.")

    dfp = dfp.sort_values(["cluster_leiden", "bin_mid"])

    fig = go.Figure()

    # one line per Leiden cluster
    for cl, sub in dfp.groupby("cluster_leiden"):
        fig.add_trace(
            go.Scatter(
                x=sub["bin_mid"],
                y=sub["pct"],
                mode="lines+markers",
                line=dict(width=line_width),
                marker=dict(size=marker_size),
                opacity=opacity,
                showlegend=False,  # ✅ no legend
                hovertemplate=(
                    f"Leiden cluster: {cl}<br>"
                    "Distance bin mid: %{x:.1f} µm<br>"
                    "Frequency: %{y:.2f}%<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=70, r=20, t=80, b=60),
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hovermode="x unified",
    )

    fig.update_xaxes(
        title="Distance to nearest plaque (µm)",
        showgrid=False,
        zeroline=False,
    )
    fig.update_yaxes(
        title="Cluster frequency (%)",
        rangemode="tozero",
        showgrid=True,
        gridcolor="rgba(0,0,0,0.08)",
        zeroline=False,
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={
            "responsive": True,
            "displayModeBar": False,
            "scrollZoom": True,
        },
        default_width="100%",
        default_height="100%",
    )

    return fig

def make_marker_enrichment_heatmap_plotly(
    expr_z: pd.DataFrame,
    *,
    title: str = "Marker gene enrichment (z-scored across clusters)",
    filename: str = "expression_per_cluster.html",
    out_dir: str = "frontend/public/plots",
    z_clip: float = 3.0,                # clip colors to [-z_clip, z_clip]
    show_values: bool = False,          # set True if you want numbers in cells
) -> go.Figure:
    """
    Interactive Plotly heatmap for expr_z from analyze_leiden_spatial().

    expr_z: DataFrame indexed by cluster, columns=genes, values=z-scores.
    """

    if expr_z is None or expr_z.empty:
        raise ValueError("expr_z is empty. (No marker genes found or enrichment skipped.)")

    # Ensure numeric matrix
    mat = expr_z.copy()
    mat = mat.apply(pd.to_numeric, errors="coerce")

    # Optional: drop columns that are all NaN
    mat = mat.loc[:, mat.notna().any(axis=0)]
    if mat.empty:
        raise ValueError("expr_z has no numeric values after cleaning.")

    # Clip to keep colormap stable
    z = mat.to_numpy(dtype=float)
    z = np.clip(z, -float(z_clip), float(z_clip))

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=mat.columns.astype(str),
            y=mat.index.astype(str),
            zmin=-float(z_clip),
            zmax=float(z_clip),
            zmid=0.0,
            colorscale="RdBu",
            colorbar=dict(title="Z-score"),
            hovertemplate=(
                "Cluster: %{y}<br>"
                "Gene: %{x}<br>"
                "Z-score: %{z:.2f}<extra></extra>"
            ),
        )
    )

    if show_values:
        # overlays text values (can get crowded if many genes)
        fig.update_traces(
            text=np.round(z, 2),
            texttemplate="%{text}",
        )

    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=120, r=30, t=80, b=80),
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.update_xaxes(
        title="",
        tickangle=-45,
        showgrid=False,
        zeroline=False,
    )
    fig.update_yaxes(
        title="",
        showgrid=False,
        zeroline=False,
        autorange="reversed",  # keeps top row at top (like seaborn)
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False},
        default_width="100%",
        default_height="100%",
    )
    return fig


def make_pig_type_spearman_heatmap_plotly(
    pig_mat: pd.DataFrame,
    prop_mat: pd.DataFrame,
    pig_cols: Sequence[str],
    *,
    fdr_alpha: float = 0.01,
    show_values: bool = True,
    title: str = "PIG correlation ↔ cellular type proportion (per bins distance)",
    filename: str = "PIG_type_spearman.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    """
    Interactive Plotly version of plot_pig_comp_heatmap().

    - Computes Spearman correlations + p-values
    - BH-FDR correction
    - Masks non-significant cells (q > fdr_alpha) as NaN
    - Displays NaNs as black (like seaborn set_bad('black'))
    """

    pig_cols = [g for g in pig_cols if g in pig_mat.columns]
    if not pig_cols:
        raise ValueError("None of pig_cols are present in pig_mat columns.")
    if prop_mat.shape[1] == 0:
        raise ValueError("prop_mat has no columns (cell types).")

    # ---- correlations + p-values ----
    corrs = pd.DataFrame(index=pig_cols, columns=prop_mat.columns, dtype=float)
    pvals = pd.DataFrame(index=pig_cols, columns=prop_mat.columns, dtype=float)

    for g in pig_cols:
        y = pd.to_numeric(pig_mat[g], errors="coerce").to_numpy(dtype=float)
        for ct in prop_mat.columns:
            x = pd.to_numeric(prop_mat[ct], errors="coerce").to_numpy(dtype=float)
            if len(y) >= 2:
                r, p = spearmanr(y, x, nan_policy="omit")
            else:
                r, p = (np.nan, np.nan)
            corrs.loc[g, ct] = float(r) if np.isfinite(r) else np.nan
            pvals.loc[g, ct] = float(p) if np.isfinite(p) else np.nan

    # ---- FDR correction ----
    mask = np.isfinite(pvals.to_numpy())
    flat = pvals.to_numpy()[mask]
    q = pvals.copy()

    if flat.size > 0:
        _, qvals, _, _ = multipletests(flat, method="fdr_bh")
        q.to_numpy()[mask] = qvals
    else:
        q[:] = np.nan

    sig = (q <= float(fdr_alpha))
    corrs_masked = corrs.where(sig)  # non-sig -> NaN

    # ---- Plotly heatmap ----
    z = corrs_masked.to_numpy(dtype=float)

    # Text annotations only for significant cells
    text = None
    if show_values:
        text = np.where(np.isfinite(z), np.round(z, 2).astype(str), "")

    # Make NaNs appear black: use a separate "background" heatmap layer in black,
    # then overlay the coolwarm heatmap with NaNs transparent.
    # Layer 1: black background
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=np.zeros_like(z),
            x=corrs_masked.columns.astype(str),
            y=corrs_masked.index.astype(str),
            colorscale=[[0, "black"], [1, "black"]],
            showscale=False,
            hoverinfo="skip",
        )
    )

    # Layer 2: coolwarm-like heatmap for significant cells only
    # (Plotly's RdBu is close; we reverse it to match coolwarm orientation)
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=corrs_masked.columns.astype(str),
            y=corrs_masked.index.astype(str),
            zmin=-1,
            zmax=1,
            zmid=0,
            colorscale="RdBu",
            reversescale=True,  # closer to seaborn coolwarm
            colorbar=dict(title="Spearman ρ"),
            text=text,
            texttemplate="%{text}" if show_values else None,
            hovertemplate=(
                "PIG: %{y}<br>"
                "Cell type: %{x}<br>"
                "Spearman ρ: %{z:.2f}<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=140, r=30, t=80, b=80),
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.update_xaxes(title="Cellular type", tickangle=-35)
    fig.update_yaxes(title="PIG Gene", autorange="reversed")

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False},
        default_width="100%",
        default_height="100%",
    )

    return fig


def make_half_distance_plotly(
    df: pd.DataFrame,
    *,
    gene_col: str = "gene",
    slope_col: str = "slope",
    qval_col: str | None = "qval",
    filter_negative: bool = True,
    sort: str = "half",  # "half" | "abs_half_desc" | "qval_then_half"
    top_n: int | None = None,
    log_scale: bool = False,
    tissue_radius_um: float | None = None,
    annotate: bool = True,
    title: str | None = None,
    filename: str = "distances_to_halve_expression.html",
    out_dir: str = "frontend/public/plots",
) -> tuple[go.Figure, pd.DataFrame]:
    """
    Interactive Plotly version of plot_half_distance().

    Returns
    -------
    fig : plotly.graph_objects.Figure
    plot_df : DataFrame indexed by gene with column 'half_dist_um'
    """

    # -----------------------------
    # Filtering
    # -----------------------------
    if filter_negative:
        work = df.loc[
            df[slope_col] < 0,
            [gene_col, slope_col] + ([qval_col] if qval_col else []),
        ].copy()
    else:
        work = df[[gene_col, slope_col] + ([qval_col] if qval_col else [])].copy()

    # -----------------------------
    # Compute half-distance
    # -----------------------------
    work["half_dist_um"] = np.log(2) / work[slope_col].abs()
    work.replace([np.inf, -np.inf], np.nan, inplace=True)
    work.dropna(subset=["half_dist_um"], inplace=True)

    # -----------------------------
    # Sorting
    # -----------------------------
    if sort == "half":
        work.sort_values("half_dist_um", ascending=True, inplace=True)
    elif sort == "abs_half_desc":
        work["abs_half"] = work["half_dist_um"].abs()
        work.sort_values("abs_half", ascending=False, inplace=True)
    elif sort == "qval_then_half":
        if qval_col is None or qval_col not in work.columns:
            raise ValueError("qval_then_half requires qval_col.")
        work.sort_values([qval_col, "half_dist_um"], ascending=[True, True], inplace=True)
    else:
        raise ValueError("Invalid sort option.")

    if top_n is not None:
        work = work.head(int(top_n))

    plot_df = work.set_index(gene_col)

    n = len(plot_df)
    if n == 0:
        raise ValueError("No rows to plot after filtering.")

    # -----------------------------
    # Plotly figure
    # -----------------------------
    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=plot_df["half_dist_um"],
            y=plot_df.index.astype(str),
            orientation="h",
            marker=dict(color="hsl(180, 60%, 35%)"),
            hovertemplate="Gene: %{y}<br>d₁/₂: %{x:,.1f} µm<extra></extra>",
        )
    )

    # Reference line for tissue radius
    if tissue_radius_um is not None:
        fig.add_vline(
            x=tissue_radius_um,
            line_dash="dash",
            line_width=1,
            annotation_text=f"radius = {tissue_radius_um:,.0f} µm",
            annotation_position="top right",
        )

    # -----------------------------
    # Layout
    # -----------------------------
    fig.update_layout(
        title=title or f"Half-distance expression for {n} genes",
        xaxis_title="Distance to halve expression (µm)",
        yaxis_title="Gene",
        yaxis=dict(autorange="reversed"),  # shortest at top
        autosize=True,
        margin=dict(l=140, r=30, t=80, b=60),
        template="simple_white",

        # transparent background
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    if log_scale:
        fig.update_xaxes(type="log")

    # -----------------------------
    # Annotations
    # -----------------------------
    if annotate:
        for gene, v in plot_df["half_dist_um"].items():
            if not np.isfinite(v):
                continue
            fig.add_annotation(
                x=v,
                y=str(gene),
                text=f"{v:,.0f} µm",
                showarrow=False,
                xanchor="left",
                xshift=6,
                font=dict(size=11),
            )

    # -----------------------------
    # Save HTML
    # -----------------------------
    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False, "scrollZoom": True},
        default_width="100%",
        default_height="100%",
    )

    fig.show()
    return fig, plot_df


def plot_model_performance_interactive(
    results_df: pd.DataFrame,
    *,
    title: str = "Model performance (R²)",
    filename: str = "model_performance_predict_dist.html",
    out_dir: str = "frontend/public/plots",
):
    """Interactive bar plot comparing R² scores across models (train vs test)."""

    results_melted = results_df.melt(
        id_vars="model",
        value_vars=["train_r2", "test_r2"],
        var_name="Dataset",
        value_name="R²",
    )

    fig = px.bar(
        results_melted,
        x="model",
        y="R²",
        color="Dataset",
        barmode="group",
        title=title,
    )

    fig.update_layout(
        autosize=True,
        template="simple_white",
        margin=dict(l=60, r=20, t=70, b=60),
        paper_bgcolor="rgba(0,0,0,0)",  # transparent background
        plot_bgcolor="rgba(0,0,0,0)",   # transparent plot area
        legend_title_text="",
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"responsive": True, "displayModeBar": False},
        default_width="100%",
        default_height="100%",
    )

    fig.show()
    return fig

