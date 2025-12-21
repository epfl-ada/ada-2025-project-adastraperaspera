from __future__ import annotations

import base64
from collections.abc import Sequence, Mapping
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
from collections.abc import Iterable
from anndata import AnnData
from typing import Any
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import to_hex


def _tab20_hex(n: int) -> list[str]:
    cmap = cm.get_cmap("tab20", max(n, 1))
    return [mcolors.to_hex(cmap(i)) for i in range(max(n, 1))]


vis_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.dirname(vis_dir)
src_dir = os.path.dirname(scripts_dir)
figures_dir = os.path.join(src_dir, "data", "figures")


def interactive_comp_pig_regression_grid(agg, PIGS, bin_order, figures_dir):

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
    sem_expr: pd.DataFrame | None = None,
) -> go.Figure:
    """
    Interactive version of 'plot_gene_trends' with optional SEM error bars.

    mean_expr: wide matrix (rows=bins, cols=genes) of means.
    sem_expr:  wide matrix (rows=bins, cols=genes) of SEM (optional).
    genes: list of genes to display (subset of columns in mean_expr/sem_expr).
    """
    if mean_expr.empty:
        raise ValueError("mean_expr is empty; verify inputs")

    idx = mean_expr.index
    if isinstance(idx, pd.IntervalIndex):

        order = sorted(idx, key=lambda iv: iv.mid)
        M = mean_expr.loc[order]
        S = sem_expr.loc[order] if sem_expr is not None else None
        bin_labels = [f"{int(round(iv.left))}-{int(round(iv.right))}" for iv in order]
    else:

        M = mean_expr.copy()
        S = sem_expr.copy() if sem_expr is not None else None
        bin_labels = [str(x) for x in M.index]

    genes_present = [g for g in genes if g in M.columns]
    if not genes_present:
        raise ValueError(
            "None of the requested genes were found in 'mean_expr' columns."
        )
    if S is not None:

        genes_present = [g for g in genes_present if g in S.columns]
        if not genes_present:
            raise ValueError(
                "Requested genes not present in both mean_expr and sem_expr."
            )

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

    use_gl = bool(use_webgl and S is None)

    fig = go.Figure()
    for g in genes_present:
        sub = long[long["gene"] == g]

        err = None
        if S is not None:

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
                error_y=err,
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

    grad = expr.diff().abs().sum().sort_values(ascending=False)
    top_genes = grad.head(top_n).index
    sub_df = expr[top_genes]

    if zscore:

        sub_df = sub_df.apply(
            lambda x: x.sparse.to_dense() if pd.api.types.is_sparse(x) else x
        )
        sub_df = (sub_df - sub_df.mean()) / (sub_df.std(ddof=0).replace(0, np.nan))

    if isinstance(mean_expr.index, pd.IntervalIndex):
        x_labels = [f"{b.left:.0f}-{b.right:.0f}" for b in mean_expr.index]
    else:
        x_labels = [str(x) for x in mean_expr.index]

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

    traces = []
    vis_map = {}
    x_ranges = {}

    def _clean(x):
        x = np.asarray(x, dtype=float)
        return x[np.isfinite(x)]

    raw_vals = _clean(pd.concat([df[g] for g in genes], axis=0, ignore_index=True))
    log_vals = _clean(np.log1p(raw_vals))
    x_ranges["raw"] = (float(np.nanmin(raw_vals)), float(np.nanmax(raw_vals)))
    x_ranges["log1p"] = (float(np.nanmin(log_vals)), float(np.nanmax(log_vals)))

    for g in genes:

        x_raw = _clean(df[g].to_numpy())

        x_log = _clean(np.log1p(df[g].to_numpy()))

        h_raw = go.Histogram(
            x=x_raw,
            nbinsx=bins,
            histnorm="probability density",
            name=f"{g} - hist (raw)",
            opacity=0.45,
            showlegend=False,
        )

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

        h_log = go.Histogram(
            x=x_log,
            nbinsx=bins,
            histnorm="probability density",
            name=f"{g} - hist (log1p)",
            opacity=0.45,
            showlegend=False,
        )

        log_kde_trace = None
        if x_log.size > 5:
            xl = np.linspace(
                max(x_ranges["log1p"][0], np.min(x_log)),
                min(x_ranges["log1p"][1], np.max(x_log)),
                kde_points,
            )
            try:
                kde_l = gaussian_kde(x_log)
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

        start_idx = len(traces)
        g_raw_idxs = [start_idx]
        traces.append(h_raw)
        if raw_kde_trace is not None:
            g_raw_idxs.append(len(traces))
            traces.append(raw_kde_trace)

        g_log_idxs = [len(traces)]
        traces.append(h_log)
        if log_kde_trace is not None:
            g_log_idxs.append(len(traces))
            traces.append(log_kde_trace)

        vis_map[(g, "raw")] = g_raw_idxs
        vis_map[(g, "log1p")] = g_log_idxs

    fig = go.Figure(data=traces)
    for t in fig.data:
        t.visible = False

    init_gene = genes[0]
    init_scale = "log1p"
    for idx in vis_map[(init_gene, init_scale)]:
        fig.data[idx].visible = True

    def visibility_for(g, scale):
        vis = [False] * len(traces)
        for idx in vis_map[(g, scale)]:
            vis[idx] = True
        return vis

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

    fig.layout.meta = dict(xrange_raw=x_ranges["raw"], xrange_log=x_ranges["log1p"])
    fig.update_layout(
        title={
            "text": f"{title} - {init_gene} ({init_scale})",
            "x": 0.5,
            "xanchor": "center",
            "y": 0.97,
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
                yanchor="top",
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
    metrics: list[str] | None = None,
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

    if metrics is None:
        non_metric = {gene_col, *p_col_candidates, *fdr_col_candidates}
        metrics = [
            c
            for c in df.select_dtypes(include=[np.number]).columns
            if c not in non_metric
        ]
    if not metrics:
        raise ValueError("Aucune métrique numérique détectée. Fournis `metrics=[...]`.")

    p_col = next((c for c in p_col_candidates if c in df.columns), None)
    fdr_col = next((c for c in fdr_col_candidates if c in df.columns), None)

    traces = []
    vis_map = {}
    top_n = int(min(top_n, len(df)))

    for m in metrics:
        if m not in df.columns:
            continue
        top = df.nlargest(top_n, m).copy().sort_values(m, ascending=True)

        hover = f"<b>%{{y}}</b><br>{m}: %{{x:.4g}}"
        custom = None
        if p_col or fdr_col:
            hover += f"<br>{p_col or 'p'}: %{{customdata[0]:.2e}}" if p_col else ""
            hover += (
                f"<br>{fdr_col or 'FDR'}: %{{customdata[1]:.2e}}" if fdr_col else ""
            )
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
            marker=dict(
                color=np.where(top[m] >= 0, "rgb(31,120,180)", "rgb(227,26,28)")
            ),
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
        margin=dict(l=140, r=30, t=110, b=50),
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

    btypes = agg["broad_type"].unique().tolist()
    fig = go.Figure()

    traces_per_gene = 2 * len(btypes)
    line_idxs_per_gene: list[list[int]] = []

    for gi, gene in enumerate(PIGS):
        sub = agg[agg["gene"] == gene]
        line_idxs_for_gene: list[int] = []
        for bi, bt in enumerate(btypes):
            dsub = sub[sub["broad_type"] == bt]
            if dsub.empty:

                continue

            xcats = dsub["distance_bin"].astype(str)

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
                    visible=(gi == 0),
                    hoverinfo="skip",
                )
            )

            fig.add_trace(
                go.Scatter(
                    x=xcats,
                    y=dsub["mean"],
                    mode="lines+markers",
                    name=bt,
                    legendgroup=bt,
                    showlegend=(gi == 0),
                    visible=(gi == 0),
                    hovertemplate=(
                        f"Gene: {gene}<br>Type: {bt}<br>Bin: %{{x}}<br>"
                        "Mean expr: %{y:.3f}<extra></extra>"
                    ),
                )
            )

            line_idxs_for_gene.append(len(fig.data) - 1)

        line_idxs_per_gene.append(line_idxs_for_gene)

    buttons = []
    total_traces = len(fig.data)
    for gi, gene in enumerate(PIGS):
        vis = [False] * total_traces
        showlegend = [False] * total_traces

        start = gi * traces_per_gene
        for idx in range(traces_per_gene):
            k = start + idx
            if k < total_traces:
                vis[k] = True

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
                x=1.1,
                xanchor="right",
                y=1.22,
                yanchor="top",
                pad=dict(l=2, r=2, t=2, b=2),
                direction="down",
                showactive=True,
            )
        ],
        margin=dict(l=60, r=20, t=120, b=60),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,
            xanchor="center",
            x=0.5,
        ),
    )

    fig.update_layout(
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0)
    )

    fig.show()
    fig.write_html(os.path.join(figures_dir, "pig_by_distance_interactive.html"))
    logging.info(
        f"Saved to {os.path.join(figures_dir, 'pig_by_distance_interactive.html')}"
    )


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
    use_ci95: bool = True,
):
    """
    Interactive line plot with ALL genes overlaid (different colors) and CI bands.
    Saves to frontend/public/plots/PIG_expression_vs_distance.html with transparent background.
    """

    df = summary_df.copy()

    pig_genes = [g for g in pig_genes if g in df[gene_col].unique()]
    if not pig_genes:
        raise ValueError("None of the requested genes are present in summary_df.")

    def _bin_label(b):
        if isinstance(b, pd.Interval):
            return f"{int(round(b.left))}-{int(round(b.right))}"
        return str(b)

    df["bin_label"] = df[distance_col].apply(_bin_label)

    scale = 1.96 if use_ci95 else 1.0

    fig = go.Figure()

    for gene in pig_genes:
        sub = df[df[gene_col] == gene].sort_values("bin_label")

        m = sub[mean_col].to_numpy(float)
        s = sub[sem_col].to_numpy(float) * scale
        x = sub["bin_label"].tolist()

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

    top = pad_h // 2
    bottom = pad_h - top
    left = pad_w // 2
    right = pad_w - left
    return np.pad(
        img,
        ((top, bottom), (left, right), (0, 0)),
        mode="constant",
        constant_values=255,
    )


def _to_data_uri(img: np.ndarray) -> str:
    pil = Image.fromarray(img.astype(np.uint8))
    from io import BytesIO

    buf = BytesIO()
    pil.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def make_wt_tg_age_grid_scatter_from_csv(
    *,
    csv_paths: dict[str, str | Path],
    age_map: dict[str, tuple[str, str]],
    title: str | None = "WT vs TG by age (interactive)",
    filename: str = "wt_tg_age_grid_scatter.html",
    out_dir: str = "frontend/public/plots",
    max_points: int | None = 150_000,
    x_candidates=("x_centroid", "x"),
    y_candidates=("y_centroid", "y"),
    reverse_y: bool = True,
    swap_xy: bool = False,
    uniform_color: str = "gold",
    marker_size: float = 1.8,
    marker_opacity: float = 0.65,
    flip_h_keys: Iterable[str] = ("wt5", "tg17"),
    row_label_wt: str = "Wild type",
    row_label_tg: str = "Transgenic",
    row_label_font_size: int = 18,
) -> go.Figure:
    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)

    flip_h_keys = set(flip_h_keys)

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

        if key in flip_h_keys:
            xmin, xmax = float(x.min()), float(x.max())
            x = (xmin + xmax) - x

        return x, y

    age_labels = list(age_map.keys())
    if len(age_labels) != 3:
        raise ValueError("age_map should contain exactly 3 ages for a 2×3 grid.")

    fig = make_subplots(
        rows=2,
        cols=3,
        column_titles=[f"{a} months" for a in age_labels],
        horizontal_spacing=0.01,
        vertical_spacing=0.05,
    )

    xmins, xmaxs, ymins, ymaxs = [], [], [], []

    marker_common = dict(
        size=marker_size,
        opacity=marker_opacity,
        color=uniform_color,
    )

    for j, age in enumerate(age_labels, start=1):
        wt_key, tg_key = age_map[age]

        x_wt, y_wt = prep_df(wt_key)
        xmins.append(float(x_wt.min()))
        xmaxs.append(float(x_wt.max()))
        ymins.append(float(y_wt.min()))
        ymaxs.append(float(y_wt.max()))

        fig.add_trace(
            go.Scattergl(
                x=x_wt,
                y=y_wt,
                mode="markers",
                marker=marker_common,
                showlegend=False,
                hovertemplate=f"{row_label_wt}<br>Age: {age}<br>x=%{{x:.2f}}<br>y=%{{y:.2f}}<extra></extra>",
            ),
            row=1,
            col=j,
        )

        x_tg, y_tg = prep_df(tg_key)
        xmins.append(float(x_tg.min()))
        xmaxs.append(float(x_tg.max()))
        ymins.append(float(y_tg.min()))
        ymaxs.append(float(y_tg.max()))

        fig.add_trace(
            go.Scattergl(
                x=x_tg,
                y=y_tg,
                mode="markers",
                marker=marker_common,
                showlegend=False,
                hovertemplate=f"{row_label_tg}<br>Age: {age}<br>x=%{{x:.2f}}<br>y=%{{y:.2f}}<extra></extra>",
            ),
            row=2,
            col=j,
        )

    xr = [min(xmins), max(xmaxs)]
    yr = [min(ymins), max(ymaxs)]

    for r in (1, 2):
        for c in (1, 2, 3):
            fig.update_xaxes(
                row=r,
                col=c,
                range=xr,
                showgrid=False,
                zeroline=False,
                visible=False,
            )
            fig.update_yaxes(
                row=r,
                col=c,
                range=yr,
                showgrid=False,
                zeroline=False,
                visible=False,
                scaleanchor=f"x{'' if (r == 1 and c == 1) else ((r-1)*3 + c)}",
                scaleratio=1,
                autorange="reversed" if reverse_y else True,
            )

    fig.update_layout(
        annotations=list(fig.layout.annotations)
        + [
            dict(
                text=f"<b>{row_label_wt}</b>",
                x=0.01,
                y=0.97,
                xref="paper",
                yref="paper",
                xanchor="left",
                yanchor="middle",
                showarrow=False,
                font=dict(size=row_label_font_size),
            ),
            dict(
                text=f"<b>{row_label_tg}</b>",
                x=0.01,
                y=0.50,
                xref="paper",
                yref="paper",
                xanchor="left",
                yanchor="middle",
                showarrow=False,
                font=dict(size=row_label_font_size),
            ),
        ]
    )

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
    swap_xy: bool = True,
    reverse_y: bool = True,
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
        autosize=True,
        margin=dict(l=30, r=20, t=80, b=30),
        dragmode="pan",
        template="plotly_white",
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
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

    if "geometry" not in df.columns:
        raise ValueError("df must contain a 'geometry' column with shapely objects.")
    if "is_convex" not in df.columns:
        raise ValueError("df must contain an 'is_convex' boolean column.")

    P = df.copy()

    if "plaque_id" not in P.columns:
        P["plaque_id"] = np.arange(1, len(P) + 1, dtype=int)
    if "area" not in P.columns:
        P["area"] = P["geometry"].map(lambda g: getattr(g, "area", np.nan))

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

    if brain_geom is not None and brain_geom.is_valid:
        g = rot(brain_geom)
        x, y = map(list, g.exterior.xy)
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

            x, y = map(list, poly.exterior.xy)

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
                    hovertemplate=(
                        f"Plaque {pid}<br>" f"Area: {area:,.0f} µm²<extra></extra>"
                    ),
                    showlegend=False,
                )
            )

    rng = np.random.default_rng(seed)
    n = min(len(P), int(sample_hulls))
    if n > 0:
        sample = P.sample(n=n, random_state=seed)
        for row in sample.itertuples(index=False):
            g = getattr(row, "geometry", None)
            if g is None:
                continue

            hull = rot(g.convex_hull)
            x, y = map(list, hull.exterior.xy)

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

    if clip_quantiles is not None:
        loq, hiq = clip_quantiles
        lo = np.quantile(x, loq)
        hi = np.quantile(x, hiq)
        x_plot = x[(x >= lo) & (x <= hi)]
    else:
        x_plot = x

    q5, q25, q50, q75, q95 = np.percentile(x, [5, 25, 50, 75, 95])

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

    if show_kde and counts.sum() > 0:

        if bandwidth is None:

            bandwidth = max(1.0, n_bins / 30.0)

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
        height=None,
        width=980,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        full_html=True,
        config={"displayModeBar": False},
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


def get_leiden_color_map(adata_by_mouse, order, key="leiden"):
    """
    Return {cluster_label: hex_color} using Scanpy's stored palette:
      adata.uns[f"{key}_colors"].
    Works even if adata.obs[key] is not categorical yet.
    """
    for mouse in order:
        ad = adata_by_mouse.get(mouse)
        if ad is None:
            continue

        if key not in ad.obs:
            continue

        if not pd.api.types.is_categorical_dtype(ad.obs[key]):
            ad.obs[key] = ad.obs[key].astype("category")

        color_key = f"{key}_colors"
        if color_key not in ad.uns:
            raise ValueError(
                f"{mouse}: missing adata.uns['{color_key}'].\n"
                f"Run once: sc.pl.umap(adata, color='{key}', show=False) "
                f"or set adata.uns['{color_key}'] manually."
            )

        levels = list(ad.obs[key].cat.categories.astype(str))
        colors = list(ad.uns[color_key])

        if len(levels) != len(colors):
            raise ValueError(
                f"{mouse}: mismatch between {key} categories ({len(levels)}) "
                f"and {color_key} ({len(colors)})."
            )

        return dict(zip(levels, colors))

    raise ValueError("Could not extract colors from any mouse in `order`.")


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
    leiden_color_map = get_leiden_color_map(adata_by_mouse, order)

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
                    marker=dict(
                        size=marker_size,
                        opacity=marker_opacity,
                        color=leiden_color_map[k],
                    ),
                    showlegend=False,
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

        fig.update_xaxes(
            visible=False, showgrid=False, zeroline=False, row=row, col=col
        )
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
        height=300 * n_rows,
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

    return fig, leiden_color_map


def make_leiden_spatial_grid_plotly(
    df_by_mouse: Mapping[str, pd.DataFrame],
    order: Sequence[str],
    *,
    n_cols: int = 3,
    sample_for_scatter: int | None = 20_000,
    random_state: int = 0,
    suptitle: str = "Spatial map of Leiden clusters across mice",
    legend_title: str = "Leiden cluster",
    palette_name: str = "tab20",
    point_size: float = 4,
    point_alpha: float = 0.7,
    invert_y: bool = True,
    figsize_px: tuple[int, int] = (1100, 750),
    filename: str | None = None,
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
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

    if palette_name != "tab20":
        raise ValueError(
            "This implementation currently matches Matplotlib 'tab20' exactly. Use palette_name='tab20'."
        )

    colors = _tab20_hex(len(hue_order))
    cluster_to_color = {cl: col for cl, col in zip(hue_order, colors)}

    plot_cols = max(1, n_cols)
    n = len(order)
    n_rows = (n + plot_cols - 1) // plot_cols

    fig = make_subplots(
        rows=n_rows,
        cols=plot_cols,
        subplot_titles=[str(m) for m in order] + [""] * (n_rows * plot_cols - n),
        horizontal_spacing=0.12,
        vertical_spacing=0.16,
    )

    for i, mouse in enumerate(order):
        r, c = divmod(i, plot_cols)
        row, col = r + 1, c + 1

        df = df_by_mouse.get(mouse)
        if df is None or not required.issubset(df.columns) or len(df) == 0:

            fig.add_annotation(
                x=0.5,
                y=0.5,
                xref=f"x{i+1} domain",
                yref=f"y{i+1} domain",
                text=f"{mouse}<br>(no data)",
                showarrow=False,
            )
            continue

        plot_df = df
        if sample_for_scatter is not None and sample_for_scatter < len(df):
            plot_df = df.sample(sample_for_scatter, random_state=random_state)

        for j, cl in enumerate(hue_order):
            sub = plot_df[plot_df["cluster_leiden"] == cl]
            if sub.empty:
                continue

            fig.add_trace(
                go.Scattergl(
                    x=sub["x_centroid"],
                    y=sub["y_centroid"],
                    mode="markers",
                    name=str(cl),
                    legendgroup=str(cl),
                    showlegend=(i == 0),
                    marker=dict(
                        size=point_size, color=cluster_to_color[cl], opacity=point_alpha
                    ),
                    hovertemplate=(
                        f"Mouse: {mouse}<br>"
                        f"Cluster: {cl}<br>"
                        "x: %{x:.1f}<br>y: %{y:.1f}<extra></extra>"
                    ),
                ),
                row=row,
                col=col,
            )

        fig.update_xaxes(
            title_text="X coordinate (µm)",
            row=row,
            col=col,
            showgrid=False,
            zeroline=False,
        )
        fig.update_yaxes(
            title_text="Y coordinate (µm)",
            row=row,
            col=col,
            showgrid=False,
            zeroline=False,
            autorange="reversed" if invert_y else True,
            scaleanchor=f"x{(i+1) if (i>0) else ''}",
            scaleratio=1,
        )

    fig.update_layout(
        title=suptitle,
        width=figsize_px[0],
        height=figsize_px[1],
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=60, r=240, t=90, b=70),
        legend=dict(
            title=legend_title,
            orientation="v",
            x=1.02,
            y=1.0,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0,0,0,0)",
            borderwidth=0,
        ),
    )

    if filename:
        out_path = Path(out_dir) / filename
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(
            str(out_path),
            include_plotlyjs="cdn",
            full_html=True,
            config={"responsive": True, "displayModeBar": False, "scrollZoom": True},
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

    for cl, sub in dfp.groupby("cluster_leiden"):
        fig.add_trace(
            go.Scatter(
                x=sub["bin_mid"],
                y=sub["pct"],
                mode="lines+markers",
                line=dict(width=line_width),
                marker=dict(size=marker_size),
                opacity=opacity,
                showlegend=False,
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
    z_clip: float = 3.0,
    show_values: bool = False,
) -> go.Figure:
    """
    Interactive Plotly heatmap for expr_z from analyze_leiden_spatial().

    expr_z: DataFrame indexed by cluster, columns=genes, values=z-scores.
    """

    if expr_z is None or expr_z.empty:
        raise ValueError(
            "expr_z is empty. (No marker genes found or enrichment skipped.)"
        )

    mat = expr_z.copy()
    mat = mat.apply(pd.to_numeric, errors="coerce")

    mat = mat.loc[:, mat.notna().any(axis=0)]
    if mat.empty:
        raise ValueError("expr_z has no numeric values after cleaning.")

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
                "Cluster: %{y}<br>" "Gene: %{x}<br>" "Z-score: %{z:.2f}<extra></extra>"
            ),
        )
    )

    if show_values:

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
        autorange="reversed",
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
    title: str = "PIG correlation <-> cellular type proportion (per bins distance)",
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

    mask = np.isfinite(pvals.to_numpy())
    flat = pvals.to_numpy()[mask]
    q = pvals.copy()

    if flat.size > 0:
        _, qvals, _, _ = multipletests(flat, method="fdr_bh")
        q.to_numpy()[mask] = qvals
    else:
        q[:] = np.nan

    sig = q <= float(fdr_alpha)
    corrs_masked = corrs.where(sig)

    z = corrs_masked.to_numpy(dtype=float)

    text = None
    if show_values:
        text = np.where(np.isfinite(z), np.round(z, 2).astype(str), "")

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

    fig.add_trace(
        go.Heatmap(
            z=z,
            x=corrs_masked.columns.astype(str),
            y=corrs_masked.index.astype(str),
            zmin=-1,
            zmax=1,
            zmid=0,
            colorscale="RdBu",
            reversescale=True,
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
    sort: str = "half",
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

    if filter_negative:
        work = df.loc[
            df[slope_col] < 0,
            [gene_col, slope_col] + ([qval_col] if qval_col else []),
        ].copy()
    else:
        work = df[[gene_col, slope_col] + ([qval_col] if qval_col else [])].copy()

    work["half_dist_um"] = np.log(2) / work[slope_col].abs()
    work.replace([np.inf, -np.inf], np.nan, inplace=True)
    work.dropna(subset=["half_dist_um"], inplace=True)

    if sort == "half":
        work.sort_values("half_dist_um", ascending=True, inplace=True)
    elif sort == "abs_half_desc":
        work["abs_half"] = work["half_dist_um"].abs()
        work.sort_values("abs_half", ascending=False, inplace=True)
    elif sort == "qval_then_half":
        if qval_col is None or qval_col not in work.columns:
            raise ValueError("qval_then_half requires qval_col.")
        work.sort_values(
            [qval_col, "half_dist_um"], ascending=[True, True], inplace=True
        )
    else:
        raise ValueError("Invalid sort option.")

    if top_n is not None:
        work = work.head(int(top_n))

    plot_df = work.set_index(gene_col)

    n = len(plot_df)
    if n == 0:
        raise ValueError("No rows to plot after filtering.")

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

    if tissue_radius_um is not None:
        fig.add_vline(
            x=tissue_radius_um,
            line_dash="dash",
            line_width=1,
            annotation_text=f"radius = {tissue_radius_um:,.0f} µm",
            annotation_position="top right",
        )

    fig.update_layout(
        title=title or f"Half-distance expression for {n} genes",
        xaxis_title="Distance to halve expression (µm)",
        yaxis_title="Gene",
        yaxis=dict(autorange="reversed"),
        autosize=True,
        margin=dict(l=140, r=30, t=80, b=60),
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    if log_scale:
        fig.update_xaxes(type="log")

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
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
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


def plot_residual_vs_distance_interactive(
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    model_name: str = "Model",
    out_dir: str = "frontend/public/plots",
    filename: str = "residuals_diagnostics.html",
) -> go.Figure:
    """
    Interactive scatter of residual (y_true - y_pred) vs true distance,
    with rolling-median smoother and Pearson r printed (and shown in title).
    Transparent background + responsive HTML export.
    """

    y_true_s = pd.to_numeric(pd.Series(y_true), errors="coerce")
    y_pred_s = pd.to_numeric(pd.Series(y_pred), errors="coerce")

    residuals = y_true_s - y_pred_s

    mask = ~(y_true_s.isna() | residuals.isna())
    y_true_s = y_true_s[mask]
    residuals = residuals[mask]

    if len(y_true_s) == 0:
        raise ValueError("No finite values for residual vs distance scatter.")

    r = float(np.corrcoef(y_true_s.to_numpy(), residuals.to_numpy())[0, 1])
    print(f"Pearson r = {r:.2f}")

    order = np.argsort(y_true_s.to_numpy())
    x_sorted = y_true_s.to_numpy()[order]
    y_sorted = residuals.to_numpy()[order]
    window = max(20, len(x_sorted) // 50)

    y_smooth = None
    if window > 5:
        y_smooth = pd.Series(y_sorted).rolling(window, center=True).median().to_numpy()

    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x=y_true_s.to_numpy(),
            y=residuals.to_numpy(),
            mode="markers",
            marker=dict(size=6, opacity=0.5),
            name=f"Cells (n={len(y_true_s)})",
            hovertemplate=(
                "True distance: %{x:.1f} µm<br>" "Residual: %{y:.1f} µm<extra></extra>"
            ),
        )
    )

    if y_smooth is not None:
        fig.add_trace(
            go.Scatter(
                x=x_sorted,
                y=y_smooth,
                mode="lines",
                line=dict(width=2, color="black"),
                name="Rolling median residual",
                hoverinfo="skip",
            )
        )

    fig.add_hline(y=0, line_dash="dash", line_width=1, line_color="red")

    fig.update_layout(
        title=f"{model_name}: Residuals vs distance (Pearson r = {r:.2f})",
        xaxis_title="True distance to plaque (µm)",
        yaxis_title="True − predicted (µm)",
        autosize=True,
        template="simple_white",
        margin=dict(l=70, r=20, t=80, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
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

    fig.show()
    return fig


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
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
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


def plot_true_pred_kde_interactive(
    y_true: pd.Series,
    y_pred: np.ndarray,
    *,
    model_name: str = "Model",
    show_qq: bool = True,
    kde_points: int = 400,
    qq_points: int | None = 200,
    qq_ref_line: bool = True,
    qq_log: bool = False,
    out_dir: str = "frontend/public/plots",
    filename: str = "true_vs_predicted.html",
):
    """
    Interactive KDE(true vs predicted) + optional Q-Q plot.
    Always stacked vertically (responsive-safe).
    """

    y_true_arr = pd.to_numeric(pd.Series(y_true), errors="coerce").to_numpy()
    y_pred_arr = pd.to_numeric(pd.Series(y_pred), errors="coerce").to_numpy()

    y_true_arr = y_true_arr[np.isfinite(y_true_arr)]
    y_pred_arr = y_pred_arr[np.isfinite(y_pred_arr)]

    if y_true_arr.size == 0 or y_pred_arr.size == 0:
        raise ValueError("y_true and y_pred must contain finite values.")

    rows = 2 if show_qq else 1
    fig = make_subplots(
        rows=rows,
        cols=1,
        vertical_spacing=0.12,
        subplot_titles=(
            ["KDE: true vs predicted", "Q–Q plot"]
            if show_qq
            else ["KDE: true vs predicted"]
        ),
    )

    common_min = float(min(y_true_arr.min(), y_pred_arr.min()))
    common_max = float(max(y_true_arr.max(), y_pred_arr.max()))
    xs = np.linspace(common_min, common_max, kde_points)

    kde_true = gaussian_kde(y_true_arr)(xs)
    kde_pred = gaussian_kde(y_pred_arr)(xs)

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=kde_true,
            mode="lines",
            name="True distance",
            line=dict(width=2),
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=kde_pred,
            mode="lines",
            name="Predicted distance",
            line=dict(width=2, dash="dash"),
        ),
        row=1,
        col=1,
    )

    fig.update_xaxes(title_text="Distance to plaque (µm)", row=1, col=1)
    fig.update_yaxes(title_text="Density", row=1, col=1)

    if show_qq:
        n = min(y_true_arr.size, y_pred_arr.size)
        if qq_points is not None:
            n = min(n, int(qq_points))
        if n < 2:
            raise ValueError("Need ≥ 2 points for Q-Q plot.")

        ps = (np.arange(1, n + 1) - 0.5) / n
        q_true = np.quantile(y_true_arr, ps)
        q_pred = np.quantile(y_pred_arr, ps)

        fig.add_trace(
            go.Scatter(
                x=q_true,
                y=q_pred,
                mode="markers",
                marker=dict(size=6, opacity=0.7),
                name="Quantiles",
            ),
            row=2,
            col=1,
        )

        if qq_ref_line:
            lo = float(min(q_true.min(), q_pred.min()))
            hi = float(max(q_true.max(), q_pred.max()))
            fig.add_trace(
                go.Scatter(
                    x=[lo, hi],
                    y=[lo, hi],
                    mode="lines",
                    line=dict(dash="dash", width=1),
                    showlegend=False,
                ),
                row=2,
                col=1,
            )

        fig.update_xaxes(title_text="True quantiles", row=2, col=1)
        fig.update_yaxes(title_text="Predicted quantiles", row=2, col=1)

        if qq_log:
            if (q_true <= 0).any() or (q_pred <= 0).any():
                raise ValueError("qq_log=True requires positive values.")
            fig.update_xaxes(type="log", row=2, col=1)
            fig.update_yaxes(type="log", row=2, col=1)

    fig.update_layout(
        height=700 if show_qq else 420,
        title_text=f"{model_name}: True vs Predicted Distance",
        template="simple_white",
        margin=dict(t=80, l=60, r=40, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
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


def _bivariate_palette_rowwise(
    *,
    n_resid: int,
    n_dist: int,
    dist_lo: float,
    dist_hi: float,
) -> np.ndarray:
    """
    Returns array shape (n_resid, n_dist, 3) in [0,1].

    Residual controls base hue row-wise (low/mid/high residual),
    distance controls brightness col-wise (near->far).
    """

    base = np.array(
        [
            [0.25, 0.55, 0.85],
            [0.55, 0.75, 0.55],
            [0.90, 0.45, 0.45],
        ],
        dtype=float,
    )[:n_resid]

    factors = (
        np.array([dist_hi]) if n_dist == 1 else np.linspace(dist_lo, dist_hi, n_dist)
    )

    out = np.zeros((n_resid, n_dist, 3), dtype=float)
    for r in range(n_resid):
        for d in range(n_dist):
            f = factors[d]
            mix = 0.35 + 0.65 * f
            out[r, d] = base[r] * mix + (1.0 - mix) * np.array([1.0, 1.0, 1.0])
            out[r, d] = np.clip(out[r, d], 0.0, 1.0)
    return out


def plot_bivariate_resid_distance_spatial_interactive(
    cells_df: pd.DataFrame,
    y_true: pd.Series,
    y_pred: np.ndarray,
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
    legend_count_fontsize: int = 11,
    legend_count_min: int | None = None,
    palette_dist_lo: float = 0.25,
    palette_dist_hi: float = 0.95,
    colors_rgb: np.ndarray | None = None,
    plaque_x_col: str = "plaque_x",
    plaque_y_col: str = "plaque_y",
    show_plaques: bool = True,
    plaque_color: str = "red",
    plaque_edgecolor: str = "red",
    plaque_size: float = 8,
    plaque_alpha: float = 0.95,
    plaque_linewidth: float = 0.8,
    plaque_round_decimals: int | None = None,
    plaque_label: str = "Plaque",
    show_marker_legend: bool = True,
    marker_legend_fontsize: int = 12,
    marker_legend_framealpha: float = 0.90,
    layout: str = "vertical",
    title: str = "Residual vs distance bivariate map",
    out_dir: str = "frontend/public/plots",
    filename: str = "residuals_vs_distance.html",
) -> go.Figure:
    if n_bins != 3:
        raise ValueError("This implementation expects n_bins=3.")

    if layout not in {"vertical", "horizontal"}:
        raise ValueError("layout must be 'vertical' or 'horizontal'")

    df = cells_df.copy()

    resid = np.asarray(y_true, dtype=float) - np.asarray(y_pred, dtype=float)
    df["_resid_"] = resid
    df["_dist_"] = pd.to_numeric(df[dist_col], errors="coerce")
    df["_x_"] = pd.to_numeric(df[x_col], errors="coerce")
    df["_y_"] = pd.to_numeric(df[y_col], errors="coerce")

    mask = (
        np.isfinite(df["_resid_"].to_numpy())
        & np.isfinite(df["_dist_"].to_numpy())
        & np.isfinite(df["_x_"].to_numpy())
        & np.isfinite(df["_y_"].to_numpy())
    )
    df = df.loc[mask].copy()

    df["_resid_bin_"] = pd.qcut(df["_resid_"], q=n_bins, duplicates="drop")
    df["_dist_bin_"] = pd.qcut(df["_dist_"], q=n_bins, duplicates="drop")

    if (
        df["_resid_bin_"].cat.categories.size != 3
        or df["_dist_bin_"].cat.categories.size != 3
    ):
        raise ValueError(
            "qcut produced fewer than 3 bins (ties). Consider jitter or fixed edges."
        )

    resid_cats = df["_resid_bin_"].cat.categories
    dist_cats = df["_dist_bin_"].cat.categories
    df["_resid_idx_"] = df["_resid_bin_"].cat.codes
    df["_dist_idx_"] = df["_dist_bin_"].cat.codes

    if colors_rgb is None:
        colors_rgb = _bivariate_palette_rowwise(
            n_resid=3, n_dist=3, dist_lo=palette_dist_lo, dist_hi=palette_dist_hi
        )

    colors_rgb = np.asarray(colors_rgb, dtype=float)
    if colors_rgb.shape != (3, 3, 3):
        raise ValueError("colors_rgb must have shape (3, 3, 3).")

    flat_hex = [
        to_hex(colors_rgb[r, d], keep_alpha=False) for r in range(3) for d in range(3)
    ]

    counts = (
        df.groupby(["_resid_idx_", "_dist_idx_"], observed=True)
        .size()
        .reindex(pd.MultiIndex.from_product([range(3), range(3)]), fill_value=0)
        .to_numpy()
        .reshape(3, 3)
    )

    if layout == "vertical":
        fig = make_subplots(
            rows=2,
            cols=1,
            row_heights=[0.38, 0.62],
            vertical_spacing=0.2,
            subplot_titles=(
                "Bivariate legend (Residual × Distance)",
                title,
            ),
        )
        legend_pos = dict(row=1, col=1)
        spatial_pos = dict(row=2, col=1)
        fig_width, fig_height = 700, 820
    else:
        fig = make_subplots(
            rows=1,
            cols=2,
            column_widths=[0.33, 0.67],
            horizontal_spacing=0.06,
            subplot_titles=(
                "Bivariate legend (Residual × Distance)",
                title,
            ),
        )
        legend_pos = dict(row=1, col=1)
        spatial_pos = dict(row=1, col=2)
        fig_width, fig_height = 900, 520

    z = np.arange(9).reshape(3, 3)

    eps = 1e-6
    cs: list[list[float | str]] = []
    for i, h in enumerate(flat_hex):
        t = i / 8.0
        cs.append([max(0.0, t - eps), h])
        cs.append([min(1.0, t + eps), h])

    fig.add_trace(
        go.Heatmap(z=z, colorscale=cs, showscale=False, hoverinfo="skip"),
        **legend_pos,
    )

    dist_labels = [
        f"{dist_cats[j].left:.0f}-{dist_cats[j].right:.0f} µm" for j in range(3)
    ]
    resid_labels_topdown = [str(c) for c in resid_cats][::-1]

    fig.update_xaxes(
        **legend_pos,
        tickmode="array",
        tickvals=[0, 1, 2],
        ticktext=dist_labels,
        tickangle=45,
        title_text="Distance bin",
        showgrid=False,
        zeroline=False,
    )
    fig.update_yaxes(
        **legend_pos,
        tickmode="array",
        tickvals=[0, 1, 2],
        ticktext=resid_labels_topdown,
        title_text="Residual bin",
        showgrid=False,
        zeroline=False,
        autorange="reversed",
    )

    if legend_show_counts:
        ann = list(fig.layout.annotations)
        for r in range(3):
            for d in range(3):
                n = int(counts[r, d])
                if (legend_count_min is not None) and (n < legend_count_min):
                    continue
                ann.append(
                    dict(
                        x=d,
                        y=2 - r,
                        xref="x1",
                        yref="y1",
                        text=legend_count_fmt.format(n),
                        showarrow=False,
                        font=dict(size=legend_count_fontsize, color="white"),
                    )
                )
        fig.update_layout(annotations=ann)

    ridx = df["_resid_idx_"].to_numpy()
    didx = df["_dist_idx_"].to_numpy()
    color_idx = (ridx * 3 + didx).astype(int)
    marker_colors = [flat_hex[i] for i in color_idx]

    fig.add_trace(
        go.Scattergl(
            x=df["_x_"],
            y=df["_y_"],
            mode="markers",
            marker=dict(size=point_size, color=marker_colors, opacity=point_alpha),
            showlegend=False,
            customdata=np.c_[df["_resid_"], df["_dist_"]],
            hovertemplate=(
                "x: %{x:.1f} µm<br>"
                "y: %{y:.1f} µm<br>"
                "Residual: %{customdata[0]:.2f}<br>"
                "Nearest-plaque dist: %{customdata[1]:.1f} µm<extra></extra>"
            ),
        ),
        **spatial_pos,
    )

    if (
        show_plaques
        and (plaque_x_col in cells_df.columns)
        and (plaque_y_col in cells_df.columns)
    ):
        plaques = cells_df[[plaque_x_col, plaque_y_col]].copy()
        plaques["_px_"] = pd.to_numeric(plaques[plaque_x_col], errors="coerce")
        plaques["_py_"] = pd.to_numeric(plaques[plaque_y_col], errors="coerce")
        plaques = plaques.loc[
            np.isfinite(plaques["_px_"]) & np.isfinite(plaques["_py_"])
        ]

        if plaque_round_decimals is not None:
            plaques["_px_"] = plaques["_px_"].round(plaque_round_decimals)
            plaques["_py_"] = plaques["_py_"].round(plaque_round_decimals)

        plaques = plaques.drop_duplicates(subset=["_px_", "_py_"])

        fig.add_trace(
            go.Scatter(
                x=plaques["_px_"],
                y=plaques["_py_"],
                mode="markers",
                marker=dict(
                    symbol="star",
                    size=plaque_size,
                    color=plaque_color,
                    opacity=plaque_alpha,
                    line=dict(color=plaque_edgecolor, width=plaque_linewidth),
                ),
                name=plaque_label,
                showlegend=bool(show_marker_legend),
                hoverinfo="skip",
            ),
            **spatial_pos,
        )

    fig.update_xaxes(**spatial_pos, title_text="X (µm)", showgrid=False, zeroline=False)
    fig.update_yaxes(
        **spatial_pos,
        title_text="Y (µm)",
        showgrid=False,
        zeroline=False,
        autorange="reversed",
        scaleanchor="x2",
        scaleratio=1,
    )

    fig.update_layout(
        title=f"{model_name}: Residual bin × distance bin (bivariate) - spatial map",
        width=fig_width,
        height=fig_height,
        template="simple_white",
        margin=dict(l=60, r=30, t=90, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            font=dict(size=marker_legend_fontsize),
            bgcolor=f"rgba(255,255,255,{marker_legend_framealpha})",
        ),
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
    fig.show()
    return fig


def plot_top_gene_importances_interactive(
    importance_df: pd.DataFrame,
    top_n: int = 20,
    *,
    title: str = "Top predictive genes across models (normalized importance)",
    filename: str = "top_gene_importance.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    """Interactive heatmap of top predictive genes across models (Plotly)."""

    normed = importance_df.groupby("model", group_keys=False).apply(
        lambda d: d.assign(norm_importance=d["importance"] / d["importance"].max())
    )

    top_genes = normed.groupby("model", group_keys=False).apply(
        lambda d: d.nlargest(top_n, "norm_importance")
    )

    pivot = top_genes.pivot_table(
        index="gene", columns="model", values="norm_importance", fill_value=0
    )

    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

    fig = go.Figure(
        data=go.Heatmap(
            z=pivot.to_numpy(dtype=float),
            x=pivot.columns.astype(str),
            y=pivot.index.astype(str),
            colorscale="Mako",
            colorbar=dict(title="Normalized Importance"),
            hovertemplate="Gene: %{y}<br>Model: %{x}<br>Norm importance: %{z:.3f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=title,
        autosize=True,
        template="simple_white",
        margin=dict(l=140, r=30, t=70, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
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

    if "geometry" not in df.columns:
        raise ValueError("df must contain a 'geometry' column with shapely objects.")
    if "is_convex" not in df.columns:
        raise ValueError("df must contain an 'is_convex' boolean column.")

    P = df.copy()

    if "plaque_id" not in P.columns:
        P["plaque_id"] = np.arange(1, len(P) + 1, dtype=int)
    if "area" not in P.columns:
        P["area"] = P["geometry"].map(lambda g: getattr(g, "area", np.nan))

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

    if brain_geom is not None and brain_geom.is_valid:
        g = rot(brain_geom)
        x, y = map(list, g.exterior.xy)
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

            x, y = map(list, poly.exterior.xy)

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
                    hovertemplate=(
                        f"Plaque {pid}<br>" f"Area: {area:,.0f} µm²<extra></extra>"
                    ),
                    showlegend=False,
                )
            )

    rng = np.random.default_rng(seed)
    n = min(len(P), int(sample_hulls))
    if n > 0:
        sample = P.sample(n=n, random_state=seed)
        for row in sample.itertuples(index=False):
            g = getattr(row, "geometry", None)
            if g is None:
                continue

            hull = rot(g.convex_hull)
            x, y = map(list, hull.exterior.xy)

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

    fig.update_layout(
        title=title or "Detected Aβ plaques after normalization",
        autosize=True,
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

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_path, include_plotlyjs="cdn")

    return fig


def make_cell_to_plaque_distance_map_plotly(
    *,
    cells_df: pd.DataFrame,
    plaques_gdf: pd.DataFrame | None = None,
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
    dist_col: str = "distance_to_plaque",
    clip_quantiles: tuple[float, float] = (0.01, 0.99),
    vmin: float | None = None,
    vmax: float | None = None,
    max_points: int | None = None,
    point_size: float = 4,
    point_alpha: float = 0.85,
    plaque_edgecolor: str = "cyan",
    plaque_linewidth: float = 1.2,
    invert_y: bool = False,
    title: str = "Cell–plaque distance map (µm)",
    figsize_px: tuple[int, int] = (820, 700),
    filename: str = "cell_to_plaque_distance_map.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    """
    Interactive spatial map of distance from each cell to nearest plaque.
    """

    C = cells_df.copy()

    if max_points is not None and len(C) > max_points:
        C = C.sample(n=max_points, random_state=0)

    dvals = C[dist_col].astype(float)

    lo = np.quantile(dvals, clip_quantiles[0]) if vmin is None else vmin
    hi = np.quantile(dvals, clip_quantiles[1]) if vmax is None else vmax

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=C[x_col],
            y=C[y_col],
            mode="markers",
            marker=dict(
                size=point_size,
                color=C[dist_col],
                colorscale="Plasma",
                cmin=lo,
                cmax=hi,
                opacity=point_alpha,
                colorbar=dict(
                    title="Distance to plaque (µm)",
                    ticks="outside",
                ),
            ),
            hovertemplate=(
                "Cell<br>" f"Distance: %{{marker.color:.1f}} µm<extra></extra>"
            ),
            showlegend=False,
        )
    )

    def iter_polygons(g):
        if isinstance(g, Polygon):
            yield g
        elif isinstance(g, MultiPolygon):
            for sub in g.geoms:
                if isinstance(sub, Polygon):
                    yield sub

    if plaques_gdf is not None and len(plaques_gdf):
        for row in plaques_gdf.itertuples(index=False):
            geom = getattr(row, "geometry", None)
            if geom is None:
                continue

            for poly in iter_polygons(geom):
                if poly.is_empty or not poly.is_valid:
                    continue

                x, y = map(list, poly.exterior.xy)

                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        line=dict(
                            color=plaque_edgecolor,
                            width=plaque_linewidth,
                        ),
                        hoverinfo="skip",
                        showlegend=False,
                    )
                )

    fig.update_layout(
        title=title,
        autosize=True,
        template="simple_white",
        margin=dict(l=60, r=40, t=60, b=50),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.update_xaxes(
        title="X (µm)",
        showgrid=False,
        zeroline=False,
        scaleanchor="y",
        scaleratio=1,
    )

    fig.update_yaxes(
        title="Y (µm)",
        showgrid=False,
        zeroline=False,
        autorange="reversed" if invert_y else True,
    )

    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_path, include_plotlyjs="cdn")

    return fig


def _infer_grid_shape(key_to_pos: dict[str, tuple[int, int]]) -> tuple[int, int]:
    rs = [r for r, _ in key_to_pos.values()]
    cs = [c for _, c in key_to_pos.values()]
    return (max(rs) + 1, max(cs) + 1)


def _load_png_rgb(path: Path) -> np.ndarray:
    img = Image.open(path).convert("RGB")
    return np.array(img)


def _pad_to_max(img: np.ndarray, *, target_h: int, target_w: int) -> np.ndarray:
    h, w = img.shape[:2]
    if h == target_h and w == target_w:
        return img
    out = np.zeros((target_h, target_w, 3), dtype=img.dtype)

    out[:h, :w, :] = img
    return out


def make_image_grid_interactive(
    *,
    files: dict[str, Path],
    key_to_pos: dict[str, tuple[int, int]],
    col_ticks: Iterable[str],
    row_ticks: Iterable[str],
    flip_h_keys: Iterable[str] = (),
    title: str | None = None,
    filename: str = "image_grid.html",
    out_dir: str = "frontend/public/plots",
) -> go.Figure:
    """
    Interactive image grid (Plotly) with unified row/col labels.

    Notes
    -----
    - Each image is horizontally flipped if its key is in flip_h_keys.
    - Each image is vertically cropped to remove the bottom 10% strip.
    - Images are padded to a common (max_h, max_w) before plotting.
    - Transparent background, responsive HTML export.
    """
    flip_h_keys = set(flip_h_keys)

    imgs: dict[str, np.ndarray] = {}
    shapes: list[tuple[int, int]] = []

    for key, path in files.items():
        img = _load_png_rgb(Path(path))
        if key in flip_h_keys:
            img = np.fliplr(img)

        h = img.shape[0]
        keep_h = max(1, int(round(h * 0.9)))
        img = img[:keep_h, ...]

        imgs[key] = img
        shapes.append(img.shape[:2])

    if not shapes:
        raise ValueError("No images loaded. `files` is empty?")

    max_h = max(h for h, _ in shapes)
    max_w = max(w for _, w in shapes)

    n_rows, n_cols = _infer_grid_shape(key_to_pos)

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        horizontal_spacing=0.02,
        vertical_spacing=0.02,
    )

    for key, (r0, c0) in key_to_pos.items():
        if key not in imgs:
            continue
        img = _pad_to_max(imgs[key], target_h=max_h, target_w=max_w)

        r = r0 + 1
        c = c0 + 1

        fig.add_trace(go.Image(z=img), row=r, col=c)

        fig.update_xaxes(visible=False, row=r, col=c)
        fig.update_yaxes(
            visible=False,
            row=r,
            col=c,
            scaleanchor=f"x{(r-1)*n_cols + c}",
            scaleratio=1,
        )

    col_ticks = list(col_ticks)
    row_ticks = list(row_ticks)

    for j in range(n_cols):
        x_center = (j + 0.5) / n_cols
        fig.add_annotation(
            x=x_center,
            y=1.02,
            xref="paper",
            yref="paper",
            text=str(col_ticks[j]) if j < len(col_ticks) else "",
            showarrow=False,
            font=dict(size=12),
        )

    for i in range(n_rows):
        y_center = 1 - (i + 0.5) / n_rows
        fig.add_annotation(
            x=-0.02,
            y=y_center,
            xref="paper",
            yref="paper",
            text=str(row_ticks[i]) if i < len(row_ticks) else "",
            showarrow=False,
            xanchor="right",
            font=dict(size=12),
        )

    fig.add_annotation(
        x=0.5,
        y=-0.06,
        xref="paper",
        yref="paper",
        text="Age (months)",
        showarrow=False,
        font=dict(size=13),
    )
    fig.add_annotation(
        x=-0.09,
        y=0.5,
        xref="paper",
        yref="paper",
        text="Type",
        showarrow=False,
        textangle=-90,
        font=dict(size=13),
    )

    fig.update_layout(
        title=title,
        autosize=True,
        margin=dict(l=90, r=20, t=80 if title else 55, b=80),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
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


def plot_leiden_logit_slopes_interactive(
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
    fig_width_min: float = 900,
    fig_width_per_bar: float = 35,
    fig_height: float = 600,
    ypad_scale: float = 0.001,
    rotation: int = 90,
    title: str | None = None,
    out_dir: str = "frontend/public/plots",
    filename: str = "slopes_types.html",
) -> tuple[go.Figure, pd.DataFrame]:
    """
    Interactive split bar plot of significant Leiden logistic slopes.
    Returns (fig, sig_df).
    """

    required = {pval_col, slope_col, cluster_col}
    missing = required - set(logit_df.columns)
    if missing:
        raise KeyError(f"logit_df is missing required columns: {sorted(missing)}")

    df = logit_df.sort_values(pval_col, ascending=True).reset_index(drop=True)
    sig = df.loc[df[pval_col] < pval_threshold].copy()

    def _map_cluster(x: Any) -> str:
        try:
            return my_label_to_type[int(x)]
        except Exception:
            return str(x)

    sig["cell_type"] = sig[cluster_col].apply(_map_cluster)

    neg = sig.loc[sig[slope_col] < 0].sort_values(slope_col, ascending=True)
    pos = sig.loc[sig[slope_col] > 0].sort_values(slope_col, ascending=True)

    nL, nR = len(neg), len(pos)
    xL = np.arange(nL)
    xR = np.arange(nR) + nL + gap

    fig_width = max(fig_width_min, fig_width_per_bar * (nL + nR + gap))

    fig = go.Figure()

    fig.add_bar(
        x=xL,
        y=neg[slope_col],
        marker_color=blue,
        name="Negative slope",
        hovertemplate="Cell type: %{customdata}<br>Slope: %{y:.4f}<extra></extra>",
        customdata=neg["cell_type"],
    )

    fig.add_bar(
        x=xR,
        y=pos[slope_col],
        marker_color=red,
        name="Positive slope",
        hovertemplate="Cell type: %{customdata}<br>Slope: %{y:.4f}<extra></extra>",
        customdata=pos["cell_type"],
    )

    fig.add_hline(y=0, line_width=1, line_color="black")

    max_abs = float(sig[slope_col].abs().max()) if len(sig) else 1.0
    ypad = ypad_scale * max_abs if max_abs > 0 else 0.1

    annotations = []

    for x, y, txt in zip(xL, neg[slope_col], neg["cell_type"], strict=False):
        annotations.append(
            dict(
                x=x,
                y=y - ypad,
                text=str(txt),
                showarrow=False,
                textangle=rotation,
                xanchor="center",
                yanchor="top",
                font=dict(size=11),
            )
        )

    for x, y, txt in zip(xR, pos[slope_col], pos["cell_type"], strict=False):
        annotations.append(
            dict(
                x=x,
                y=y + ypad,
                text=str(txt),
                showarrow=False,
                textangle=rotation,
                xanchor="center",
                yanchor="bottom",
                font=dict(size=11),
            )
        )

    fig.update_layout(
        title=title,
        annotations=annotations,
        xaxis=dict(showticklabels=False),
        yaxis=dict(title=slope_col),
        autosize=False,
        width=fig_width,
        height=fig_height,
        margin=dict(l=80, r=30, t=80 if title else 50, b=40),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
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

    return fig, sig


def make_cluster_frequency_distance_to_plaque_plotly(
    freq_df: pd.DataFrame,
    *,
    expanded_types: Mapping | None = None,
    cluster_col: str = "cluster_leiden",
    x_col: str = "bin_mid",
    y_col: str = "pct",
    title: str = "Cluster frequency (%) vs. distance to nearest plaque",
    filename: str = "cluster_frequency_distance_to_plaque.html",
    out_dir: str = "frontend/public/plots",
) -> "px.Figure":
    """
    Interactive version of the seaborn lineplot:
      x = distance bin midpoint
      y = % of cells in each bin
      hue = cluster name (expanded_types mapping if provided)

    Expects freq_df to already contain columns: [cluster_col, x_col, y_col].
    """

    dfp = freq_df.copy()

    if expanded_types is not None:
        dfp["cluster_name"] = dfp[cluster_col].map(expanded_types)

        dfp["cluster_name"] = dfp["cluster_name"].astype(object)
        miss = dfp["cluster_name"].isna()
        if miss.any():
            dfp.loc[miss, "cluster_name"] = dfp.loc[miss, cluster_col].astype(str)
    else:
        dfp["cluster_name"] = dfp[cluster_col].astype(str)

    dfp[x_col] = pd.to_numeric(dfp[x_col], errors="coerce")
    dfp[y_col] = pd.to_numeric(dfp[y_col], errors="coerce")
    dfp = dfp.dropna(subset=[x_col, y_col, "cluster_name"])

    dfp = dfp.sort_values(["cluster_name", x_col])

    fig = px.line(
        dfp,
        x=x_col,
        y=y_col,
        color="cluster_name",
        markers=True,
        title=title,
        labels={
            x_col: "Distance to plaque",
            y_col: "% of cells in each bin",
            "cluster_name": "Cluster",
        },
    )

    fig.update_traces(opacity=0.7)

    fig.update_layout(
        autosize=True,
        template="simple_white",
        margin=dict(l=70, r=30, t=80, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend_title_text="Cluster",
        legend=dict(
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
        ),
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

def plot_resid_vs_distance_combined_plotly(models, dist_col, title):
    """
    Interactive Plotly version of plot_resid_vs_distance_combined.

    Writes:
      frontend/public/plots/residual_v_dist_true_perm_fake.html
    """

    out_html = Path("frontend/public/plots/residual_v_dist_true_perm_fake.html")
    out_html.parent.mkdir(parents=True, exist_ok=True)

    fig = go.Figure()

    # Matplotlib default color cycle
    mpl_colors = plt.rcParams["axes.prop_cycle"].by_key().get(
        "color", ["C0", "C1", "C2", "C3"]
    )

    for i, m in enumerate(models):
        name = m["name"]
        dfi = m["df"]

        if dist_col not in dfi.columns:
            raise ValueError(f"Model '{name}' df missing '{dist_col}'. Have: {list(dfi.columns)}")
        if "oof_resid" not in dfi.columns:
            raise ValueError(f"Model '{name}' df missing 'oof_resid'. Have: {list(dfi.columns)}")

        x = dfi[dist_col].to_numpy(dtype=float)
        y = dfi["oof_resid"].to_numpy(dtype=float)

        color = mpl_colors[i % len(mpl_colors)]

        # scatter
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                mode="markers",
                name=name,
                marker=dict(size=6, color=color, opacity=0.6),
                hovertemplate=(
                    f"{dist_col}=%{{x:.3f}}<br>"
                    "Residual=%{y:.3f}"
                    "<extra></extra>"
                ),
            )
        )

        # linear fit (like _plot_scatter_with_fit)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() >= 2:
            coef = np.polyfit(x[mask], y[mask], 1)
            xfit = np.array([x[mask].min(), x[mask].max()])
            yfit = coef[0] * xfit + coef[1]

            fig.add_trace(
                go.Scatter(
                    x=xfit,
                    y=yfit,
                    mode="lines",
                    line=dict(color=color, width=2),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    # layout: transparent background + matplotlib-like axes
    fig.update_layout(
        title=title,
        xaxis_title=dist_col,
        yaxis_title="OOF Residual (mean obs - mean pred)",
        legend=dict(bgcolor="rgba(0,0,0,0)", borderwidth=0),
        margin=dict(l=70, r=40, t=70, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width = 900,
        autosize = True
    )

    fig.update_xaxes(showline=True, mirror=False, linecolor="black", zeroline=False)
    fig.update_yaxes(showline=True, mirror=False, linecolor="black", zeroline=False)

    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config=dict(responsive=True, displayModeBar=True, displaylogo=False),
        auto_open=False,
    )

    return fig

def plot_pred_vs_obs_combined_plotly(models, target_gene=None, title=""):
    out_html = Path("frontend/public/plots/true_permitted_fake_pred_vs_obs.html")
    out_html.parent.mkdir(parents=True, exist_ok=True)

    fig = go.Figure()

    mpl_colors = plt.rcParams["axes.prop_cycle"].by_key().get(
        "color", ["C0", "C1", "C2", "C3"]
    )

    for i, m in enumerate(models):
        name = m["name"]
        dfi = m["df"]
        r2 = m.get("r2", None)

        if "oof_pred" not in dfi.columns:
            raise ValueError(f"Model '{name}' df is missing required column 'oof_pred'. Have: {list(dfi.columns)}")

        y_col = target_gene if (target_gene is not None and target_gene in dfi.columns) else "y"
        if y_col not in dfi.columns:
            raise ValueError(f"Model '{name}' df is missing y column '{y_col}'. Have: {list(dfi.columns)}")

        obs = dfi[y_col].to_numpy(dtype=float)
        pred = dfi["oof_pred"].to_numpy(dtype=float)

        color = mpl_colors[i % len(mpl_colors)]
        label = name if r2 is None else f"{name} (mean R²={r2:.3f})"

        # scatter
        fig.add_trace(
            go.Scattergl(
                x=obs,
                y=pred,
                mode="markers",
                name=label,
                marker=dict(size=6, color=color, opacity=0.6),
                hovertemplate="Observed=%{x:.3f}<br>Predicted=%{y:.3f}<extra></extra>",
            )
        )

        # fit line
        mask = np.isfinite(obs) & np.isfinite(pred)
        if mask.sum() >= 2:
            coef = np.polyfit(obs[mask], pred[mask], 1)
            xfit = np.array([obs[mask].min(), obs[mask].max()])
            yfit = coef[0] * xfit + coef[1]

            fig.add_trace(
                go.Scatter(
                    x=xfit,
                    y=yfit,
                    mode="lines",
                    line=dict(color=color, width=2),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    fig.update_layout(
        title=title,
        xaxis_title="Observed (mean across genes)",
        yaxis_title="OOF Predicted (mean across genes)",
        # Plotly-compatible legend styling (matplotlib-like, no box)
        legend=dict(
            bgcolor="rgba(0,0,0,0)",
            bordercolor="rgba(0,0,0,0)",
            borderwidth=0,
        ),
        margin=dict(l=70, r=40, t=70, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=900
    )

    fig.update_xaxes(showline=True, mirror=True, linecolor="black", zeroline=False)
    fig.update_yaxes(showline=True, mirror=True, linecolor="black", zeroline=False)

    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config=dict(responsive=True, displayModeBar=True, displaylogo=False),
        auto_open=False,
    )

    return fig


def plot_tile_heatmap_plotly(df_in, x_col, y_col, val_col, title, n_tiles_x=10, n_tiles_y=10):
    """
    Interactive Plotly version of plot_tile_heatmap.
    Same call signature as original.
    Writes: frontend/public/plots/perm_vs_true_tiles.html
    """

    # --- compute tile ids (expects your existing tile_ids helper) ---
    tx, ty = tile_ids(df_in, x_col, y_col, n_tiles_x, n_tiles_y)

    tmp = df_in.copy()
    tmp["tx"] = tx
    tmp["ty"] = ty

    # mean per tile
    mat = tmp.groupby(["ty", "tx"])[val_col].mean().unstack("tx")

    # Ensure full tile grid exists (fill missing tiles with NaN for consistent shape)
    mat = mat.reindex(index=range(n_tiles_y), columns=range(n_tiles_x))

    Z = mat.to_numpy(dtype=float)

    out_html = Path("frontend/public/plots/perm_vs_true_tiles.html")
    out_html.parent.mkdir(parents=True, exist_ok=True)

    fig = go.Figure(
        data=go.Heatmap(
            z=Z,
            colorscale="Viridis",  # close to matplotlib imshow default look
            colorbar=dict(title=f"mean({val_col})"),
            hovertemplate=(
                "tile x=%{x}<br>"
                "tile y=%{y}<br>"
                f"mean({val_col})=%{{z:.4g}}"
                "<extra></extra>"
            ),
        )
    )

    # Match matplotlib-ish axes/labels
    fig.update_layout(
        title=title,
        xaxis_title="tile x",
        yaxis_title="tile y",
        margin=dict(l=70, r=40, t=70, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=900,
        autosize=True
    )

    # Make it look like imshow (cells are pixels, y increasing downward in image)
    # If you prefer origin="lower" behavior, remove autorange="reversed".
    fig.update_yaxes(autorange="reversed")

    # Visible axis lines like matplotlib
    fig.update_xaxes(showline=True, mirror=False, linecolor="black", zeroline=False)
    fig.update_yaxes(showline=True, mirror=False, linecolor="black", zeroline=False)

    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config=dict(responsive=True, displayModeBar=True, displaylogo=False),
        auto_open=False,
    )

    return fig

def plot_spatial_scatter_plotly(
    df_base,
    x_col,
    y_col,
    val_col,
    title,
    s=6,
    *,
    out_html: str | Path = "frontend/public/plots/permuted_neighbors_spatial.html",
    show: bool = False,
):
    """
    Interactive Plotly version of plot_spatial_scatter.

    Preserves:
      - scatter coloring by val_col
      - colorbar
      - transparent background
      - matplotlib-like default style
    """

    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x=df_base[x_col],
            y=df_base[y_col],
            mode="markers",
            marker=dict(
                size=s,
                color=df_base[val_col],
                colorscale="Viridis",  # matplotlib default-like
                showscale=True,
                colorbar=dict(title=val_col),
            ),
            hovertemplate=(
                f"{x_col}=%{{x:.3f}}<br>"
                f"{y_col}=%{{y:.3f}}<br>"
                f"{val_col}=%{{marker.color:.3f}}"
                "<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title="Permuted minus true closest",
        yaxis_title=y_col,
        margin=dict(l=70, r=40, t=70, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width= 900,
        autosize=True
    )

    fig.update_xaxes(
        showline=True,
        mirror=False,
        linecolor="black",
        zeroline=False,
    )
    fig.update_yaxes(
        showline=True,
        mirror=False,
        linecolor="black",
        zeroline=False,
        scaleanchor="x",   # preserves aspect if coordinates are spatial
        scaleratio=1,
    )

    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)

    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config=dict(responsive=True, displayModeBar=True, displaylogo=False),
        auto_open=False,
    )

    if show:
        fig.show()

    return fig


def plot_true_vs_fake_far_neighbors_across_genes(
    df,
    *,
    target_genes,
    groups,
    x_col,
    y_col,
    dist_col,
    optional_covs,
    neigh_cols,
    k_far=100,
    n_splits=5,
    alpha=1.0,
    confidence=0.95,
    maroon="maroon",
    capsize=6,
    figsize=(9, 6),
    pct_eps=1e-12,
    # ---- new ----
    out_html: str | Path = "frontend/public/plots/fake_neighbors.html",
    show: bool = False,
):
    """
    Plotly interactive version of plot_true_vs_fake_far_neighbors_across_genes.

    Produces:
      - res_per_gene
      - res_summary
      - interactive bar plot with CI + min/max gene annotations
    """

    # ------------------------------------------------------------------
    # Build fake features ONCE
    # ------------------------------------------------------------------
    df_fake, fake_all_cols = build_fake_far_neighbor_means(
        df, x_col, y_col, neigh_cols, k=k_far
    )

    true2fake = {
        c: f"fake_far_neigh_mean_{c.replace('neigh_mean_', '')}"
        for c in neigh_cols
    }

    missing_fake = [true2fake[c] for c in neigh_cols if true2fake[c] not in df_fake.columns]
    if missing_fake:
        raise RuntimeError(f"Missing fake columns in df_fake: {missing_fake[:10]}")

    # ------------------------------------------------------------------
    # Per-gene OOF R²
    # ------------------------------------------------------------------
    rows = []
    for gene in target_genes:
        current_true_neighbors = [c for c in neigh_cols if gene not in c]
        current_fake_neighbors = [true2fake[c] for c in current_true_neighbors]

        feats_true = [dist_col] + list(optional_covs) + current_true_neighbors
        feats_fake = [dist_col] + list(optional_covs) + current_fake_neighbors

        r2_true, _ = _oof_r2_spatial_blocks_ridge(
            df, gene, feats_true, groups,
            n_splits=n_splits, alpha=alpha
        )
        r2_fake, _ = _oof_r2_spatial_blocks_ridge(
            df_fake, gene, feats_fake, groups,
            n_splits=n_splits, alpha=alpha
        )

        rows.append(
            dict(
                gene=gene,
                R2_true_neighbors=float(r2_true),
                R2_fake_far_neighbors=float(r2_fake),
                gap_true_minus_fake=100.0 * (r2_true - r2_fake) / (abs(r2_true) + pct_eps),
            )
        )

    res_per_gene = pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Summary (mean ± CI)
    # ------------------------------------------------------------------
    m_true, lo_true, hi_true = mean_ci_t(
        res_per_gene["R2_true_neighbors"].to_numpy(), confidence=confidence
    )
    m_fake, lo_fake, hi_fake = mean_ci_t(
        res_per_gene["R2_fake_far_neighbors"].to_numpy(), confidence=confidence
    )

    avg_gap = float(res_per_gene["gap_true_minus_fake"].mean())

    res_summary = pd.DataFrame(
        [
            dict(condition="true_neighbors", mean_R2=m_true,
                 ci95_low=lo_true, ci95_high=hi_true, n_genes=len(res_per_gene)),
            dict(condition="fake_far_neighbors", mean_R2=m_fake,
                 ci95_low=lo_fake, ci95_high=hi_fake, n_genes=len(res_per_gene)),
        ]
    )

    # ------------------------------------------------------------------
    # Plotly plotting
    # ------------------------------------------------------------------
    BLUE = "#1f77b4"    # matplotlib tab:blue
    ORANGE = "#ff7f0e"  # matplotlib tab:orange

    labels = [
        "100 closest neighbors",
        f"100 farthest neighbors (avg drop = {avg_gap:.1f})",
    ]

    means = np.array([m_true, m_fake])
    err_plus = np.array([hi_true - m_true, hi_fake - m_fake])
    err_minus = np.array([m_true - lo_true, m_fake - lo_fake])
    x = np.array([0.0, 1.0])

    fig = go.Figure()

    # Bars + CI
    fig.add_trace(
        go.Bar(
            x=x,
            y=means,
            width=0.65,
            marker=dict(color=[BLUE, BLUE]),
            showlegend=False,
            error_y=dict(
                type="data",
                symmetric=False,
                array=err_plus,
                arrayminus=err_minus,
                color="black",
                thickness=2,
                width=capsize,
                visible=True,
            ),
            hovertemplate="mean=%{y:.4f}<extra></extra>",
        )
    )

    # ------------------------------------------------------------------
    # Min / max gene points (maroon)
    # ------------------------------------------------------------------
    def _minmax(col, xi):
        sub = res_per_gene[["gene", col]]
        imin = sub[col].idxmin()
        imax = sub[col].idxmax()
        idxs = [imin, imax] if imin != imax else [imin]

        xs, ys, txt = [], [], []
        for i in idxs:
            g = sub.loc[i, "gene"]
            r2 = float(sub.loc[i, col])
            xs.append(xi)
            ys.append(r2)
            txt.append(f"{g} ({ceil_to_decimals(r2, 2):.2f})")
        return xs, ys, txt

    x0, y0, t0 = _minmax("R2_true_neighbors", 0.0)
    x1, y1, t1 = _minmax("R2_fake_far_neighbors", 1.0)

    fig.add_trace(
        go.Scatter(
            x=x0 + x1,
            y=y0 + y1,
            mode="markers+text",
            marker=dict(color=maroon, size=10),
            text=t0 + t1,
            textposition="top right",
            textfont=dict(color=maroon, size=9),
            hovertemplate="%{text}<extra></extra>",
            showlegend=False,
        )
    )

    # Axes & layout
    fig.update_xaxes(
        tickmode="array",
        tickvals=x,
        ticktext=labels,
        tickangle=20,
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
        range=[-0.6, 1.6],
    )

    fig.update_yaxes(
        title=f"OOF R² (spatial-block GroupKFold), mean across {len(res_per_gene)} genes",
        showgrid=True,
        gridcolor="rgba(0,0,0,0.25)",
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
    )

    fig.update_layout(
        title=f"Closest vs farthest-neighbor features (mean ± {int(confidence*100)}% CI)",
        margin=dict(l=80, r=30, t=80, b=160),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=int(figsize[0] * 90),
        height=int(figsize[1] * 90),
    )

    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config=dict(responsive=True, displayModeBar=True, displaylogo=False),
        auto_open=False,
    )

    if show:
        fig.show()

    return res_per_gene, res_summary


def plot_permutation_effect_across_genes(
    df,
    target_genes,
    groups,
    *,
    dist_col,
    neigh_cols,
    optional_covs=None,
    seed=0,
    n_splits=5,
    alpha=1.0,
    confidence=0.95,
    maroon="maroon",
    capsize=6,
    figsize=(9, 6),
    pct_eps=1e-12,
    # --- new: output ---
    out_html: str | Path = "frontend/public/plots/full_model_neighbor_permute.html",
    show: bool = False,
):
    """
    Drop-in replacement with Plotly output.

    Returns:
      - res_per_gene
      - res_summary
    And writes:
      - frontend/public/plots/full_model_neighbor_permute.html
    """

    if optional_covs is None:
        optional_covs = []

    # --- local helper: pooled OOF R^2 under GroupKFold ---
    def _oof_r2_spatial_blocks_ridge(df_in, target, feature_cols, groups_in, *, n_splits, alpha):
        from sklearn.model_selection import GroupKFold
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.linear_model import Ridge
        from sklearn.metrics import r2_score

        if len(groups_in) != len(df_in):
            raise ValueError("groups must be aligned with df rows (same length)")

        use = df_in[[target] + list(feature_cols)].copy()
        use["__group__"] = np.asarray(groups_in)

        use = use.dropna()
        if len(use) == 0:
            raise ValueError(f"No rows left after dropna for target={target} and features={feature_cols}")

        y = use[target].to_numpy(dtype=float)
        X = use[list(feature_cols)].to_numpy(dtype=float)
        g = use["__group__"].to_numpy()

        n_groups = len(np.unique(g))
        if n_groups < n_splits:
            raise ValueError(
                f"Not enough unique groups for GroupKFold: n_groups={n_groups} < n_splits={n_splits}"
            )

        gkf = GroupKFold(n_splits=n_splits)

        oof_pred = np.full(len(use), np.nan, dtype=float)
        fold_r2s = []

        for tr, te in gkf.split(X, y, groups=g):
            model = Pipeline([
                ("scaler", StandardScaler(with_mean=True, with_std=True)),
                ("ridge", Ridge(alpha=alpha, random_state=0)),
            ])
            model.fit(X[tr], y[tr])
            pred = model.predict(X[te])

            oof_pred[te] = pred
            fold_r2s.append(r2_score(y[te], pred))

        if not np.isfinite(oof_pred).all():
            raise RuntimeError("OOF predictions contain NaNs; check data, groups, and split coverage.")

        r2_oof = r2_score(y, oof_pred)
        return float(r2_oof), fold_r2s

    # permute once (same permuted dataset used for all genes)
    df_perm = permute_within_groups(df, neigh_cols, groups, seed=seed)

    rows = []
    for gene in target_genes:
        # exclude any neighbor columns that correspond to the current gene
        current_neighbors = [c for c in neigh_cols if gene not in c]
        full_feats = [dist_col] + list(optional_covs) + current_neighbors

        r2_true, _ = _oof_r2_spatial_blocks_ridge(
            df, gene, full_feats, groups, n_splits=n_splits, alpha=alpha
        )
        r2_perm, _ = _oof_r2_spatial_blocks_ridge(
            df_perm, gene, full_feats, groups, n_splits=n_splits, alpha=alpha
        )

        pct_drop = 100.0 * (r2_true - r2_perm) / (np.abs(r2_true) + pct_eps)

        rows.append({
            "gene": gene,
            "R2_true": float(r2_true),
            "R2_permuted": float(r2_perm),
            "drop_abs": float(r2_true - r2_perm),
            "drop_pct": float(pct_drop),
        })

    res_per_gene = pd.DataFrame(rows)

    # summaries + CI across genes
    m_true, lo_true, hi_true = mean_ci_t(res_per_gene["R2_true"].to_numpy(), confidence=confidence)
    m_perm, lo_perm, hi_perm = mean_ci_t(res_per_gene["R2_permuted"].to_numpy(), confidence=confidence)

    avg_pct_drop = float(np.mean(res_per_gene["drop_pct"].to_numpy()))

    res_summary = pd.DataFrame([
        {"condition": "true_neighbors",     "mean_R2": m_true, "ci95_low": lo_true, "ci95_high": hi_true, "n_genes": len(res_per_gene)},
        {"condition": "permuted_neighbors", "mean_R2": m_perm, "ci95_low": lo_perm, "ci95_high": hi_perm, "n_genes": len(res_per_gene)},
    ])

    # ---------- Plotly plotting ----------
    # Preserve matplotlib default bar colors: tab:blue then tab:orange
    BLUE = "#1f77b4"
    ORANGE = "#ff7f0e"

    labels = [
        "true neighbors",
        f"permuted within a tile (avg drop = {avg_pct_drop:.1f}%)",
    ]

    means = np.array([m_true, m_perm], dtype=float)
    ci_lows = np.array([lo_true, lo_perm], dtype=float)
    ci_highs = np.array([hi_true, hi_perm], dtype=float)

    err_plus = ci_highs - means
    err_minus = means - ci_lows

    x = np.array([0.0, 1.0], dtype=float)

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=x,
            y=means,
            width=0.65,
            marker=dict(color=[BLUE, BLUE]),
            showlegend=False,
            hovertemplate="mean=%{y:.4f}<extra></extra>",
            error_y=dict(
                type="data",
                symmetric=False,
                array=err_plus,
                arrayminus=err_minus,
                color="black",
                thickness=2,
                width=capsize,  # approximates matplotlib capsize in px
                visible=True,
            ),
        )
    )

    # maroon min/max per condition, with labels
    def _minmax_points(condition_col: str, xi: float):
        sub = res_per_gene[["gene", condition_col]].copy()
        imin = sub[condition_col].idxmin()
        imax = sub[condition_col].idxmax()
        idxs = [imin, imax] if imax != imin else [imin]

        pts_x, pts_y, pts_text = [], [], []
        for idx in idxs:
            gene = sub.loc[idx, "gene"]
            r2 = float(sub.loc[idx, condition_col])
            r2_disp = float(ceil_to_decimals(r2, 2))
            pts_x.append(xi)
            pts_y.append(r2)
            pts_text.append(f"{gene} ({r2_disp:.2f})")
        return pts_x, pts_y, pts_text

    x0, y0, t0 = _minmax_points("R2_true", 0.0)
    x1, y1, t1 = _minmax_points("R2_permuted", 1.0)

    sx = x0 + x1
    sy = y0 + y1
    st = t0 + t1

    if len(sx):
        fig.add_trace(
            go.Scatter(
                x=sx,
                y=sy,
                mode="markers+text",
                marker=dict(color=maroon, size=10),
                text=st,
                textposition="top right",
                textfont=dict(color=maroon, size=9),
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            )
        )

    fig.update_xaxes(
        tickmode="array",
        tickvals=x,
        ticktext=labels,
        tickangle=20,
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
        range=[-0.6, 1.6],
    )

    fig.update_yaxes(
        title=f"OOF R² (spatial-block GroupKFold), mean across {len(res_per_gene)} genes",
        showgrid=True,
        gridcolor="rgba(0,0,0,0.25)",
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
    )

    fig.update_layout(
        title=f"Permutation test within tiles (OOF R² mean ± {int(confidence*100)}% CI)",
        margin=dict(l=80, r=30, t=80, b=160),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=int(figsize[0] * 90),
        height=int(figsize[1] * 90),
    )

    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config={"responsive": True, "displayModeBar": True, "displaylogo": False},
        auto_open=False,
    )

    if show:
        fig.show()

    return res_per_gene, res_summary

def mean_ci_t(x: np.ndarray, confidence: float = 0.95) -> tuple[float, float, float]:
    """
    Mean and two-sided t-based CI for the mean.
    Returns (mean, lo, hi).
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return (np.nan, np.nan, np.nan)
    m = float(np.mean(x))
    if n == 1:
        return (m, m, m)
    sem = stats.sem(x)
    df = n - 1
    tcrit = stats.t.ppf((1 + confidence) / 2, df)
    half = float(tcrit * sem)
    return (m, m - half, m + half)


def ceil_to_decimals(x: float, decimals: int = 2) -> float:
    p = 10**decimals
    return float(np.ceil(x * p) / p)


def plot_model_comparison_with_extremes(
    res_per_gene: pd.DataFrame,
    confidence: float = 0.95,
    figsize=(10, 6),
    bar_alpha: float = 1.0,
    capsize: int = 6,
    point_color: str = "maroon",
    annotate: bool = True,
    annotate_fontsize: int = 9,
    text_dx: float = 0.02,
    text_dy: float = 0.002,
    # output (added, but optional; call style still matches original)
    out_html: str | Path = "frontend/public/plots/interaction_model_performance.html",
    show: bool = False,
):
    """
    Plotly interactive version.

    res_per_gene must contain columns:
      - gene
      - model
      - R2_spatialBlockCV

    Produces:
      - res_summary dataframe (mean across genes + t-based CI)
      - bar plot with black error bars
      - scatter points for min/max gene per model, annotated
      - transparent background
      - writes to interaction_model_performance.html
    """
    required = {"gene", "model", "R2_spatialBlockCV"}
    missing = required - set(res_per_gene.columns)
    if missing:
        raise ValueError(f"res_per_gene missing columns: {missing}")

    # ---- summary across genes with t-based CI ----
    summary_rows = []
    for model_name, sub in res_per_gene.groupby("model"):
        gene_vals = sub["R2_spatialBlockCV"].to_numpy()
        m, lo, hi = mean_ci_t(gene_vals, confidence=confidence)
        summary_rows.append(
            dict(
                model=model_name,
                mean_R2_across_genes=m,
                ci95_low=lo,
                ci95_high=hi,
                n_genes=len(gene_vals),
            )
        )

    res_summary = (
        pd.DataFrame(summary_rows)
        .sort_values("mean_R2_across_genes", ascending=False)
        .reset_index(drop=True)
    )

    # ---- min/max gene per model ----
    extremes = {}
    for model_name in res_summary["model"]:
        sub = res_per_gene[res_per_gene["model"] == model_name].copy()
        i_min = sub["R2_spatialBlockCV"].idxmin()
        i_max = sub["R2_spatialBlockCV"].idxmax()
        rows = [sub.loc[i_min]]
        if i_max != i_min:
            rows.append(sub.loc[i_max])
        extremes[model_name] = rows

    # ---- plot (numeric x like matplotlib) ----
    n = len(res_summary)
    x = np.arange(n, dtype=float)

    means = res_summary["mean_R2_across_genes"].to_numpy(dtype=float)
    lo = res_summary["ci95_low"].to_numpy(dtype=float)
    hi = res_summary["ci95_high"].to_numpy(dtype=float)

    err_plus = hi - means
    err_minus = means - lo

    # preserve matplotlib default blue
    BAR_COLOR = "#1f77b4"

    fig = go.Figure()

    # Bars
    fig.add_trace(
        go.Bar(
            x=x,
            y=means,
            marker=dict(color=BAR_COLOR),
            opacity=bar_alpha,
            width=0.7,
            showlegend=False,
            hovertemplate="mean=%{y:.4f}<extra></extra>",
        )
    )

    # Error bars (black, like mpl errorbar)
    # Plotly error bars support width/thickness; capsize in px is approximated via 'width'
    fig.update_traces(
        error_y=dict(
            type="data",
            symmetric=False,
            array=err_plus,
            arrayminus=err_minus,
            color="black",
            thickness=2,
            width=capsize,   # approximate capsize
            visible=True,
        ),
        selector=dict(type="bar"),
    )

    # Extremes points + annotations
    px, py, ptext = [], [], []
    for xi, model_name in enumerate(res_summary["model"]):
        for row in extremes[model_name]:
            gene = row["gene"]
            r2 = float(row["R2_spatialBlockCV"])
            r2_disp = float(ceil_to_decimals(r2, 2))

            px.append(float(xi))
            py.append(r2)
            ptext.append(f"{gene} ({r2_disp:.2f})")

    if px:
        fig.add_trace(
            go.Scatter(
                x=px,
                y=py,
                mode=("markers+text" if annotate else "markers"),
                marker=dict(color=point_color, size=10),
                text=ptext if annotate else None,
                textposition="top right",
                textfont=dict(size=annotate_fontsize, color=point_color),
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            )
        )

    # X tick labels (rotated)
    fig.update_xaxes(
        tickmode="array",
        tickvals=x,
        ticktext=res_summary["model"].tolist(),
        tickangle=25,
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
    )

    n_genes_hint = int(res_summary["n_genes"].iloc[0]) if len(res_summary) else 0

    fig.update_yaxes(
        title=f"OOF R² (spatial-block CV), mean across {n_genes_hint} genes",
        showgrid=True,
        gridcolor="rgba(0,0,0,0.25)",
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
    )

    fig.update_layout(
        title=f"Model comparison across genes (OOF R² mean ± {int(confidence*100)}% CI) + min/max gene points",
        margin=dict(l=70, r=30, t=80, b=140),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=int(figsize[0] * 90),
        height=int(figsize[1] * 90),
    )

    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config={"responsive": True, "displayModeBar": True, "displaylogo": False},
        auto_open=False,
    )

    if show:
        fig.show()

    return res_summary

def _apply_color_scale(A: np.ndarray, *, color_scale: str, eps: float, symlog_linthresh: float) -> np.ndarray:
    A = np.asarray(A, dtype=float)

    if color_scale == "linear":
        return A

    if color_scale == "log1p":
        return np.log1p(A)

    if color_scale == "log10":
        return np.log10(A + eps)

    if color_scale == "lognorm":
        # For display, LogNorm is effectively log scaling.
        # We'll use log10(A+eps) as a close visual analogue.
        return np.log10(A + eps)

    if color_scale == "symlognorm":
        # SymLog-like: linear near 0, log further out
        # Here A is non-negative (counts), so it's basically log-like with a linear region.
        # We mimic by: log10(1 + A/linthresh)
        return np.log10(1.0 + (A / symlog_linthresh))

    raise ValueError("color_scale must be one of: linear, log1p, log10, lognorm, symlognorm")


def plot_dist_signature_heatmaps_grid_plotly(
    df,
    dist_col,
    sig_col,
    y_cols,                 # list of 16 gene columns
    q=5,
    y_is_log1p=True,        # your current y: log1p(count)
    aggfunc="mean",
    nrows=4,
    ncols=4,
    figsize=(18, 14),
    cmap="viridis",
    # --- knobs ---
    color_scale="log1p",    # {"linear", "log1p", "log10", "lognorm", "symlognorm"}
    norm_scope="per_gene",  # {"per_gene", "global"}
    share_colorbar=False,   # if True, forces norm_scope="global"
    eps=1e-6,
    symlog_linthresh=0.1,
    # --- annotation ---
    annotate=True,
    annot_decimals=2,
    annot_fontsize=7,
    label_mode="index",     # {"index", "interval"}
    # --- output ---
    out_html: str | Path = "frontend/public/plots/interaction_model.html",
    show: bool = False,
):
    """
    Interactive Plotly version of the matplotlib function.
    Writes a single HTML with a 4x4 grid of heatmaps.
    Background is transparent.
    """

    if len(y_cols) != nrows * ncols:
        raise ValueError(f"Expected {nrows*ncols} genes, got {len(y_cols)}")

    if share_colorbar:
        norm_scope = "global"

    # ---- 1) Compute common bin edges once ----
    base = df[[dist_col, sig_col]].dropna()
    _, dist_edges = pd.qcut(base[dist_col], q=q, duplicates="drop", retbins=True)
    _, sig_edges = pd.qcut(base[sig_col], q=q, duplicates="drop", retbins=True)

    dist_cats = pd.cut(base[dist_col], bins=dist_edges, include_lowest=True).cat.categories
    sig_cats = pd.cut(base[sig_col], bins=sig_edges, include_lowest=True).cat.categories

    if label_mode == "interval":
        xlabels = [str(c) for c in dist_cats]
        ylabels = [str(c) for c in sig_cats]
    elif label_mode == "index":
        xlabels = [f"D{i+1}" for i in range(len(dist_cats))]
        ylabels = [f"S{i+1}" for i in range(len(sig_cats))]
    else:
        raise ValueError("label_mode must be 'index' or 'interval'")

    # ---- 2) Build mean tables (natural scale) + colored arrays ----
    mean_tables: list[pd.DataFrame] = []
    colored_arrays: list[np.ndarray] = []

    for gene in y_cols:
        tmp = df[[dist_col, sig_col, gene]].dropna().copy()
        tmp["dist_bin"] = pd.cut(tmp[dist_col], bins=dist_edges, include_lowest=True)
        tmp["sig_bin"] = pd.cut(tmp[sig_col], bins=sig_edges, include_lowest=True)

        y_nat = np.expm1(tmp[gene].to_numpy()) if y_is_log1p else tmp[gene].to_numpy()
        tmp["_y_nat"] = y_nat

        mean_tbl = tmp.pivot_table(
            index="sig_bin",
            columns="dist_bin",
            values="_y_nat",
            aggfunc=aggfunc,
            observed=True,
        ).reindex(index=sig_cats, columns=dist_cats)

        mean_tables.append(mean_tbl)

        A = mean_tbl.to_numpy(dtype=float)
        A_col = _apply_color_scale(A, color_scale=color_scale, eps=eps, symlog_linthresh=symlog_linthresh)
        colored_arrays.append(A_col)

    # ---- 3) normalization (global vs per gene) ----
    global_cmin = global_cmax = None
    if norm_scope == "global":
        stacked = np.concatenate([a.ravel() for a in colored_arrays])
        stacked = stacked[np.isfinite(stacked)]
        if stacked.size == 0:
            global_cmin, global_cmax = 0.0, 1.0
        else:
            global_cmin, global_cmax = float(stacked.min()), float(stacked.max())
            if global_cmin == global_cmax:
                global_cmax = global_cmin + 1e-12

    # ---- 4) plotly subplot grid ----
    # viridis: use plotly's built-in viridis colorscale
    colorscale = pc.sequential.Viridis if str(cmap).lower() == "viridis" else pc.sequential.Viridis

    subplot_titles = list(y_cols)
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.03,
        vertical_spacing=0.07,
    )

    # We'll show a single shared colorbar only if share_colorbar=True.
    show_scale_first = bool(share_colorbar)

    n_x = len(dist_cats)
    n_y = len(sig_cats)

    for idx, (gene, mean_tbl, A_col) in enumerate(zip(y_cols, mean_tables, colored_arrays)):
        r = idx // ncols + 1
        c = idx % ncols + 1

        # per-gene scaling
        if norm_scope == "per_gene":
            finite = A_col[np.isfinite(A_col)]
            if finite.size == 0:
                cmin, cmax = 0.0, 1.0
            else:
                cmin, cmax = float(finite.min()), float(finite.max())
                if cmin == cmax:
                    cmax = cmin + 1e-12
        else:
            cmin, cmax = global_cmin, global_cmax

        # Provide hover in natural scale (mean_tbl), while coloring uses A_col
        V_nat = mean_tbl.to_numpy(dtype=float)

        # customdata holds natural values for hover
        customdata = V_nat

        heat = go.Heatmap(
            z=A_col,
            zmin=cmin,
            zmax=cmax,
            colorscale=colorscale,
            showscale=(show_scale_first and idx == 0),
            colorbar=dict(
                title=(
                    f"{aggfunc}(mean count)" if color_scale == "linear"
                    else f"{color_scale}({aggfunc}(mean count))"
                ),
                len=0.85,
            ) if (show_scale_first and idx == 0) else None,
            customdata=customdata,
            hovertemplate=(
                "Dist=%{x}<br>"
                "Sig=%{y}<br>"
                f"{gene}=%{{customdata:.{annot_decimals}f}}<extra></extra>"
            ),
        )
        fig.add_trace(heat, row=r, col=c)

        # Axes ticks / labels
        fig.update_xaxes(
            tickmode="array",
            tickvals=list(range(n_x)),
            ticktext=xlabels,
            tickangle=45,
            row=r,
            col=c,
            showgrid=False,
            zeroline=False,
            showline=True,
            mirror=True,
            linecolor="rgba(0,0,0,0.65)",
        )
        fig.update_yaxes(
            tickmode="array",
            tickvals=list(range(n_y)),
            ticktext=ylabels,
            row=r,
            col=c,
            showgrid=False,
            zeroline=False,
            showline=True,
            mirror=True,
            linecolor="rgba(0,0,0,0.65)",
            autorange="reversed",  # match imshow default top-to-bottom index order? (if you want origin="lower", remove)
        )

        # ---- cell gridlines (minor-grid look) as shapes (on top) ----
        # Draw lines at cell boundaries: -0.5, 0.5, ..., n-0.5
        # This mimics your dashed minor grid (Plotly doesn't do dashed cell grid perfectly)
        # Keeping alpha low like mpl (0.25).
        for xx in np.arange(-0.5, n_x, 1.0):
            fig.add_shape(
                type="line",
                x0=xx, x1=xx,
                y0=-0.5, y1=n_y - 0.5,
                xref=f"x{idx+1}" if idx > 0 else "x",
                yref=f"y{idx+1}" if idx > 0 else "y",
                line=dict(color="rgba(0,0,0,0.25)", width=1),
                layer="above",
            )
        for yy in np.arange(-0.5, n_y, 1.0):
            fig.add_shape(
                type="line",
                x0=-0.5, x1=n_x - 0.5,
                y0=yy, y1=yy,
                xref=f"x{idx+1}" if idx > 0 else "x",
                yref=f"y{idx+1}" if idx > 0 else "y",
                line=dict(color="rgba(0,0,0,0.25)", width=1),
                layer="above",
            )

        # ---- annotations (numbers) ----
        if annotate:
            for i in range(n_y):
                for j in range(n_x):
                    val = V_nat[i, j]
                    if not np.isfinite(val):
                        continue
                    fig.add_annotation(
                        x=j,
                        y=i,
                        xref=f"x{idx+1}" if idx > 0 else "x",
                        yref=f"y{idx+1}" if idx > 0 else "y",
                        text=f"{val:.{annot_decimals}f}",
                        showarrow=False,
                        font=dict(size=annot_fontsize, color="black"),
                    )

    # ---- global labels + title + transparent bg ----
    fig.update_layout(
        title=dict(
            text=f"Interaction maps (color_scale={color_scale}, norm_scope={norm_scope}, q={q})",
            x=0.5,
            y=0.98,
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=60, r=40, t=90, b=70),
        width=int(figsize[0] * 90),   # inches -> px rough
        height=int(figsize[1] * 90),
    )

    # "supxlabel"/"supylabel" equivalents as figure annotations
    fig.add_annotation(
        text="Distance bin",
        x=0.5, y=-0.02, xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=14),
    )
    fig.add_annotation(
        text="Signature bin",
        x=0.02, y=0.5, xref="paper", yref="paper",
        showarrow=False,
        textangle=-90,
        font=dict(size=14),
    )

    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)

    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config={"responsive": True, "displayModeBar": True, "displaylogo": False},
        auto_open=False,
    )

    if show:
        fig.show()

    return dist_edges, sig_edges

    # ---------- helpers (same idea as before, with target_gene support) ----------

def _guess_col(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Couldn't find {what} column. Tried {candidates}. Have {list(df.columns)}")


def extremes_by_model(
    leakage_all: pd.DataFrame,
    model_order: list[str],
    value_col: str,
) -> pd.DataFrame:
    model_col = _guess_col(
        leakage_all,
        ["model", "model_name", "estimator", "method", "regressor", "classifier"],
        "model",
    )
    gene_col = _guess_col(
        leakage_all,
        [
            "gene", "gene_name", "gene_id", "symbol", "ensg", "ENSG", "Gene", "GeneID",
            "target_gene", "target", "y_gene",
        ],
        "gene",
    )
    if value_col not in leakage_all.columns:
        raise ValueError(f"leakage_all must contain '{value_col}'. Have {list(leakage_all.columns)}")

    rows = []
    for m in model_order:
        sub = leakage_all.loc[leakage_all[model_col] == m, [gene_col, value_col]].dropna()
        if sub.empty:
            rows.append((m, np.nan, np.nan, np.nan, np.nan))
            continue
        imin = sub[value_col].idxmin()
        imax = sub[value_col].idxmax()
        rows.append(
            (
                m,
                sub.loc[imin, gene_col], float(sub.loc[imin, value_col]),
                sub.loc[imax, gene_col], float(sub.loc[imax, value_col]),
            )
        )

    return (
        pd.DataFrame(rows, columns=["model", "min_gene", "min_value", "max_gene", "max_value"])
        .set_index("model")
    )


# ---------- symlog-like transform for Plotly (workaround) ----------

def _symlog_forward(y: np.ndarray, linthresh: float = 1.0) -> np.ndarray:
    """
    Map y -> symlog space:
      sign(y) * ( log10(|y|/linthresh + 1) )
    Linear-ish around 0; log away from 0.
    """
    y = np.asarray(y, dtype=float)
    s = np.sign(y)
    a = np.abs(y)
    return s * np.log10(a / linthresh + 1.0)


def _symlog_ticks(linthresh: float, y_min: float, y_max: float) -> tuple[list[float], list[str]]:
    """
    Build tick positions in symlog-space, labeled with original y values.
    """
    # candidate ticks in original space (symmetric)
    # include 0, +/- linthresh, and log decades beyond linthresh
    candidates = [0.0, -linthresh, linthresh]

    # add decades up to max magnitude
    max_mag = max(abs(y_min), abs(y_max), linthresh)
    if max_mag > linthresh:
        top_dec = int(np.ceil(np.log10(max_mag)))
        # generate 1, 10, 100... * linthresh
        for p in range(0, top_dec + 1):
            v = (10**p) * linthresh
            candidates.extend([-v, v])

    # filter to range and unique-sort
    candidates = sorted(set([v for v in candidates if y_min <= v <= y_max]))
    tickvals = [_symlog_forward(v, linthresh) for v in candidates]
    ticktext = [f"{v:g}" for v in candidates]
    return tickvals, ticktext


def _symlog_limits_with_padding(values: list[float], linthresh: float, upper_mult: float = 2.2, lower_mult: float = 2.2):
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if v.size == 0:
        return (-1.0, 1.0), (-1.0, 1.0)

    lo = float(v.min())
    hi = float(v.max())

    # pad multiplicatively in absolute space (like your helper intent)
    # ensure we don't collapse when near 0
    lo_pad = lo * lower_mult if lo < 0 else lo / lower_mult
    hi_pad = hi * upper_mult if hi > 0 else hi / upper_mult

    # if all values are positive/negative, keep 0 in view a bit
    if lo >= 0:
        lo_pad = min(0.0, lo_pad)
    if hi <= 0:
        hi_pad = max(0.0, hi_pad)

    # convert to symlog space for plot axis range
    return (lo_pad, hi_pad), (_symlog_forward(lo_pad, linthresh), _symlog_forward(hi_pad, linthresh))


# ---------- main plot (call signature matches original) ----------

def plot_mean_relative_gap_with_extremes(
    agg: pd.DataFrame,
    leakage_all: pd.DataFrame,
    model_order: list[str],
    *,
    rel_gap_col: str = "Relative leakage gap, %",
    title: str = "Relative leakage gap, averaged over PIGs",
    yscale: str = "symlog",
    symlog_linthresh: float = 1.0,
    annotate_extremes: bool = True,
    fontsize: int = 8,
    figsize=(12, 5),
    out_html: str | Path = "frontend/public/plots/relative_variance_gap.html",
    show: bool = False,
):
    models = list(model_order)
    n = len(models)
    x = np.arange(n, dtype=float)

    # Align agg to model order
    if agg.index.name is None and "model" in agg.columns:
        agg2 = agg.set_index("model")
    else:
        agg2 = agg.copy()
    agg2 = agg2.reindex(models)

    if "rel_gap_pct_mean" not in agg2.columns or "rel_gap_pct_sem" not in agg2.columns:
        raise ValueError(
            "agg must contain columns ['rel_gap_pct_mean','rel_gap_pct_sem'] "
            f"but has {list(agg2.columns)}"
        )

    mean = agg2["rel_gap_pct_mean"].to_numpy(dtype=float)
    sem = agg2["rel_gap_pct_sem"].to_numpy(dtype=float)

    # Matplotlib default bar color (tab:blue)
    BAR_COLOR = "#1f77b4"
    SCATTER_COLOR = "black"

    # y-limits include mean±sem and extremes
    y_for_limits = []
    y_for_limits += list(mean - sem)
    y_for_limits += list(mean + sem)

    # extremes
    annotations = []
    ex_gap = None
    if annotate_extremes:
        ex_gap = extremes_by_model(leakage_all, models, rel_gap_col)
        y_for_limits += list(ex_gap["min_value"].to_numpy())
        y_for_limits += list(ex_gap["max_value"].to_numpy())

    fig = go.Figure()

    # Choose y mapping depending on scale
    if yscale == "log":
        if np.any(mean <= 0):
            raise ValueError("yscale='log' requires strictly positive mean relative gaps. Use yscale='symlog' instead.")
        y_plot = mean
        yerr_plot = sem
        yaxis_type = "log"
        tickvals = None
        ticktext = None
        y_range = None  
        hover_y = mean

    elif yscale == "symlog":
        # transform mean and (mean±sem) into symlog space; errorbars become asymmetric
        yaxis_type = "linear"  # we are in transformed space
        lo = mean - sem
        hi = mean + sem
        y_plot = _symlog_forward(mean, symlog_linthresh)
        yerr_plus = _symlog_forward(hi, symlog_linthresh) - y_plot
        yerr_minus = y_plot - _symlog_forward(lo, symlog_linthresh)

        # axis limits + ticks computed in original space then mapped
        (ymin0, ymax0), (ymin_t, ymax_t) = _symlog_limits_with_padding(
            y_for_limits, symlog_linthresh, upper_mult=2.2, lower_mult=2.2
        )
        y_range = [ymin_t, ymax_t]
        tickvals, ticktext = _symlog_ticks(symlog_linthresh, ymin0, ymax0)
        hover_y = mean

    else:
        raise ValueError("yscale must be one of {'log','symlog'}.")

    # Bars with error bars
    if yscale == "symlog":
        error_y = dict(
            type="data",
            symmetric=False,
            array=yerr_plus,
            arrayminus=yerr_minus,
            visible=True,
            thickness=1.2,
            width=6,
        )
    else:
        error_y = dict(type="data", array=yerr_plot, visible=True, thickness=1.2, width=6)

    fig.add_trace(
        go.Bar(
            x=x,
            y=y_plot,
            width=0.6,
            marker=dict(color=BAR_COLOR),
            error_y=error_y,
            name="Mean relative leakage gap",
            customdata=np.stack([hover_y], axis=1),
            hovertemplate="mean=%{customdata[0]:.3g}%<extra></extra>",
            showlegend=False,
        )
    )

    # Extremes points + annotations (plotted in same y-space)
    if annotate_extremes and ex_gap is not None:
        px, py, ptxt = [], [], []
        for i, m in enumerate(models):
            for which in ("min", "max"):
                v = ex_gap.loc[m, f"{which}_value"]
                g = ex_gap.loc[m, f"{which}_gene"]
                if pd.notna(v):
                    px.append(float(x[i]))
                    if yscale == "symlog":
                        py.append(float(_symlog_forward(v, symlog_linthresh)))
                    else:
                        py.append(float(v))
                    ptxt.append(f"{g}, {float(v):.1f}")

                    annotations.append(
                        dict(
                            x=float(x[i]),
                            y=(float(_symlog_forward(v, symlog_linthresh)) if yscale == "symlog" else float(v)),
                            xref="x",
                            yref="y",
                            text=f"{g}, {float(v):.1f}",
                            showarrow=True,
                            arrowhead=2,
                            arrowsize=1,
                            arrowwidth=1,
                            ax=0,
                            ay=-20 if which == "max" else 20,
                            font=dict(size=fontsize, color="black"),
                            bgcolor="rgba(255,255,255,0.70)",
                            bordercolor="rgba(0,0,0,0.25)",
                            borderwidth=1,
                        )
                    )

        if px:
            fig.add_trace(
                go.Scatter(
                    x=px,
                    y=py,
                    mode="markers",
                    marker=dict(size=7, color=SCATTER_COLOR),
                    showlegend=False,
                    customdata=np.array(ptxt, dtype=object),
                    hovertemplate="%{customdata}<extra></extra>",
                )
            )

    # Layout / axes styling (matplotlib-ish)
    fig.update_layout(
        title=title,
        margin=dict(l=70, r=30, t=70, b=150),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        annotations=annotations,
        width=int(figsize[0] * 90),
        height=int(figsize[1] * 90),
    )

    # X ticks
    fig.update_xaxes(
        tickmode="array",
        tickvals=x,
        ticktext=models,
        tickangle=25,
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
        range=[-0.6, (n - 1) + 0.6],
    )

    # Y axis
    yaxis_kwargs = dict(
        title="Mean relative leakage gap (%)",
        showgrid=True,
        gridcolor="rgba(0,0,0,0.25)",
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
        type=yaxis_type,
    )
    if yscale == "symlog":
        yaxis_kwargs["range"] = y_range
        yaxis_kwargs["tickmode"] = "array"
        yaxis_kwargs["tickvals"] = tickvals
        yaxis_kwargs["ticktext"] = ticktext

    fig.update_yaxes(**yaxis_kwargs)

    # Write HTML
    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config={"responsive": True, "displayModeBar": True, "displaylogo": False},
        auto_open=False,
    )

    if show:
        fig.show()

    return fig

def _guess_col(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Couldn't find {what} column. Tried {candidates}. Have {list(df.columns)}")


def extremes_by_model(
    leakage_all: pd.DataFrame,
    model_order: list[str],
    value_col: str,
) -> pd.DataFrame:
    """
    Auto-detect model/gene columns and compute per-model min/max (gene, value).
    Returns index=model with columns: min_gene, min_value, max_gene, max_value
    """

    # ✅ Updated: include your 'target_gene' column (and a few other common variants)
    model_col = _guess_col(
        leakage_all,
        ["model", "model_name", "estimator", "method", "regressor", "classifier"],
        "model",
    )
    gene_col = _guess_col(
        leakage_all,
        [
            "gene",
            "gene_name",
            "gene_id",
            "symbol",
            "ensg",
            "ENSG",
            "Gene",
            "GeneID",
            "target_gene",   # <-- YOUR COLUMN
            "target",        # common in ML pipelines
            "y_gene",        # sometimes used
        ],
        "gene",
    )

    if value_col not in leakage_all.columns:
        raise ValueError(f"leakage_all must contain '{value_col}'. Have {list(leakage_all.columns)}")

    rows = []
    for m in model_order:
        sub = leakage_all.loc[leakage_all[model_col] == m, [gene_col, value_col]].dropna()
        if sub.empty:
            rows.append((m, np.nan, np.nan, np.nan, np.nan))
            continue

        imin = sub[value_col].idxmin()
        imax = sub[value_col].idxmax()

        rows.append(
            (
                m,
                sub.loc[imin, gene_col],
                float(sub.loc[imin, value_col]),
                sub.loc[imax, gene_col],
                float(sub.loc[imax, value_col]),
            )
        )

    return (
        pd.DataFrame(rows, columns=["model", "min_gene", "min_value", "max_gene", "max_value"])
        .set_index("model")
    )



def _linear_limits_with_padding(values: list[float], pad_frac: float = 0.16) -> tuple[float, float]:
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if v.size == 0:
        return (-0.1, 0.1)
    lo = float(v.min())
    hi = float(v.max())
    if lo == hi:
        pad = 1.0 if hi == 0 else abs(hi) * pad_frac
        return (lo - pad, hi + pad)
    pad = (hi - lo) * pad_frac
    return (lo - pad, hi + pad)


def plot_mean_r2_with_extremes(
    agg: pd.DataFrame,
    leakage_all: pd.DataFrame,
    model_order: list[str],
    *,
    title: str = "Mean variance explained, averaged over PIGs",
    annotate_extremes: bool = True,
    fontsize: int = 8,
    figsize=(12, 5),
    # extra: output control (kept optional; call still matches original)
    out_html: str | Path = "frontend/public/plots/variance_spatial.html",
    show: bool = False,
):
    """
    Plotly interactive version of your matplotlib plot.

    Call signature matches the original function.
    Writes: frontend/public/plots/variance_spatial.html
    Transparent background, preserves style choices as closely as Plotly allows.
    """

    models = list(model_order)
    n = len(models)

    # Ensure agg is aligned
    # If agg is indexed by model names -> great
    # Otherwise, if it has a 'model' column, use it.
    if agg.index.name is None and "model" in agg.columns:
        agg2 = agg.set_index("model")
    else:
        agg2 = agg.copy()

    agg2 = agg2.reindex(models)

    required = [
        "R2_randomCV_mean",
        "R2_randomCV_sem",
        "R2_spatialBlockCV_mean",
        "R2_spatialBlockCV_sem",
    ]
    missing = [c for c in required if c not in agg2.columns]
    if missing:
        raise ValueError(f"agg missing columns {missing}. Have {list(agg2.columns)}")

    # Matplotlib geometry
    x = np.arange(n, dtype=float)
    w = 0.38
    xr = x - w / 2
    xs = x + w / 2

    rand_mean = agg2["R2_randomCV_mean"].to_numpy(dtype=float)
    rand_sem = agg2["R2_randomCV_sem"].to_numpy(dtype=float)
    spat_mean = agg2["R2_spatialBlockCV_mean"].to_numpy(dtype=float)
    spat_sem = agg2["R2_spatialBlockCV_sem"].to_numpy(dtype=float)

    fig = go.Figure()

    RAND_COLOR = "#1f77b4"
    SPAT_COLOR = "#ff7f0e"
    # Bars (keep Plotly default colors unless you want to hard-pin them)
    fig.add_trace(
        go.Bar(
            x=xr,
            y=rand_mean,
            width=w,
            name="Random cross-validation",
            marker=dict(color=RAND_COLOR),
            error_y=dict(type="data", array=rand_sem, visible=True, thickness=1.2, width=6),
            hovertemplate="Random CV<br>mean=%{y:.4f}<extra></extra>",
        )
    )
    fig.add_trace(
        go.Bar(
            x=xs,
            y=spat_mean,
            width=w,
            name="Spatial block cross-validation",
            marker=dict(color=SPAT_COLOR), 
            error_y=dict(type="data", array=spat_sem, visible=True, thickness=1.2, width=6),
            hovertemplate="Spatial block CV<br>mean=%{y:.4f}<extra></extra>",
        )
    )

    # y-limits should include mean±sem and extremes
    y_for_limits = []
    y_for_limits += list(rand_mean - rand_sem)
    y_for_limits += list(rand_mean + rand_sem)
    y_for_limits += list(spat_mean - spat_sem)
    y_for_limits += list(spat_mean + spat_sem)

    annotations = []

    if annotate_extremes:
        ex_rand = extremes_by_model(leakage_all, models, "R2_randomCV")
        ex_spat = extremes_by_model(leakage_all, models, "R2_spatialBlockCV")

        y_for_limits += list(ex_rand["min_value"].to_numpy())
        y_for_limits += list(ex_rand["max_value"].to_numpy())
        y_for_limits += list(ex_spat["min_value"].to_numpy())
        y_for_limits += list(ex_spat["max_value"].to_numpy())

        SCATTER_COLOR = "black"

        # Random extremes points at xr[i]
        rx, ry, rtxt = [], [], []
        sx, sy, stxt = [], [], []

        for i, m in enumerate(models):
            # Random min/max
            for which in ("min", "max"):
                v = ex_rand.loc[m, f"{which}_value"]
                g = ex_rand.loc[m, f"{which}_gene"]
                if pd.notna(v):
                    rx.append(float(xr[i]))
                    ry.append(float(v))
                    rtxt.append(f"{g}, {float(v):.4f}")

                    annotations.append(
                        dict(
                            x=float(xr[i]),
                            y=float(v),
                            xref="x",
                            yref="y",
                            text=f"{g}, {float(v):.4f}",
                            showarrow=True,
                            arrowhead=2,
                            arrowsize=1,
                            arrowwidth=1,
                            ax=-20,
                            ay=-20 if which == "max" else 20,
                            font=dict(size=fontsize, color="black"),
                            bgcolor="rgba(255,255,255,0.70)",
                            bordercolor="rgba(0,0,0,0.25)",
                            borderwidth=1,
                        )
                    )

            # Spatial min/max
            for which in ("min", "max"):
                v = ex_spat.loc[m, f"{which}_value"]
                g = ex_spat.loc[m, f"{which}_gene"]
                if pd.notna(v):
                    sx.append(float(xs[i]))
                    sy.append(float(v))
                    stxt.append(f"{g}, {float(v):.4f}")

                    annotations.append(
                        dict(
                            x=float(xs[i]),
                            y=float(v),
                            xref="x",
                            yref="y",
                            text=f"{g}, {float(v):.4f}",
                            showarrow=True,
                            arrowhead=2,
                            arrowsize=1,
                            arrowwidth=1,
                            ax=20,
                            ay=-20 if which == "max" else 20,
                            font=dict(size=fontsize, color="black"),
                            bgcolor="rgba(255,255,255,0.70)",
                            bordercolor="rgba(0,0,0,0.25)",
                            borderwidth=1,
                        )
                    )

        if rx:
            fig.add_trace(
                go.Scatter(
                    x=rx,
                    y=ry,
                    mode="markers",
                    marker=dict(size=7, color=SCATTER_COLOR),
                    showlegend=False,
                    hovertemplate="Random extreme<br>%{y:.4f}<extra></extra>",
                )
            )
        if sx:
            fig.add_trace(
                go.Scatter(
                    x=sx,
                    y=sy,
                    mode="markers",
                    marker=dict(size=7, color=SCATTER_COLOR),
                    showlegend=False,
                    hovertemplate="Spatial extreme<br>%{y:.4f}<extra></extra>",
                )
            )

    ylo, yhi = _linear_limits_with_padding(y_for_limits, pad_frac=0.16)

    # Tick labels (rotated like matplotlib)
    fig.update_layout(
        title=title,
        barmode="overlay",  # we’re manually positioning bars with numeric x
        legend=dict(x=1.0, y=1.0, xanchor="right", yanchor="top"),
        margin=dict(l=70, r=30, t=70, b=150),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        annotations=annotations,
        width=int(figsize[0] * 90),   # rough mapping inches->px
        height=int(figsize[1] * 90),
    )

    fig.update_xaxes(
        tickmode="array",
        tickvals=x,
        ticktext=models,
        tickangle=25,
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
        range=[-0.6, (n - 1) + 0.6],
    )

    fig.update_yaxes(
        title="Mean OOF R² across genes",
        range=[ylo, yhi],
        showgrid=True,
        gridcolor="rgba(0,0,0,0.25)",
        showline=True,
        mirror=True,
        linecolor="black",
        zeroline=False,
    )

    # Write HTML
    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(out_html),
        full_html=True,
        include_plotlyjs="cdn",
        config={"responsive": True, "displayModeBar": True, "displaylogo": False},
        auto_open=False,
    )

    if show:
        fig.show()

    return fig