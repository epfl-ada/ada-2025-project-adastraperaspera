from __future__ import annotations
from collections.abc import Sequence
import logging
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import gaussian_kde
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from typing import Optional, List
from pathlib import Path
from typing import Iterable
import base64
import numpy as np
from PIL import Image
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Get current directory
vis_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.dirname(vis_dir)
src_dir = os.path.dirname(scripts_dir)
figures_dir = os.path.join(src_dir, "data", "figures")


def plot_gene_trends_interactive(
    mean_expr: pd.DataFrame,
    genes: List[str],
    ylabel: str = "Mean expression (log1p normalized)",
    xlabel: str = "Distance to plaque (µm, binned)",
    title: str = "Spatial gene expression gradients",
    line_width: int = 2,
    height: int = 480,
    width: int = 820,
    use_webgl: bool = True,
    *,
    sem_expr: Optional[pd.DataFrame] = None,  # NEW: optional SEM matrix (same shape as mean_expr)
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
                    "Mean: %{y:.3f}" + ("<br>SEM: %{customdata:.3f}" if S is not None else "") +
                    "<extra></extra>"
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
    fig.update_xaxes(
        tickangle=45,
        categoryorder="array",
        categoryarray=bin_labels
    )

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
                    y=pd.concat([
                        dsub["mean"] + 1.96 * dsub["sem"],
                        (dsub["mean"] - 1.96 * dsub["sem"])[::-1],
                    ]),
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
    use_ci95: bool = True,   # multiply SEM by 1.96
):
    """
    Interactive line plot with gene selector and CI bands (±1.96×SEM ≈ 95% CI).
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

    buttons = []
    traces_per_gene = 2  # main line + shaded band

    for i, gene in enumerate(pig_genes):
        sub = df[df[gene_col] == gene].sort_values("bin_label")

        m = sub[mean_col].to_numpy(float)
        s = sub[sem_col].to_numpy(float) * scale
        x = sub["bin_label"].tolist()

        # mean line
        fig.add_trace(go.Scatter(
            x=x, y=m, mode="lines+markers",
            name=f"{gene}",
            visible=(i == 0),
            line=dict(width=2),
            marker=dict(size=7),
        ))

        # SEM (shaded band)
        fig.add_trace(go.Scatter(
            x=x + x[::-1],
            y=(m + s).tolist() + (m - s)[::-1].tolist(),
            fill="toself",
            fillcolor="rgba(31, 119, 180, 0.18)",
            line=dict(width=0),
            hoverinfo="skip",
            name=f"{gene} CI",
            visible=(i == 0),
        ))

        # button to toggle visibility
        vis = [False] * (len(pig_genes) * traces_per_gene)
        vis[i*traces_per_gene:(i+1)*traces_per_gene] = [True, True]

        buttons.append(dict(label=gene, method="update",
                            args=[{"visible": vis},
                                  {"title": f"{title}<br><sup>{gene}</sup>"}]))

    fig.update_layout(
        updatemenus=[dict(
            buttons=buttons,
            direction="down",
            x=0.5, xanchor="center",
            y=1.15, yanchor="top"
        )],
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        template="plotly_white",
        margin=dict(t=120, l=60, r=20, b=60),
    )

    fig.show()



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
    csv_paths: dict[str, str | Path],                 # e.g. {"wt2": "...csv", "tg2": "...csv", ...}
    age_map: dict[str, tuple[str, str]],              # age_label -> (wt_key, tg_key)
    title: str | None = "WT vs TG by age (interactive)",
    filename: str = "wt_tg_age_grid_scatter.html",
    out_dir: str = "frontend/public/plots",
    max_points: int | None = 150_000,                 # downsample for browser performance
    x_candidates=("x_centroid", "x"),
    y_candidates=("y_centroid", "y"),
    # If your tissue coordinates behave like images (origin top-left), this makes it look right:
    reverse_y: bool = True,
    # If x/y are truly swapped in your CSV, set this to True:
    swap_xy: bool = False,
    # styling
    marker_size: float = 1.8,
    marker_opacity: float = 0.65,
    color_col: str | None = None,                     # e.g. "prediction" (numeric) or None
    row_label_wt: str = "Wild type",
    row_label_tg: str = "Transgenic",
    row_label_font_size: int = 18,
) -> go.Figure:
    out_path = Path(out_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)

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

    def prep_df(key: str) -> tuple[pd.Series, pd.Series, pd.Series | None]:
        df = pd.read_csv(csv_paths[key])

        if max_points is not None and len(df) > max_points:
            df = df.sample(n=max_points, random_state=0)

        xcol, ycol = pick_xy(df)

        x = df[xcol].astype(float)
        y = df[ycol].astype(float)

        if swap_xy:
            x, y = y, x

        c = None
        if color_col is not None and color_col in df.columns:
            c = df[color_col]
        return x, y, c

    # fixed age order based on insertion order in age_map
    age_labels = list(age_map.keys())
    if len(age_labels) != 3:
        # not required, but your target layout is 3 columns; raise early if inconsistent
        raise ValueError("age_map should contain exactly 3 ages for a 2×3 grid.")

    # ---- create subplot grid ----
    fig = make_subplots(
        rows=2,
        cols=3,
        column_titles=[f"{a} months" if a.replace('.', '', 1).isdigit() else f"{a} months" for a in age_labels],
        horizontal_spacing=0.02,
        vertical_spacing=0.06,
    )

    # Track global ranges so all panels match (no jumping/unequal zoom)
    xmins, xmaxs, ymins, ymaxs = [], [], [], []

    # ---- add traces ----
    for j, age in enumerate(age_labels, start=1):
        wt_key, tg_key = age_map[age]

        # WT (row 1)
        x_wt, y_wt, c_wt = prep_df(wt_key)
        xmins.append(float(x_wt.min())); xmaxs.append(float(x_wt.max()))
        ymins.append(float(y_wt.min())); ymaxs.append(float(y_wt.max()))

        marker_wt = dict(size=marker_size, opacity=marker_opacity)
        if c_wt is not None and pd.api.types.is_numeric_dtype(c_wt):
            marker_wt["color"] = c_wt
            marker_wt["showscale"] = (j == 3)  # show colorbar only on last column

        fig.add_trace(
            go.Scattergl(
                x=x_wt,
                y=y_wt,
                mode="markers",
                marker=marker_wt,
                showlegend=False,
                hovertemplate=f"{row_label_wt}<br>Age: {age}<br>x=%{{x:.2f}}<br>y=%{{y:.2f}}<extra></extra>",
            ),
            row=1, col=j
        )

        # TG (row 2)
        x_tg, y_tg, c_tg = prep_df(tg_key)
        xmins.append(float(x_tg.min())); xmaxs.append(float(x_tg.max()))
        ymins.append(float(y_tg.min())); ymaxs.append(float(y_tg.max()))

        marker_tg = dict(size=marker_size, opacity=marker_opacity)
        if c_tg is not None and pd.api.types.is_numeric_dtype(c_tg):
            marker_tg["color"] = c_tg
            marker_tg["showscale"] = False  # already shown (if any) on WT last col

        fig.add_trace(
            go.Scattergl(
                x=x_tg,
                y=y_tg,
                mode="markers",
                marker=marker_tg,
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
                scaleanchor=f"x{'' if (r == 1 and c == 1) else ( (r-1)*3 + c )}",
                scaleratio=1,
                autorange="reversed" if reverse_y else True,
            )

    # ---- row labels (bigger + bold) ----
    # Add annotations on the left side, vertically centered per row
    fig.update_layout(
    annotations=list(fig.layout.annotations) + [
        dict(
            text=f"<b>{row_label_wt}</b>",
            x=0.01, y=0.97,                 # ⬅ inside the plot
            xref="paper", yref="paper",
            xanchor="left", yanchor="middle",
            showarrow=False,
            font=dict(size=row_label_font_size),
        ),
        dict(
            text=f"<b>{row_label_tg}</b>",
            x=0.01, y=0.50,                 # ⬅ inside the plot
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
        width=None,
        height=720,
        margin=dict(l=30, r=20, t=80 if title else 40, b=30),
        dragmode="pan",
        template="plotly_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig.write_html(str(out_path), include_plotlyjs="cdn")
    return fig