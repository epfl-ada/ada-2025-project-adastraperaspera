from collections.abc import Sequence
import logging
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import gaussian_kde

# Get current directory
vis_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.dirname(vis_dir)
src_dir = os.path.dirname(scripts_dir)
figures_dir = os.path.join(src_dir, "data", "figures")


import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from typing import Optional, List

def plot_gene_trends_interactive(
    mean_expr: pd.DataFrame,
    genes: List[str],
    ylabel: str = "Mean expression (log₁₊ normalized)",
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

    # We’ll add 2 traces per broad type per gene: line (mean) + band (±SEM).
    # Use legendgroup and showlegend only once to avoid duplicate legend entries.
    traces_per_gene = 2 * len(btypes)  # band + line for each type

    for gi, gene in enumerate(PIGS):
        sub = agg[agg["gene"] == gene]
        for bi, bt in enumerate(btypes):
            dsub = sub[sub["broad_type"] == bt]
            if dsub.empty:
                # Skip empty groups (avoid blank traces)
                continue

            xcats = dsub["distance_bin"].astype(str)

            # Error band (invisible in legend)
            fig.add_trace(
                go.Scatter(
                    x=pd.concat([xcats, xcats[::-1]]),
                    y=pd.concat([dsub["mean"] + dsub["sem"], (dsub["mean"] - dsub["sem"])[::-1]]),
                    mode="lines",
                    fill="toself",
                    line=dict(width=0),
                    name=f"{bt} ± SEM",
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

    # Dropdown buttons (toggle visibility blocks per gene)
    buttons = []
    total_traces = len(fig.data)
    for gi, gene in enumerate(PIGS):
        vis = [False] * total_traces
        # find how many traces actually exist for this gene
        # we added traces sequentially: groups of traces per gene in the same order
        start = gi * traces_per_gene
        # But some btypes may be missing for a gene; safer approach: toggle by name in data slice
        # Simple approach if your dataset is dense:
        for idx in range(traces_per_gene):
            k = start + idx
            if k < total_traces:
                vis[k] = True
        buttons.append(
            dict(
                label=gene,
                method="update",
                args=[
                    {"visible": vis},
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
    title: str = "Plaque-Induced Gene Expression vs Distance (mean ± SEM, log₁₊)",
    x_label: str = "Distance to Plaque (µm, binned)",
    y_label: str = "Mean log₁₊ Expression",
    use_ci95: bool = True,   # multiply SEM by 1.96
):
    """
    Interactive line plot with gene selector and SEM bands.
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