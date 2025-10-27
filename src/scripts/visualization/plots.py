import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import plotly.graph_objects as go
from typing import List, Optional
import pandas as pd
from typing import Sequence, Optional
import plotly.graph_objects as go
from scipy.stats import gaussian_kde

def add_plaques_to_plotly(fig, plaques_gdf, name="Plaques", line_color="lime", line_width=2):
    def _add_ring(ring, fig):
        xs, ys = ring.xy  # -> array('d', ...)
        fig.add_trace(go.Scatter(
            x=list(xs), y=list(ys),  # conversion indispensable
            mode="lines",
            line=dict(color=line_color, width=line_width),
            name=name,
            hoverinfo="skip",
            showlegend=False  # évite de répéter la légende pour chaque polygone
        ))

    for geom in plaques_gdf.geometry.dropna():
        if geom.is_empty:
            continue
        gtype = geom.geom_type
        if gtype == "Polygon":
            _add_ring(geom.exterior, fig)
            for interior in geom.interiors:  # trous éventuels
                _add_ring(interior, fig)
        elif gtype == "MultiPolygon":
            for poly in geom.geoms:
                _add_ring(poly.exterior, fig)
                for interior in poly.interiors:
                    _add_ring(interior, fig)
        elif gtype in ("LineString", "LinearRing"):
            _add_ring(geom, fig)
        else:
            # fallback: tenter d'extraire un exterior s'il existe
            if hasattr(geom, "exterior") and geom.exterior is not None:
                _add_ring(geom.exterior, fig)

    return fig

def plot_plaques_dist(cells_with_dist, plaques_gdf):
    fig = px.scatter(
    cells_with_dist,
    x="x_centroid", y="y_centroid",
    color="distance_to_plaque",
    color_continuous_scale="plasma",
    title="Cell-to-plaque distance map",
    width=900, height=800,
    render_mode="webgl"  # accélère l’affichage si beaucoup de points
    )

    fig.update_traces(marker=dict(size=3), selector=dict(mode='markers'))
    fig = add_plaques_to_plotly(fig, plaques_gdf, name="Plaques", line_color="lime", line_width=2)

    fig.update_yaxes(scaleanchor="x", scaleratio=1)

    # export
    fig.write_html("src/data/figures/cell_to_plaque_map_interactive.html")
    #fig.write_image("src/data/figures/cell_to_plaque_map_interactive.png", scale=2)
    fig.show()

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
) -> go.Figure:
    """
    Version Plotly interactive de 'plot_gene_trends'.

    Paramètres
    ----------
    mean_expr : DataFrame
        Index = bins (IntervalIndex, Categorical ou labels), colonnes = gènes.
    genes : list[str]
        Sous-ensemble de gènes à tracer (ignore ceux absents).
    """
    if mean_expr.empty:
        raise ValueError("mean_expr is empty; check your inputs")

    # Crée des labels propres pour l'axe X (gère IntervalIndex -> 'L-R')
    if isinstance(mean_expr.index, pd.IntervalIndex):
        bin_labels = [f"{b.left:.0f}-{b.right:.0f}" for b in mean_expr.index]
    else:
        bin_labels = [str(x) for x in mean_expr.index]

    # Assemble en format long pour Plotly
    df = mean_expr.copy()
    df["__bin__"] = bin_labels
    long = df.melt(id_vars="__bin__", var_name="gene", value_name="mean_expr")

    # Filtre sur les gènes demandés et garde seulement les présents
    genes_present = [g for g in genes if g in mean_expr.columns]
    if not genes_present:
        raise ValueError("None of the requested genes were found in 'mean_expr' columns.")

    long = long[long["gene"].isin(genes_present)]

    # Trace (line + markers)
    if use_webgl:
        # Scattergl (via graph_objects) -> plus fluide si beaucoup de points/traces
        fig = go.Figure()
        for g in genes_present:
            sub = long[long["gene"] == g]
            fig.add_trace(
                go.Scattergl(
                    x=sub["__bin__"], y=sub["mean_expr"],
                    mode="lines+markers",
                    name=g,
                    line=dict(width=line_width),
                )
            )
    else:
        fig = px.line(
            long, x="__bin__", y="mean_expr", color="gene",
            markers=True, title=title, height=height, width=width
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
    fig.update_xaxes(tickangle=45)
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
    Version Plotly interactive de 'plot_mean_heatmap'.

    - Sélectionne les 'top_n' gènes qui varient le plus (somme des |diff| entre bins).
    - Z-score optionnel par gène pour comparer les profils.
    """
    non_gene_cols = {
        "cell_id","x_centroid","y_centroid","cell_area","nucleus_area",
        "total_counts","transcript_counts","distance_to_plaque","distance_bin",
    }
    gene_cols = [c for c in mean_expr.columns if c not in non_gene_cols]
    if not gene_cols:
        raise ValueError("No gene columns detected in 'mean_expr'.")

    expr = mean_expr[gene_cols].copy()

    # Sélection des gènes avec plus forte variation spatiale
    grad = expr.diff().abs().sum().sort_values(ascending=False)
    top_genes = grad.head(top_n).index
    sub_df = expr[top_genes]

    # Z-score optionnel
    if zscore:
        # convertit sparse -> dense si besoin, puis zscore colonne par colonne
        sub_df = sub_df.apply(
            lambda x: x.sparse.to_dense() if pd.api.types.is_sparse(x) else x
        )
        sub_df = (sub_df - sub_df.mean()) / (sub_df.std(ddof=0).replace(0, np.nan))

    # Labels d'axe X
    if isinstance(mean_expr.index, pd.IntervalIndex):
        x_labels = [f"{b.left:.0f}-{b.right:.0f}" for b in mean_expr.index]
    else:
        x_labels = [str(x) for x in mean_expr.index]

    # Plotly attend un array 2D -> transpose pour avoir (gene x bin)
    fig = px.imshow(
        sub_df.T.values,
        x=x_labels,
        y=sub_df.columns,
        color_continuous_scale="RdBu" if zscore else "Magma",
        origin="upper",
        aspect="auto",
        labels=dict(color="Z-score" if zscore else "Mean log1p"),
        title=title,
        height=height, width=width,
    )
    # Centrer la palette à 0 si zscore
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
    kde_points: int = 400
):
    """
    Create an interactive HTML figure with dropdowns to choose:
      - the gene
      - the scale (raw vs log1p)

    df: DataFrame (cells x genes) with numeric columns for gene expression.
    genes: list of gene column names in df.
    """
    # --- Precompute per-gene traces (raw + log1p) ---
    traces = []      # list of go.Scatter / go.Histogram
    vis_map = {}     # (gene, scale) -> list of trace indices to set visible=True
    x_ranges = {}    # scale -> (global_min, global_max) for consistent axes

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
            x=x_raw, nbinsx=bins, histnorm="probability density",
            name=f"{g} — hist (raw)", opacity=0.45, showlegend=False
        )
        # KDE (raw)
        raw_kde_trace = None
        if x_raw.size > 5:
            xr = np.linspace(max(x_ranges["raw"][0], np.min(x_raw)),
                             min(x_ranges["raw"][1], np.max(x_raw)),
                             kde_points)
            try:
                kde = gaussian_kde(x_raw[x_raw > 0] if (x_raw > 0).sum() > 5 else x_raw)
                yr = kde(xr)
                raw_kde_trace = go.Scatter(
                    x=xr, y=yr, mode="lines",
                    name=f"{g} — kde (raw)", line=dict(width=2), showlegend=False
                )
            except Exception:
                pass

        # Histogram (log1p)
        h_log = go.Histogram(
            x=x_log, nbinsx=bins, histnorm="probability density",
            name=f"{g} — hist (log1p)", opacity=0.45, showlegend=False
        )
        # KDE (log1p)
        log_kde_trace = None
        if x_log.size > 5:
            xl = np.linspace(max(x_ranges["log1p"][0], np.min(x_log)),
                             min(x_ranges["log1p"][1], np.max(x_log)),
                             kde_points)
            try:
                kde_l = gaussian_kde(x_log)  # déjà > 0
                yl = kde_l(xl)
                log_kde_trace = go.Scatter(
                    x=xl, y=yl, mode="lines",
                    name=f"{g} — kde (log1p)", line=dict(width=2), showlegend=False
                )
            except Exception:
                pass

        # Store indices for visibility toggling
        start_idx = len(traces)
        g_raw_idxs = [start_idx]                     # raw hist
        traces.append(h_raw)
        if raw_kde_trace is not None:
            g_raw_idxs.append(len(traces))
            traces.append(raw_kde_trace)

        g_log_idxs = [len(traces)]                   # log hist
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
        gene_buttons.append(dict(
            label=g,
            method="update",
            args=[
                {"visible": visibility_for(g, init_scale)},
                {"title": f"{title} — {g} ({init_scale})",
                 "xaxis": {"title": "log1p(expression)"} if init_scale == "log1p" else {"title": "expression"}}
            ],
        ))

    # Buttons for scale
    scale_buttons = []
    for sc in ["raw", "log1p"]:
        scale_buttons.append(dict(
            label=sc,
            method="update",
            args=[
                {"visible": visibility_for(init_gene, sc)},
                {"title": f"{title} — {init_gene} ({sc})",
                 "xaxis": {"title": "log1p(expression)"} if sc == "log1p" else {"title": "expression"}}
            ],
        ))

    fig.update_layout(
        width=width, height=height, template="simple_white",
        title=f"{title} — {init_gene} ({init_scale})",
        xaxis_title="log1p(expression)",
        yaxis_title="density",
        barmode="overlay",
        legend_title=None,
        updatemenus=[
            dict(
                buttons=gene_buttons,
                direction="down", showactive=True, x=0.02, xanchor="left", y=1.15, yanchor="top",
                bgcolor="white", bordercolor="#ccc"
            ),
            dict(
                buttons=scale_buttons,
                direction="down", showactive=True, x=0.30, xanchor="left", y=1.15, yanchor="top",
                bgcolor="white", bordercolor="#ccc",
            ),
        ],
        margin=dict(l=60, r=20, t=90, b=60),
    )

    # Consistent x ranges per scale (switch via relayout on button click)
    # We’ll attach ranges to layout meta for clarity (optional)
    fig.layout.meta = dict(xrange_raw=x_ranges["raw"], xrange_log=x_ranges["log1p"])
    fig.update_layout(
    title={
        "text": f"{title} — {init_gene} ({init_scale})",
        "x": 0.5,                # center horizontally
        "xanchor": "center",     
        "y": 0.97,               # slightly below top edge
        "yanchor": "top",
    },
    updatemenus=[
        dict(
            buttons=gene_buttons,
            direction="down",
            showactive=True,
            x=0.0, xanchor="left",
            y=1.12, yanchor="top",   # a bit below title
            bgcolor="white", bordercolor="#ccc"
        ),
        dict(
            buttons=scale_buttons,
            direction="down",
            showactive=True,
            x=0.25, xanchor="left",
            y=1.12, yanchor="top",
            bgcolor="white", bordercolor="#ccc"
        ),
    ],
    margin=dict(l=60, r=20, t=100, b=60),
    )

    return fig

def plot_top_spatial_genes_interactive(
    stats_df: pd.DataFrame,
    top_n: int = 20,
    metrics: Optional[List[str]] = None,   # ex: ["spearman_r","slope"]
    gene_col: str = "gene",
    p_col_candidates = ("p_value","pval","p"),
    fdr_col_candidates = ("fdr","q_value","adj_p","qval"),
    title: str = "Top 20 genes by spatial metric",
    width: int = 820,
    height: int = 650
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

    # Détection auto des métriques si non fournies
    if metrics is None:
        non_metric = {gene_col, *p_col_candidates, *fdr_col_candidates}
        metrics = [
            c for c in df.select_dtypes(include=[np.number]).columns
            if c not in non_metric
        ]
    if not metrics:
        raise ValueError("Aucune métrique numérique détectée. Fournis `metrics=[...]`.")

    # Colonnes p/FDR pour le hover (si présentes)
    p_col = next((c for c in p_col_candidates if c in df.columns), None)
    fdr_col = next((c for c in fdr_col_candidates if c in df.columns), None)

    # Construit 1 trace/barplot par métrique (on bascule la visibilité via dropdown)
    traces = []
    vis_map = {}
    top_n = int(min(top_n, len(df)))

    for m in metrics:
        if m not in df.columns:
            continue
        top = df.nlargest(top_n, m).copy().sort_values(m, ascending=True)  # pour empiler vers le haut

        hover = f"<b>%{{y}}</b><br>{m}: %{{x:.4g}}"
        custom = None
        if p_col or fdr_col:
            hover += f"<br>{p_col or 'p'}: %{{customdata[0]:.2e}}" if p_col else ""
            hover += f"<br>{fdr_col or 'FDR'}: %{{customdata[1]:.2e}}" if fdr_col else ""
            custom = np.stack([
                top[p_col].to_numpy() if p_col else np.full(len(top), np.nan),
                top[fdr_col].to_numpy() if fdr_col else np.full(len(top), np.nan),
            ], axis=1)

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

    # Boutons du menu (sélecteur de “gène” dans ton vocabulaire, ici = métrique)
    buttons = []
    for m in metrics:
        vis = [i == vis_map[m] for i in range(len(traces))]
        buttons.append(dict(
            label=m,
            method="update",
            args=[
                {"visible": vis},
                {"title": f"{title} — {m}", "xaxis": {"title": m}},
            ],
        ))

    fig.update_layout(
        width=width, height=height, template="simple_white",
        title=dict(text=f"{title} — {init_metric}", x=0.5, xanchor="center", y=0.96, yanchor="top"),
        xaxis_title=init_metric,
        yaxis_title="Gene",
        margin=dict(l=140, r=30, t=110, b=50),  # top ↑ pour loger le menu
        showlegend=False,
        updatemenus=[dict(
            buttons=buttons,
            direction="down", showactive=True,
            x=0.02, xanchor="left",
            y=1.10, yanchor="top",       # juste sous le titre
            bgcolor="white", bordercolor="#ccc",
            pad={"r": 6, "t": 6},
        )],
    )
    return fig

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
            hover += f"<br>p: %{{customdata[0]:.2e}}"
        if fdr_col:
            hover += f"<br>FDR: %{{customdata[1]:.2e}}"
        custom = np.stack([
            top[p_col].to_numpy() if p_col else np.full(len(top), np.nan),
            top[fdr_col].to_numpy() if fdr_col else np.full(len(top), np.nan),
        ], axis=1) if (p_col or fdr_col) else None

        new_data.append({
            "x": [top[m].to_numpy()],
            "y": [top["gene"].to_numpy()],
            "customdata": [custom] if custom is not None else [None],
            "hovertemplate": [hover],
            "marker": [dict(color=np.where(top[m].to_numpy() >= 0, "rgb(31,120,180)", "rgb(227,26,28)"))],
        })
    # visibility mask stays the same; layout title & xaxis will be set by the dropdown
    return [{"data": new_data}, {}]

def draw_figures(plot_func,img_pth="std.png", *args, **kwargs):
    fig = plot_func(*args, **kwargs)
    fig.write_html(img_pth)
    fig.show()
