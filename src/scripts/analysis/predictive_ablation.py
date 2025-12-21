import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
from matplotlib import colors
from matplotlib.colors import PowerNorm
from matplotlib.patches import Rectangle
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold, KFold, train_test_split
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import KBinsDiscretizer, SplineTransformer, StandardScaler
from sklearn.tree import DecisionTreeRegressor


def plot_spatial_scatter(df_base, x_col, y_col, val_col, title, s=6):
    plt.figure()
    plt.scatter(df_base[x_col], df_base[y_col], c=df_base[val_col], s=s)
    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.colorbar(label=val_col)
    plt.tight_layout()
    plt.show()


def _fit_line_slope_intercept(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2:
        return np.nan, np.nan
    slope, intercept = np.polyfit(x[m], y[m], 1)
    return slope, intercept


def _plot_scatter_with_fit(ax, x, y, label, color, s=8, alpha=0.35, line_width=2.0):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ax.scatter(x, y, s=s, alpha=alpha, label=label, color=color)

    slope, intercept = _fit_line_slope_intercept(x, y)
    if np.isfinite(slope) and np.isfinite(intercept):
        xmin = np.nanmin(x)
        xmax = np.nanmax(x)
        xline = np.array([xmin, xmax], dtype=float)
        yline = intercept + slope * xline
        ax.plot(xline, yline, linewidth=line_width, color=color)

        dx = (xmax - xmin) if np.isfinite(xmax - xmin) else 0.0
        x_text = xmax + 0.015 * dx
        y_text = intercept + slope * xmax
        ax.text(
            x_text,
            y_text,
            f"{slope:.2f}",
            color=color,
            va="center",
            ha="left",
            fontsize=9,
        )


def plot_pred_vs_obs_combined(models, target_gene=None, title=""):
    """
    If target_gene is provided and exists in each model dataframe, that column is used as observed;
    otherwise, column "y" is used.
    """
    fig, ax = plt.subplots(figsize=(7, 6))

    colors_cycle = (
        plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0", "C1", "C2", "C3"])
    )
    for i, m in enumerate(models):
        name = m["name"]
        dfi = m["df"]
        r2 = m.get("r2", None)

        y_col = (
            target_gene
            if (target_gene is not None and target_gene in dfi.columns)
            else "y"
        )
        obs = dfi[y_col].to_numpy()
        pred = dfi["oof_pred"].to_numpy()

        lbl = f"{name}" if r2 is None else f"{name} (mean R²={r2:.3f})"
        _plot_scatter_with_fit(
            ax, obs, pred, lbl, color=colors_cycle[i % len(colors_cycle)]
        )

    ax.set_xlabel("Observed (mean across genes)")
    ax.set_ylabel("OOF Predicted (mean across genes)")
    ax.set_title(title)
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.show()


def plot_resid_vs_distance_combined(models, dist_col, title):
    fig, ax = plt.subplots(figsize=(7, 6))

    colors_cycle = (
        plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0", "C1", "C2", "C3"])
    )
    for i, m in enumerate(models):
        name = m["name"]
        dfi = m["df"]

        x = dfi[dist_col].to_numpy()
        y = dfi["oof_resid"].to_numpy()

        _plot_scatter_with_fit(
            ax, x, y, name, color=colors_cycle[i % len(colors_cycle)]
        )

    ax.set_xlabel(dist_col)
    ax.set_ylabel("OOF Residual (mean obs - mean pred)")
    ax.set_title(title)
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.show()


def tile_ids(df, x_col, y_col, n_tiles_x=10, n_tiles_y=10):
    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()
    x01 = (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x) + 1e-12)
    y01 = (y - np.nanmin(y)) / (np.nanmax(y) - np.nanmin(y) + 1e-12)
    tx = np.clip((x01 * n_tiles_x).astype(int), 0, n_tiles_x - 1)
    ty = np.clip((y01 * n_tiles_y).astype(int), 0, n_tiles_y - 1)
    return tx, ty


def plot_tile_heatmap(df_in, x_col, y_col, val_col, title, n_tiles_x=10, n_tiles_y=10):
    tx, ty = tile_ids(df_in, x_col, y_col, n_tiles_x, n_tiles_y)
    tmp = df_in.copy()
    tmp["tx"] = tx
    tmp["ty"] = ty
    mat = tmp.groupby(["ty", "tx"])[val_col].mean().unstack("tx")

    plt.figure(figsize=(7, 5.5))
    plt.imshow(mat.to_numpy(), aspect="auto")
    plt.title(title)
    plt.xlabel("tile x")
    plt.ylabel("tile y")
    plt.colorbar(label=f"mean({val_col})")
    plt.tight_layout()
    plt.show()


def true_features_for_gene(gene, neigh_cols, dist_col, optional_covs, true2fake):
    full_true, _ = _build_feature_lists_for_gene(
        gene, neigh_cols, dist_col, optional_covs, true2fake
    )
    return full_true


def fake_features_for_gene(gene, neigh_cols, dist_col, optional_covs, true2fake):
    _, full_fake = _build_feature_lists_for_gene(
        gene, neigh_cols, dist_col, optional_covs, true2fake
    )
    return full_fake


def spatial_oof_predictions(
    df, target_gene, feats, groups, n_splits=5, alpha=1.0, random_state=0
):
    """
    GroupKFold out-of-fold predictions and pooled OOF R^2.

    Returns:
      out_df: rows used after dropna (with __row_index__, oof_pred, oof_resid, fold)
      oof_r2: pooled/global OOF R^2
    """
    if len(groups) != len(df):
        raise ValueError("groups must be aligned with df rows (same length)")

    use = df[[target_gene] + list(feats)].copy()
    use["__group__"] = np.asarray(groups)
    use = use.dropna()

    if len(use) == 0:
        raise ValueError(
            f"No rows left after dropna for target_gene={target_gene} with feats={feats[:5]}..."
        )

    y = use[target_gene].to_numpy(dtype=float)
    X = use[list(feats)].to_numpy(dtype=float)
    g = use["__group__"].to_numpy()

    n_groups = len(np.unique(g))
    if n_groups < n_splits:
        raise ValueError(
            f"Not enough unique groups for GroupKFold: n_groups={n_groups} < n_splits={n_splits}"
        )

    cv = GroupKFold(n_splits=n_splits)

    oof_pred = np.full(len(use), np.nan, dtype=float)
    fold_id = np.full(len(use), -1, dtype=int)

    for k, (tr, te) in enumerate(cv.split(X, y, groups=g)):
        model = Pipeline(
            [
                ("scaler", StandardScaler(with_mean=True, with_std=True)),
                ("ridge", Ridge(alpha=alpha, random_state=random_state)),
            ]
        )
        model.fit(X[tr], y[tr])
        oof_pred[te] = model.predict(X[te])
        fold_id[te] = k

    if not np.isfinite(oof_pred).all():
        raise RuntimeError(
            "OOF predictions contain NaNs; check data/NA handling and GroupKFold coverage."
        )

    oof_r2 = float(r2_score(y, oof_pred))

    out = use.reset_index(drop=False).rename(columns={"index": "__row_index__"})
    out["oof_pred"] = oof_pred
    out["oof_resid"] = y - oof_pred
    out["fold"] = fold_id

    return out, oof_r2


def _build_feature_lists_for_gene(gene, neigh_cols, dist_col, optional_covs, true2fake):
    current_true_neighbors = [c for c in neigh_cols if gene not in c]
    full_true = [dist_col] + list(optional_covs) + current_true_neighbors

    current_fake_neighbors = [true2fake[c] for c in current_true_neighbors]
    full_fake = [dist_col] + list(optional_covs) + current_fake_neighbors

    return full_true, full_fake


def spatial_oof_predictions_multi(
    df,
    target_genes,
    groups,
    build_feats_fn,
    n_splits=5,
    alpha=1.0,
    random_state=0,
):
    """
    Runs spatial_oof_predictions for each gene.

    Returns:
      - oof_all: concatenated long-form OOF frame with columns: gene, y, oof_pred, oof_resid, fold, __row_index__, plus feats
      - r2_per_gene: DataFrame with per-gene R2
      - mean_r2: mean of per-gene R2
      - pooled_r2: pooled R2 across all gene observations
    """
    frames = []
    rows = []

    for gene in target_genes:
        feats = build_feats_fn(gene)
        out, r2 = spatial_oof_predictions(
            df,
            gene,
            feats,
            groups,
            n_splits=n_splits,
            alpha=alpha,
            random_state=random_state,
        )
        out = out.rename(columns={gene: "y"})
        out["gene"] = gene
        frames.append(out)

        rows.append({"gene": gene, "OOF_R2": float(r2)})

    r2_per_gene = pd.DataFrame(rows)
    mean_r2 = float(r2_per_gene["OOF_R2"].mean())

    oof_all = pd.concat(frames, ignore_index=True)
    pooled_r2 = float(
        r2_score(
            oof_all["y"].to_numpy(dtype=float),
            oof_all["oof_pred"].to_numpy(dtype=float),
        )
    )

    return oof_all, r2_per_gene, mean_r2, pooled_r2


def aggregate_oof_across_genes(
    oof_long,
    *,
    id_col="__row_index__",
    dist_col=None,
    require_all_genes=True,
    n_genes_expected=None,
):
    """
    Per-row aggregation by averaging across genes:
      y_mean, pred_mean, resid_mean.
    """
    if "gene" not in oof_long.columns:
        raise ValueError(
            "oof_long must have a 'gene' column (use spatial_oof_predictions_multi output)."
        )

    agg_dict = {
        "y": ("y", "mean"),
        "oof_pred": ("oof_pred", "mean"),
        "oof_resid": ("oof_resid", "mean"),
        "n_genes_used": ("gene", "nunique"),
    }
    if dist_col is not None and dist_col in oof_long.columns:
        agg_dict[dist_col] = (dist_col, "first")

    out = oof_long.groupby(id_col, as_index=False).agg(**agg_dict)

    if require_all_genes:
        if n_genes_expected is None:
            n_genes_expected = int(oof_long["gene"].nunique())
        out = out[out["n_genes_used"] == n_genes_expected].copy()

    return out


def build_fake_far_neighbor_means(df, x_col, y_col, gene_cols, k=15):
    coords = df[[x_col, y_col]].to_numpy()
    nn = NearestNeighbors(n_neighbors=min(len(df), k + 1), algorithm="auto")
    nn.fit(coords)
    dists, idxs = nn.kneighbors(coords)

    far_idxs = idxs[:, 1:][:, ::-1]
    far_idxs = far_idxs[:, :k]

    X_genes = df[gene_cols].to_numpy()
    fake_means = np.zeros((len(df), len(gene_cols)), dtype=float)
    for i in range(len(df)):
        fake_means[i] = X_genes[far_idxs[i]].mean(axis=0)

    fake_cols = [
        f"fake_far_neigh_mean_{c.replace('neigh_mean_', '')}" for c in gene_cols
    ]
    out = df.copy()
    for j, c in enumerate(fake_cols):
        out[c] = fake_means[:, j]
    return out, fake_cols


def mean_ci_t(x, confidence=0.95):
    x = np.asarray(x, dtype=float)
    n = len(x)
    m = float(np.mean(x))
    s = float(np.std(x, ddof=1))
    se = s / np.sqrt(n)
    tcrit = st.t.ppf((1 + confidence) / 2.0, df=n - 1)
    h = tcrit * se
    return m, m - h, m + h


def ceil_to_decimals(x, decimals=2):
    factor = 10**decimals
    return np.ceil(x * factor) / factor


def _oof_r2_spatial_blocks_ridge(
    df_in, target, feature_cols, groups_in, *, n_splits=5, alpha=1.0
):
    """
    Pooled OOF R^2 under GroupKFold using per-row out-of-fold predictions.
    """
    if len(groups_in) != len(df_in):
        raise ValueError("groups must be aligned with df rows (same length)")

    use = df_in[[target] + list(feature_cols)].copy()
    use["__group__"] = np.asarray(groups_in)
    use = use.dropna()

    if len(use) == 0:
        raise ValueError(f"No rows left after dropna for target={target}")

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
        model = Pipeline(
            [
                ("scaler", StandardScaler(with_mean=True, with_std=True)),
                ("ridge", Ridge(alpha=alpha, random_state=0)),
            ]
        )
        model.fit(X[tr], y[tr])
        pred = model.predict(X[te])

        oof_pred[te] = pred
        fold_r2s.append(r2_score(y[te], pred))

    if not np.isfinite(oof_pred).all():
        raise RuntimeError(
            "OOF predictions contain NaNs; check data and group splitting."
        )

    return float(r2_score(y, oof_pred)), fold_r2s


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
    k_far=15,
    n_splits=5,
    alpha=1.0,
    confidence=0.95,
    maroon="maroon",
    capsize=6,
    figsize=(9, 6),
):
    df_fake, fake_all_cols = build_fake_far_neighbor_means(
        df, x_col, y_col, neigh_cols, k=k_far
    )
    true2fake = {
        c: f"fake_far_neigh_mean_{c.replace('neigh_mean_', '')}" for c in neigh_cols
    }

    missing_fake = [
        true2fake[c] for c in neigh_cols if true2fake[c] not in df_fake.columns
    ]
    if missing_fake:
        raise RuntimeError(
            f"Missing fake columns in df_fake (unexpected): {missing_fake[:10]}"
        )

    rows = []
    for gene in target_genes:
        current_true_neighbors = [c for c in neigh_cols if gene not in c]
        current_fake_neighbors = [true2fake[c] for c in current_true_neighbors]

        feats_true = [dist_col] + list(optional_covs) + current_true_neighbors
        feats_fake = [dist_col] + list(optional_covs) + current_fake_neighbors

        r2_true, _ = _oof_r2_spatial_blocks_ridge(
            df, gene, feats_true, groups, n_splits=n_splits, alpha=alpha
        )
        r2_fake, _ = _oof_r2_spatial_blocks_ridge(
            df_fake, gene, feats_fake, groups, n_splits=n_splits, alpha=alpha
        )

        rows.append(
            {
                "gene": gene,
                "R2_true_neighbors": float(r2_true),
                "R2_fake_far_neighbors": float(r2_fake),
                "gap_true_minus_fake": float(r2_true - r2_fake),
            }
        )

    res_per_gene = pd.DataFrame(rows)

    m_true, lo_true, hi_true = mean_ci_t(
        res_per_gene["R2_true_neighbors"].to_numpy(), confidence=confidence
    )
    m_fake, lo_fake, hi_fake = mean_ci_t(
        res_per_gene["R2_fake_far_neighbors"].to_numpy(), confidence=confidence
    )

    avg_gap = float(np.mean(res_per_gene["gap_true_minus_fake"].to_numpy()))

    res_summary = pd.DataFrame(
        [
            {
                "condition": "true_neighbors",
                "mean_R2": m_true,
                "ci95_low": lo_true,
                "ci95_high": hi_true,
                "n_genes": len(res_per_gene),
            },
            {
                "condition": "fake_far_neighbors",
                "mean_R2": m_fake,
                "ci95_low": lo_fake,
                "ci95_high": hi_fake,
                "n_genes": len(res_per_gene),
            },
        ]
    )

    fig, ax = plt.subplots(figsize=figsize)

    labels = [
        "true neighbors",
        f"fake far neighbors (avg gap = {avg_gap:.3f})",
    ]

    means = np.array([m_true, m_fake], dtype=float)
    ci_lows = np.array([lo_true, lo_fake], dtype=float)
    ci_highs = np.array([hi_true, hi_fake], dtype=float)
    yerr = np.vstack([means - ci_lows, ci_highs - means])

    x = np.arange(2)
    ax.bar(x, means)

    ax.errorbar(
        x,
        means,
        yerr=yerr,
        fmt="none",
        ecolor="black",
        elinewidth=2,
        capsize=capsize,
        capthick=2,
        zorder=4,
    )

    def scatter_minmax(col, xi):
        sub = res_per_gene[["gene", col]].copy()
        imin = sub[col].idxmin()
        imax = sub[col].idxmax()

        for idx in [imin, imax]:
            gene = sub.loc[idx, "gene"]
            r2 = float(sub.loc[idx, col])
            r2_disp = float(ceil_to_decimals(r2, 2))

            ax.scatter([xi], [r2], color=maroon, s=60, zorder=5)
            ax.text(
                xi + 0.03,
                r2 + 0.002,
                f"{gene} ({r2_disp:.2f})",
                color=maroon,
                fontsize=9,
                ha="left",
                va="bottom",
                zorder=6,
            )

    scatter_minmax("R2_true_neighbors", 0)
    scatter_minmax("R2_fake_far_neighbors", 1)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(
        f"OOF $R^2$ (spatial-block GroupKFold), mean across {len(res_per_gene)} genes"
    )
    ax.set_title(
        f"True vs fake far-neighbor features (mean ± {int(confidence * 100)}% CI)"
    )

    ax.grid(True, axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show()

    return res_per_gene, res_summary


def permute_within_groups(df, cols, groups, seed=0):
    rng = np.random.default_rng(seed)
    out = df.copy()
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        if len(idx) <= 2:
            continue
        for c in cols:
            out.iloc[idx, out.columns.get_loc(c)] = rng.permutation(
                out.iloc[idx][c].to_numpy()
            )
    return out


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
):
    """
    Computes pooled OOF R^2 under spatial-block CV (GroupKFold) for true vs. permuted neighbors.
    """
    if optional_covs is None:
        optional_covs = []

    df_perm = permute_within_groups(df, neigh_cols, groups, seed=seed)

    rows = []
    for gene in target_genes:
        current_neighbors = [c for c in neigh_cols if gene not in c]
        full_feats = [dist_col] + list(optional_covs) + current_neighbors

        r2_true, _ = _oof_r2_spatial_blocks_ridge(
            df, gene, full_feats, groups, n_splits=n_splits, alpha=alpha
        )
        r2_perm, _ = _oof_r2_spatial_blocks_ridge(
            df_perm, gene, full_feats, groups, n_splits=n_splits, alpha=alpha
        )

        pct_drop = 100.0 * (r2_true - r2_perm) / (np.abs(r2_true) + pct_eps)

        rows.append(
            {
                "gene": gene,
                "R2_true": float(r2_true),
                "R2_permuted": float(r2_perm),
                "drop_abs": float(r2_true - r2_perm),
                "drop_pct": float(pct_drop),
            }
        )

    res_per_gene = pd.DataFrame(rows)

    m_true, lo_true, hi_true = mean_ci_t(
        res_per_gene["R2_true"].to_numpy(), confidence=confidence
    )
    m_perm, lo_perm, hi_perm = mean_ci_t(
        res_per_gene["R2_permuted"].to_numpy(), confidence=confidence
    )

    avg_pct_drop = float(np.mean(res_per_gene["drop_pct"].to_numpy()))

    res_summary = pd.DataFrame(
        [
            {
                "condition": "true_neighbors",
                "mean_R2": m_true,
                "ci95_low": lo_true,
                "ci95_high": hi_true,
                "n_genes": len(res_per_gene),
            },
            {
                "condition": "permuted_neighbors",
                "mean_R2": m_perm,
                "ci95_low": lo_perm,
                "ci95_high": hi_perm,
                "n_genes": len(res_per_gene),
            },
        ]
    )

    fig, ax = plt.subplots(figsize=figsize)

    labels = [
        "true neighbors",
        f"permuted neighbors (avg drop = {avg_pct_drop:.1f}%)",
    ]

    means = np.array([m_true, m_perm], dtype=float)
    ci_lows = np.array([lo_true, lo_perm], dtype=float)
    ci_highs = np.array([hi_true, hi_perm], dtype=float)
    yerr = np.vstack([means - ci_lows, ci_highs - means])
    x = np.arange(2)

    ax.bar(x, means)

    ax.errorbar(
        x,
        means,
        yerr=yerr,
        fmt="none",
        ecolor="black",
        elinewidth=2,
        capsize=capsize,
        capthick=2,
        zorder=4,
    )

    def scatter_minmax(condition_col, xi):
        sub = res_per_gene[["gene", condition_col]].copy()
        imin = sub[condition_col].idxmin()
        imax = sub[condition_col].idxmax()

        for idx in [imin, imax]:
            gene = sub.loc[idx, "gene"]
            r2 = float(sub.loc[idx, condition_col])
            r2_disp = float(ceil_to_decimals(r2, 2))

            ax.scatter([xi], [r2], color=maroon, s=60, zorder=5)
            ax.text(
                xi + 0.03,
                r2 + 0.002,
                f"{gene} ({r2_disp:.2f})",
                color=maroon,
                fontsize=9,
                ha="left",
                va="bottom",
                zorder=6,
            )

    scatter_minmax("R2_true", 0)
    scatter_minmax("R2_permuted", 1)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(
        f"OOF $R^2$ (spatial-block GroupKFold), mean across {len(res_per_gene)} genes"
    )
    ax.set_title(
        f"Permutation test within tiles (OOF $R^2$ mean ± {int(confidence * 100)}% CI)"
    )

    ax.grid(True, axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show()

    return res_per_gene, res_summary


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
):
    """
    res_per_gene must contain columns:
      - gene
      - model
      - R2_spatialBlockCV
    """
    required = {"gene", "model", "R2_spatialBlockCV"}
    missing = required - set(res_per_gene.columns)
    if missing:
        raise ValueError(f"res_per_gene missing columns: {missing}")

    summary_rows = []
    for model_name, sub in res_per_gene.groupby("model"):
        gene_means = sub["R2_spatialBlockCV"].to_numpy()
        m, lo, hi = mean_ci_t(gene_means, confidence=confidence)
        summary_rows.append(
            {
                "model": model_name,
                "mean_R2_across_genes": m,
                "ci95_low": lo,
                "ci95_high": hi,
                "n_genes": len(gene_means),
            }
        )

    res_summary = (
        pd.DataFrame(summary_rows)
        .sort_values("mean_R2_across_genes", ascending=False)
        .reset_index(drop=True)
    )

    extremes = {}
    for model_name in res_summary["model"]:
        sub = res_per_gene[res_per_gene["model"] == model_name].copy()
        i_min = sub["R2_spatialBlockCV"].idxmin()
        i_max = sub["R2_spatialBlockCV"].idxmax()

        extreme_rows = [sub.loc[i_min]]
        if i_max != i_min:
            extreme_rows.append(sub.loc[i_max])

        extremes[model_name] = extreme_rows

    fig, ax = plt.subplots(figsize=figsize)

    x = np.arange(len(res_summary))
    means = res_summary["mean_R2_across_genes"].to_numpy()
    yerr = np.vstack(
        [
            means - res_summary["ci95_low"].to_numpy(),
            res_summary["ci95_high"].to_numpy() - means,
        ]
    )

    ax.bar(x, means, alpha=bar_alpha)

    ax.errorbar(
        x,
        means,
        yerr=yerr,
        fmt="none",
        ecolor="black",
        elinewidth=2,
        capsize=capsize,
        capthick=2,
        zorder=4,
    )

    for xi, model_name in enumerate(res_summary["model"]):
        for row in extremes[model_name]:
            gene = row["gene"]
            r2 = float(row["R2_spatialBlockCV"])
            r2_disp = float(ceil_to_decimals(r2, 2))

            ax.scatter([xi], [r2], color=point_color, s=60, zorder=5)

            if annotate:
                ax.text(
                    xi + text_dx,
                    r2 + text_dy,
                    f"{gene} ({r2_disp:.2f})",
                    fontsize=annotate_fontsize,
                    ha="left",
                    va="bottom",
                    color=point_color,
                    zorder=6,
                )

    ax.set_xticks(x)
    ax.set_xticklabels(res_summary["model"], rotation=25, ha="right")
    ax.set_ylabel(
        f"OOF R² (spatial-block CV), mean across {int(res_summary['n_genes'].iloc[0])} genes"
    )
    ax.set_title(
        f"Model comparison across genes (OOF R² mean ± {int(confidence * 100)}% CI) + min/max gene points"
    )

    ax.grid(True, axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show()

    return res_summary


class DistSigInteractionOnly(BaseEstimator, TransformerMixin):
    """
    Assumes X = [dist_spline_features..., signature] where signature is the last column.
    Returns [dist_spline_features..., signature, dist_spline_features * signature].
    """

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X = np.asarray(X)
        dist_basis = X[:, :-1]
        sig = X[:, [-1]]
        inter = dist_basis * sig
        return np.hstack([dist_basis, sig, inter])


def make_dist_signature_model(
    dist_col, sig_col, n_knots=6, degree=3, alpha=1.0, knots="quantile"
):
    pre = ColumnTransformer(
        transformers=[
            (
                "dist",
                SplineTransformer(
                    n_knots=n_knots, degree=degree, knots=knots, include_bias=False
                ),
                [dist_col],
            ),
            ("sig", "passthrough", [sig_col]),
        ],
        remainder="drop",
        sparse_threshold=0.0,
    )
    return Pipeline(
        [
            ("pre", pre),
            ("scale", StandardScaler()),
            ("inter", DistSigInteractionOnly()),
            ("ridge", Ridge(alpha=alpha)),
        ]
    )


def make_spline_only_model(dist_col, n_knots=6, degree=3, alpha=1.0, knots="quantile"):
    pre = ColumnTransformer(
        transformers=[
            (
                "dist",
                SplineTransformer(
                    n_knots=n_knots, degree=degree, knots=knots, include_bias=False
                ),
                [dist_col],
            ),
        ],
        remainder="drop",
        sparse_threshold=0.0,
    )
    return Pipeline(
        [
            ("pre", pre),
            ("scale", StandardScaler()),
            ("ridge", Ridge(alpha=alpha)),
        ]
    )


def make_signature_only_model(sig_col, alpha=1.0):
    pre = ColumnTransformer(
        transformers=[
            ("sig", "passthrough", [sig_col]),
        ],
        remainder="drop",
        sparse_threshold=0.0,
    )
    return Pipeline(
        [
            ("pre", pre),
            ("scale", StandardScaler()),
            ("ridge", Ridge(alpha=alpha)),
        ]
    )


def cv_r2_spatial_blocks_estimator(
    df, target_col, feature_cols, groups, estimator, n_splits=5
):
    """
    Out-of-fold R^2 computed over all samples using per-sample OOF predictions.
    Returns:
      oof_r2: pooled OOF R^2
      fold_scores: per-fold R^2 values (diagnostic)
    """
    X = df[feature_cols]
    y = df[target_col].to_numpy()

    gkf = GroupKFold(n_splits=n_splits)

    y_oof = np.full(shape=y.shape, fill_value=np.nan, dtype=float)
    fold_scores = []

    for tr_idx, te_idx in gkf.split(X, y, groups=groups):
        est = clone(estimator)
        est.fit(X.iloc[tr_idx], y[tr_idx])

        yhat = est.predict(X.iloc[te_idx])
        yhat = np.asarray(yhat).reshape(-1)

        y_oof[te_idx] = yhat
        fold_scores.append(r2_score(y[te_idx], yhat))

    if np.isnan(y_oof).any():
        raise RuntimeError(
            "OOF prediction vector contains NaNs. This indicates some samples were never assigned to a test fold."
        )

    oof_r2 = r2_score(y, y_oof)
    return float(oof_r2), fold_scores


def plot_dist_signature_heatmaps_grid(
    df,
    dist_col,
    sig_col,
    y_cols,
    q=5,
    y_is_log1p=True,
    aggfunc="mean",
    nrows=4,
    ncols=4,
    figsize=(18, 14),
    cmap="viridis",
    color_scale="log1p",
    norm_scope="per_gene",
    share_colorbar=False,
    eps=1e-6,
    symlog_linthresh=0.1,
    annotate=True,
    annot_decimals=2,
    annot_fontsize=7,
    label_mode="index",
):
    if len(y_cols) != nrows * ncols:
        raise ValueError(f"Expected {nrows * ncols} genes, got {len(y_cols)}")

    if share_colorbar:
        norm_scope = "global"

    base = df[[dist_col, sig_col]].dropna()
    _, dist_edges = pd.qcut(base[dist_col], q=q, duplicates="drop", retbins=True)
    _, sig_edges = pd.qcut(base[sig_col], q=q, duplicates="drop", retbins=True)

    dist_cats = pd.cut(
        base[dist_col], bins=dist_edges, include_lowest=True
    ).cat.categories
    sig_cats = pd.cut(base[sig_col], bins=sig_edges, include_lowest=True).cat.categories

    if label_mode == "interval":
        xlabels = [str(c) for c in dist_cats]
        ylabels = [str(c) for c in sig_cats]
    elif label_mode == "index":
        xlabels = [f"D{i + 1}" for i in range(len(dist_cats))]
        ylabels = [f"S{i + 1}" for i in range(len(sig_cats))]
    else:
        raise ValueError("label_mode must be 'index' or 'interval'")

    mean_tables = []
    colored_arrays = []

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

        A = mean_tbl.to_numpy()

        if color_scale == "linear":
            A_col = A
        elif color_scale == "log1p":
            A_col = np.log1p(A)
        elif color_scale == "log10":
            A_col = np.log10(A + eps)
        elif color_scale == "lognorm":
            A_col = A + eps
        elif color_scale == "symlognorm":
            A_col = A
        else:
            raise ValueError(
                "color_scale must be one of: linear, log1p, log10, lognorm, symlognorm"
            )

        colored_arrays.append(A_col)

    global_norm = None
    if norm_scope == "global":
        stacked = np.array([np.ravel(a) for a in colored_arrays], dtype=float)
        vmin = np.nanmin(stacked)
        vmax = np.nanmax(stacked)

        if color_scale == "lognorm":
            vmin = max(float(vmin), eps)
            global_norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        elif color_scale == "symlognorm":
            global_norm = colors.SymLogNorm(
                linthresh=symlog_linthresh, vmin=vmin, vmax=vmax
            )
        else:
            global_norm = colors.Normalize(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    axes = axes.ravel()

    last_im = None

    for ax, gene, mean_tbl, A_col in zip(axes, y_cols, mean_tables, colored_arrays):
        norm = global_norm
        if norm_scope == "per_gene":
            vmin = np.nanmin(A_col)
            vmax = np.nanmax(A_col)
            if color_scale == "lognorm":
                vmin = max(float(vmin), eps)
                norm = colors.LogNorm(vmin=vmin, vmax=vmax)
            elif color_scale == "symlognorm":
                norm = colors.SymLogNorm(
                    linthresh=symlog_linthresh, vmin=vmin, vmax=vmax
                )
            else:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)

        im = ax.imshow(A_col, aspect="auto", cmap=cmap, norm=norm)
        last_im = im

        ax.set_title(gene)
        ax.set_xticks(range(len(dist_cats)))
        ax.set_yticks(range(len(sig_cats)))
        ax.set_xticklabels(xlabels, rotation=45, ha="right")
        ax.set_yticklabels(ylabels)

        ax.set_xticks(np.arange(-0.5, len(dist_cats), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(sig_cats), 1), minor=True)
        ax.grid(which="minor", linestyle="--", linewidth=0.5, alpha=0.25)
        ax.tick_params(which="minor", bottom=False, left=False)

        if annotate:
            V = mean_tbl.to_numpy()
            for i in range(V.shape[0]):
                for j in range(V.shape[1]):
                    val = V[i, j]
                    if np.isnan(val):
                        continue
                    ax.text(
                        j,
                        i,
                        f"{val:.{annot_decimals}f}",
                        ha="center",
                        va="center",
                        fontsize=annot_fontsize,
                    )

    fig.supxlabel("Distance bin")
    fig.supylabel("Signature bin")

    if share_colorbar and last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes, shrink=0.8)
        if color_scale == "linear":
            cbar.set_label(f"{aggfunc}(mean count)")
        elif color_scale == "log1p":
            cbar.set_label(f"log1p({aggfunc}(mean count))")
        elif color_scale == "log10":
            cbar.set_label(f"log10({aggfunc}(mean count)+eps)")
        elif color_scale == "lognorm":
            cbar.set_label(f"{aggfunc}(mean count) [LogNorm]")
        elif color_scale == "symlognorm":
            cbar.set_label(f"{aggfunc}(mean count) [SymLogNorm]")

    fig.suptitle(
        f"Interaction maps (color_scale={color_scale}, norm_scope={norm_scope}, q={q})",
        y=1.02,
    )
    plt.show()

    return dist_edges, sig_edges


def make_linear_ridge_model(feature_cols, alpha=1.0):
    return Pipeline(
        [
            ("scale", StandardScaler(with_mean=True, with_std=True)),
            ("ridge", Ridge(alpha=alpha)),
        ]
    )


def make_binned_distance_model(dist_col, n_bins=5, alpha=1.0, strategy="quantile"):
    pre = ColumnTransformer(
        transformers=[
            (
                "distb",
                KBinsDiscretizer(n_bins=n_bins, encode="onehot", strategy=strategy),
                [dist_col],
            ),
        ],
        remainder="drop",
    )
    return Pipeline(
        [
            ("pre", pre),
            ("ridge", Ridge(alpha=alpha)),
        ]
    )


def add_distance_transforms(df, dist_col, eps=1e-6):
    out = df.copy()
    d = out[dist_col].to_numpy()
    out["dist_log"] = np.log(d + eps)
    out["dist_sq"] = d**2
    return out


def plot_binned_distance_effect(df, dist_col, y_col, q=30):
    tmp = df[[dist_col, y_col]].dropna().copy()
    tmp["dist_bin"] = pd.qcut(tmp[dist_col], q=q, duplicates="drop")
    g = tmp.groupby("dist_bin", observed=True)

    x = g[dist_col].mean()
    y = g[y_col].mean()
    n = g.size()
    se = g[y_col].std() / np.sqrt(n)

    plt.figure()
    plt.plot(x, y, marker="o")
    plt.fill_between(x, y - 1.96 * se, y + 1.96 * se, alpha=0.2)
    plt.xlabel(dist_col)
    plt.ylabel(f"Mean({y_col})")
    plt.title(f"{y_col} vs {dist_col} (binned mean ± 95% CI)")
    plt.tight_layout()
    plt.show()


def plot_dist_signature_heatmaps(df, dist_col, sig_col, y_col, q=5):
    tmp = df[[dist_col, sig_col, y_col]].dropna().copy()
    tmp["dist_bin"] = pd.qcut(tmp[dist_col], q=q, duplicates="drop")
    tmp["sig_bin"] = pd.qcut(tmp[sig_col], q=q, duplicates="drop")

    mean_tbl = tmp.pivot_table(
        index="sig_bin", columns="dist_bin", values=y_col, aggfunc="mean", observed=True
    )
    n_tbl = tmp.pivot_table(
        index="sig_bin", columns="dist_bin", values=y_col, aggfunc="size", observed=True
    )

    plt.figure()
    plt.imshow(mean_tbl.values, aspect="auto")
    plt.colorbar(label=f"Mean({y_col})")
    plt.xticks(
        range(mean_tbl.shape[1]),
        [str(c) for c in mean_tbl.columns],
        rotation=45,
        ha="right",
    )
    plt.yticks(range(mean_tbl.shape[0]), [str(r) for r in mean_tbl.index])
    plt.xlabel("Distance bin")
    plt.ylabel("Signature bin")
    plt.title(f"Interaction map: Mean({y_col})")
    plt.tight_layout()
    plt.show()


def make_tile_groups(df, x_col, y_col, n_tiles_x=10, n_tiles_y=10):
    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()

    x_edges = np.linspace(np.nanmin(x), np.nanmax(x), n_tiles_x + 1)
    y_edges = np.linspace(np.nanmin(y), np.nanmax(y), n_tiles_y + 1)

    x_bin = np.clip(np.digitize(x, x_edges) - 1, 0, n_tiles_x - 1)
    y_bin = np.clip(np.digitize(y, y_edges) - 1, 0, n_tiles_y - 1)

    tile_id = x_bin + n_tiles_x * y_bin
    return tile_id.astype(int)


def _fit_predict_ridge(X_train, y_train, X_test, alpha=1.0):
    model = Pipeline(
        [
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("ridge", Ridge(alpha=alpha, random_state=0)),
        ]
    )
    model.fit(X_train, y_train)
    return model.predict(X_test)


def cv_r2_random(df, target, feature_cols, n_splits=5, alpha=1.0, seed=0):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    y = df[target].to_numpy()
    X = df[feature_cols].to_numpy()

    r2s = []
    for tr, te in kf.split(X):
        pred = _fit_predict_ridge(X[tr], y[tr], X[te], alpha=alpha)
        r2s.append(r2_score(y[te], pred))
    return float(np.mean(r2s)), r2s


def cv_r2_spatial_blocks(df, target, feature_cols, groups, n_splits=5, alpha=1.0):
    gkf = GroupKFold(n_splits=n_splits)
    y = df[target].to_numpy()
    X = df[feature_cols].to_numpy()

    r2s = []
    for tr, te in gkf.split(X, y, groups=groups):
        pred = _fit_predict_ridge(X[tr], y[tr], X[te], alpha=alpha)
        r2s.append(r2_score(y[te], pred))
    return float(np.mean(r2s)), r2s


def evaluate_leakage_for_gene(
    df,
    target_gene,
    n_tiles_x=10,
    n_tiles_y=10,
    n_splits=5,
    eps=1e-8,
    x_col='x_centroid', 
    y_col='y_centroid',
    dist_col='distance_to_plaque',
    optional_covs=[],
    neigh_cols=None,
):
    groups = make_tile_groups(
        df, x_col, y_col, n_tiles_x=n_tiles_x, n_tiles_y=n_tiles_y
    )
    current_neighbors = [c for c in neigh_cols if target_gene not in c]

    FEAT_XY = [x_col, y_col]
    FEAT_DIST = [dist_col]
    FEAT_NEIGH = current_neighbors
    FEAT_FULL = [dist_col] + optional_covs + current_neighbors

    results = {}
    for name, feats in [
        ("XY_only", FEAT_XY),
        ("Distance_only", FEAT_DIST),
        ("Neighbor_only", FEAT_NEIGH),
        ("Full_(dist+neighbors+XY)", FEAT_FULL),
    ]:
        r2_rand, _ = cv_r2_random(
            df, target_gene, feats, n_splits=n_splits, alpha=1.0, seed=0
        )
        r2_spat, _ = cv_r2_spatial_blocks(
            df, target_gene, feats, groups, n_splits=n_splits, alpha=1.0
        )
        results[name] = {"R2_randomCV": r2_rand, "R2_spatialBlockCV": r2_spat}

    out = pd.DataFrame(results).T
    out["Leakage_gap_(random - spatial)"] = (
        out["R2_randomCV"] - out["R2_spatialBlockCV"]
    )

    denom = out["R2_randomCV"].where(out["R2_randomCV"] > eps, np.nan)
    out["Relative leakage gap, %"] = (
        out["Leakage_gap_(random - spatial)"] / denom
    ) * 100
    out["Relative leakage gap, %"] = out["Relative leakage gap, %"].round(2)

    return out.sort_values("R2_spatialBlockCV", ascending=False)


SCATTER_COLOR = "maroon"
SCATTER_ZORDER = 6


def sem(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) <= 1:
        return np.nan
    return float(x.std(ddof=1) / np.sqrt(len(x)))


def _apply_figure_margins(fig, *, bottom=0.30, top=0.90, left=0.08, right=0.98):
    fig.subplots_adjust(bottom=bottom, top=top, left=left, right=right)


def _set_ylim_with_padding_linear(ax, y_vals, *, pad_frac=0.12):
    y_vals = np.asarray(
        [v for v in y_vals if v is not None and np.isfinite(v)], dtype=float
    )
    if y_vals.size == 0:
        return
    y_min = np.min(y_vals)
    y_max = np.max(y_vals)
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0
    pad = pad_frac * (y_max - y_min)
    ax.set_ylim(y_min - pad, y_max + pad)


def _set_ylim_with_padding_symlog(ax, y_vals, *, upper_mult=1.8, lower_mult=1.8):
    y_vals = np.asarray(
        [v for v in y_vals if v is not None and np.isfinite(v)], dtype=float
    )
    if y_vals.size == 0:
        return

    y_min = float(np.min(y_vals))
    y_max = float(np.max(y_vals))

    if y_max > 0:
        y_upper = y_max * upper_mult
    else:
        y_upper = y_max / upper_mult

    if y_min < 0:
        y_lower = y_min * lower_mult
    elif y_min > 0:
        y_lower = max(y_min / lower_mult, 1e-6)
    else:
        y_lower = -1.0

    ax.set_ylim(y_lower, y_upper)


def _annotate_dot(ax, x, y, text, *, y_upper_hint=None, fontsize=8):
    if (
        y_upper_hint is not None
        and np.isfinite(y_upper_hint)
        and y > 0.85 * y_upper_hint
    ):
        dy = -12
        va = "top"
    else:
        dy = 8
        va = "bottom"

    ax.annotate(
        text,
        xy=(x, y),
        xytext=(0, dy),
        textcoords="offset points",
        ha="center",
        va=va,
        fontsize=fontsize,
        clip_on=False,
    )


def build_leakage_all(
    df: pd.DataFrame,
    target_genes: list[str],
    model_order: list[str],
    *,
    n_tiles_x: int = 10,
    n_tiles_y: int = 10,
    n_splits: int = 5,
    x_col = None,
    y_col = None,
    neigh_cols = None,
) -> pd.DataFrame:
    missing = [g for g in target_genes if g not in df.columns]
    if missing:
        raise ValueError(f"Missing target genes in df: {missing}")

    rows = []
    for g in target_genes:
        tab = evaluate_leakage_for_gene(
            df, g, n_tiles_x=n_tiles_x, n_tiles_y=n_tiles_y, n_splits=n_splits, x_col=x_col, y_col=y_col, neigh_cols=neigh_cols
        ).copy()
        tab.insert(0, "target_gene", g)
        tab.insert(1, "model", tab.index)
        rows.append(tab.reset_index(drop=True))

    leakage_all = pd.concat(rows, ignore_index=True)

    expected = set(model_order)
    found = set(leakage_all["model"].unique())
    print(found)
    missing_models = expected - found
    if missing_models:
        raise ValueError(
            f"Some expected models are missing from leakage_all: {sorted(missing_models)}"
        )

    return leakage_all


def aggregate_leakage(
    leakage_all: pd.DataFrame, model_order: list[str], *, rel_gap_col: str
) -> pd.DataFrame:
    required = {"model", "target_gene", "R2_randomCV", "R2_spatialBlockCV", rel_gap_col}
    missing = required - set(leakage_all.columns)
    if missing:
        raise ValueError(f"leakage_all missing required columns: {sorted(missing)}")

    def _sem(x) -> float:
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        if len(x) <= 1:
            return np.nan
        return float(x.std(ddof=1) / np.sqrt(len(x)))

    def _bootstrap_rel_gap_sem(
        rand_vals, spat_vals, *, n_boot=2000, seed=0, eps=1e-8
    ) -> float:
        rand_vals = np.asarray(rand_vals, dtype=float)
        spat_vals = np.asarray(spat_vals, dtype=float)

        mask = np.isfinite(rand_vals) & np.isfinite(spat_vals)
        rand_vals = rand_vals[mask]
        spat_vals = spat_vals[mask]

        n = len(rand_vals)
        if n <= 1:
            return np.nan

        rng = np.random.default_rng(seed)
        idx = rng.integers(0, n, size=(n_boot, n))

        r_mean = rand_vals[idx].mean(axis=1)
        s_mean = spat_vals[idx].mean(axis=1)

        denom = np.where(r_mean > eps, r_mean, np.nan)
        g = 100.0 * (r_mean - s_mean) / denom

        return float(np.nanstd(g, ddof=1))

    rows = []
    for m, sub in leakage_all.groupby("model", as_index=False):
        rand_vals = sub["R2_randomCV"].to_numpy(dtype=float)
        spat_vals = sub["R2_spatialBlockCV"].to_numpy(dtype=float)

        r_mean = float(np.nanmean(rand_vals))
        s_mean = float(np.nanmean(spat_vals))

        r_sem = _sem(rand_vals)
        s_sem = _sem(spat_vals)

        eps = 1e-8
        rel_gap_mean = np.nan
        if np.isfinite(r_mean) and r_mean > eps and np.isfinite(s_mean):
            rel_gap_mean = 100.0 * (r_mean - s_mean) / r_mean

        rel_gap_sem = _bootstrap_rel_gap_sem(
            rand_vals, spat_vals, n_boot=2000, seed=0, eps=eps
        )

        rel_gap_by_gene_mean = float(np.nanmean(sub[rel_gap_col].to_numpy(dtype=float)))
        rel_gap_by_gene_sem = _sem(sub[rel_gap_col].to_numpy(dtype=float))

        rows.append(
            {
                "model": m,
                "n_genes": sub["target_gene"].nunique(),
                "R2_randomCV_mean": r_mean,
                "R2_randomCV_sem": r_sem,
                "R2_spatialBlockCV_mean": s_mean,
                "R2_spatialBlockCV_sem": s_sem,
                "rel_gap_pct_mean": rel_gap_mean,
                "rel_gap_pct_sem": rel_gap_sem,
                "rel_gap_pct_mean_of_gene_ratios": rel_gap_by_gene_mean,
                "rel_gap_pct_sem_of_gene_ratios": rel_gap_by_gene_sem,
            }
        )

    agg = pd.DataFrame(rows).set_index("model")
    return agg.loc[model_order]


def extremes_by_model(
    leakage_all: pd.DataFrame, model_order: list[str], value_col: str
) -> pd.DataFrame:
    required = {"model", "target_gene", value_col}
    missing = required - set(leakage_all.columns)
    if missing:
        raise ValueError(f"leakage_all missing required columns: {sorted(missing)}")

    out_rows = []
    for m in model_order:
        sub = leakage_all[leakage_all["model"] == m][
            ["target_gene", value_col]
        ].dropna()
        if sub.empty:
            out_rows.append((m, None, np.nan, None, np.nan))
            continue
        i_min = sub[value_col].idxmin()
        i_max = sub[value_col].idxmax()
        out_rows.append(
            (
                m,
                leakage_all.loc[i_min, "target_gene"],
                float(leakage_all.loc[i_min, value_col]),
                leakage_all.loc[i_max, "target_gene"],
                float(leakage_all.loc[i_max, value_col]),
            )
        )

    return pd.DataFrame(
        out_rows, columns=["model", "min_gene", "min_value", "max_gene", "max_value"]
    ).set_index("model")


def plot_mean_r2_with_extremes(
    agg: pd.DataFrame,
    leakage_all: pd.DataFrame,
    model_order: list[str],
    *,
    title: str = "Mean variance explained, averaged over PIGs",
    annotate_extremes: bool = True,
    fontsize: int = 8,
    figsize=(12, 5),
):
    models = model_order
    x = np.arange(len(models))
    w = 0.38

    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(
        x - w / 2,
        agg["R2_randomCV_mean"].to_numpy(),
        yerr=agg["R2_randomCV_sem"].to_numpy(),
        width=w,
        capsize=3,
        label="Random cross-validation",
    )
    ax.bar(
        x + w / 2,
        agg["R2_spatialBlockCV_mean"].to_numpy(),
        yerr=agg["R2_spatialBlockCV_sem"].to_numpy(),
        width=w,
        capsize=3,
        label="Spatial block cross-validation",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=25, ha="right")
    ax.set_ylabel("Mean OOF $R^2$ across genes")
    ax.set_title(title, pad=12)
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.25)
    ax.margins(x=0.08)

    y_for_limits = []
    y_for_limits += list((agg["R2_randomCV_mean"] - agg["R2_randomCV_sem"]).to_numpy())
    y_for_limits += list((agg["R2_randomCV_mean"] + agg["R2_randomCV_sem"]).to_numpy())
    y_for_limits += list(
        (agg["R2_spatialBlockCV_mean"] - agg["R2_spatialBlockCV_sem"]).to_numpy()
    )
    y_for_limits += list(
        (agg["R2_spatialBlockCV_mean"] + agg["R2_spatialBlockCV_sem"]).to_numpy()
    )

    if annotate_extremes:
        ex_rand = extremes_by_model(leakage_all, model_order, "R2_randomCV")
        ex_spat = extremes_by_model(leakage_all, model_order, "R2_spatialBlockCV")

        y_for_limits += list(ex_rand["min_value"].to_numpy())
        y_for_limits += list(ex_rand["max_value"].to_numpy())
        y_for_limits += list(ex_spat["min_value"].to_numpy())
        y_for_limits += list(ex_spat["max_value"].to_numpy())

        _set_ylim_with_padding_linear(ax, y_for_limits, pad_frac=0.16)
        y_upper_hint = ax.get_ylim()[1]

        for i, m in enumerate(models):
            xr = x[i] - w / 2
            xs = x[i] + w / 2

            if pd.notna(ex_rand.loc[m, "min_value"]):
                v = ex_rand.loc[m, "min_value"]
                ax.scatter([xr], [v], color=SCATTER_COLOR, zorder=SCATTER_ZORDER)
                _annotate_dot(
                    ax,
                    xr,
                    v,
                    f"{ex_rand.loc[m,'min_gene']}, {v:.4f}",
                    y_upper_hint=y_upper_hint,
                    fontsize=fontsize,
                )

            if pd.notna(ex_rand.loc[m, "max_value"]):
                v = ex_rand.loc[m, "max_value"]
                ax.scatter([xr], [v], color=SCATTER_COLOR, zorder=SCATTER_ZORDER)
                _annotate_dot(
                    ax,
                    xr,
                    v,
                    f"{ex_rand.loc[m,'max_gene']}, {v:.4f}",
                    y_upper_hint=y_upper_hint,
                    fontsize=fontsize,
                )

            if pd.notna(ex_spat.loc[m, "min_value"]):
                v = ex_spat.loc[m, "min_value"]
                ax.scatter([xs], [v], color=SCATTER_COLOR, zorder=SCATTER_ZORDER)
                _annotate_dot(
                    ax,
                    xs,
                    v,
                    f"{ex_spat.loc[m,'min_gene']}, {v:.4f}",
                    y_upper_hint=y_upper_hint,
                    fontsize=fontsize,
                )

            if pd.notna(ex_spat.loc[m, "max_value"]):
                v = ex_spat.loc[m, "max_value"]
                ax.scatter([xs], [v], color=SCATTER_COLOR, zorder=SCATTER_ZORDER)
                _annotate_dot(
                    ax,
                    xs,
                    v,
                    f"{ex_spat.loc[m,'max_gene']}, {v:.4f}",
                    y_upper_hint=y_upper_hint,
                    fontsize=fontsize,
                )
    else:
        _set_ylim_with_padding_linear(ax, y_for_limits, pad_frac=0.16)

    _apply_figure_margins(fig, bottom=0.32, top=0.88)
    plt.show()


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
):
    models = model_order
    x = np.arange(len(models))

    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(
        x,
        agg["rel_gap_pct_mean"].to_numpy(),
        yerr=agg["rel_gap_pct_sem"].to_numpy(),
        capsize=3,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=25, ha="right")
    ax.set_ylabel("Mean relative leakage gap (%)")
    ax.set_title(title, pad=12)
    ax.grid(True, axis="y", alpha=0.25)
    ax.margins(x=0.08)

    if yscale == "log":
        if np.any(agg["rel_gap_pct_mean"].to_numpy() <= 0):
            raise ValueError(
                "yscale='log' requires strictly positive mean relative gaps. Use yscale='symlog' instead."
            )
        ax.set_yscale("log")
    elif yscale == "symlog":
        ax.set_yscale("symlog", linthresh=symlog_linthresh)
    else:
        raise ValueError("yscale must be one of {'log','symlog'}.")

    y_for_limits = []
    y_for_limits += list((agg["rel_gap_pct_mean"] - agg["rel_gap_pct_sem"]).to_numpy())
    y_for_limits += list((agg["rel_gap_pct_mean"] + agg["rel_gap_pct_sem"]).to_numpy())

    if annotate_extremes:
        ex_gap = extremes_by_model(leakage_all, model_order, rel_gap_col)
        y_for_limits += list(ex_gap["min_value"].to_numpy())
        y_for_limits += list(ex_gap["max_value"].to_numpy())

        _set_ylim_with_padding_symlog(ax, y_for_limits, upper_mult=2.2, lower_mult=2.2)
        y_upper_hint = ax.get_ylim()[1]

        for i, m in enumerate(models):
            xm = x[i]

            if pd.notna(ex_gap.loc[m, "min_value"]):
                v = ex_gap.loc[m, "min_value"]
                ax.scatter([xm], [v], color=SCATTER_COLOR, zorder=SCATTER_ZORDER)
                _annotate_dot(
                    ax,
                    xm,
                    v,
                    f"{ex_gap.loc[m,'min_gene']}, {v:.1f}",
                    y_upper_hint=y_upper_hint,
                    fontsize=fontsize,
                )

            if pd.notna(ex_gap.loc[m, "max_value"]):
                v = ex_gap.loc[m, "max_value"]
                ax.scatter([xm], [v], color=SCATTER_COLOR, zorder=SCATTER_ZORDER)
                _annotate_dot(
                    ax,
                    xm,
                    v,
                    f"{ex_gap.loc[m,'max_gene']}, {v:.1f}",
                    y_upper_hint=y_upper_hint,
                    fontsize=fontsize,
                )
    else:
        _set_ylim_with_padding_symlog(ax, y_for_limits, upper_mult=2.2, lower_mult=2.2)

    _apply_figure_margins(fig, bottom=0.32, top=0.88)
    plt.show()


def plot_true_vs_pred_heatmaps(
    true_grid,
    pred_grid,
    *,
    suptitle="200x200 grid",
    titles=("True distance", "Decision tree inferred distance"),
    cmap="magma",
    vmin=0.0,
    vmax=None,
    vmax_percentile=99.0,
    gamma=0.5,
    cbar_label="distance to plaque",
):
    """
    Side-by-side heatmaps with a unified (shared) color scale and a single colorbar.

    Uses PowerNorm (gamma < 1) to expand low-end contrast, and a robust vmax based on a percentile.
    """
    true_grid = np.asarray(true_grid, dtype=float)
    pred_grid = np.asarray(pred_grid, dtype=float)

    combined = np.concatenate([true_grid.ravel(), pred_grid.ravel()])
    combined = combined[~np.isnan(combined)]
    if combined.size == 0:
        raise ValueError("Both grids are all-NaN; cannot plot.")

    if vmax is None:
        vmax = np.percentile(combined, vmax_percentile)

    norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)

    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(color="white")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    fig.suptitle(suptitle)

    im0 = axes[0].imshow(true_grid, origin="lower", aspect="auto", cmap=cm, norm=norm)
    axes[0].set_title(titles[0])

    im1 = axes[1].imshow(pred_grid, origin="lower", aspect="auto", cmap=cm, norm=norm)
    axes[1].set_title(titles[1])

    cbar = fig.colorbar(im1, ax=axes, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label)

    plt.show()


def plot_spatial_block_split(
    work: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    grid_size: float,
    test_blocks: set,
    x0: float | None = None,
    y0: float | None = None,
    train_mask: pd.Series | None = None,
    test_mask: pd.Series | None = None,
    figsize: tuple[int, int] = (10, 8),
    title: str = "Spatial Dataset Split",
    show: bool = True,
):
    if x0 is None:
        x0 = float(work[x_col].min())
    if y0 is None:
        y0 = float(work[y_col].min())

    if not {"bx", "by", "block_id"}.issubset(work.columns):
        x = work[x_col].to_numpy()
        y = work[y_col].to_numpy()
        bx = np.floor((x - x0) / grid_size).astype(int)
        by = np.floor((y - y0) / grid_size).astype(int)

        work = work.copy()
        work["bx"] = bx
        work["by"] = by

        block_id = bx * (by.max() + 1) + by
        work["block_id"] = block_id

    if test_mask is None:
        test_mask = work["block_id"].isin(test_blocks)
    if train_mask is None:
        train_mask = ~work["block_id"].isin(test_blocks)

    block_meta = work.drop_duplicates("block_id")[["block_id", "bx", "by"]].copy()
    block_meta["is_test"] = block_meta["block_id"].isin(test_blocks)

    fig, ax = plt.subplots(figsize=figsize)

    for r in block_meta.itertuples(index=False):
        tx0 = x0 + r.bx * grid_size
        ty0 = y0 + r.by * grid_size

        if r.is_test:
            rect = Rectangle((tx0, ty0), grid_size, grid_size, fill=True, alpha=0.20)
        else:
            rect = Rectangle(
                (tx0, ty0), grid_size, grid_size, fill=False, linewidth=0.8
            )

        ax.add_patch(rect)

    ax.scatter(
        work.loc[train_mask, x_col],
        work.loc[train_mask, y_col],
        s=8,
        alpha=0.35,
        label="Train cells",
    )

    ax.scatter(
        work.loc[test_mask, x_col],
        work.loc[test_mask, y_col],
        s=10,
        alpha=0.70,
        label="Test cells",
    )

    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_aspect("equal", adjustable="box")
    ax.legend()

    if show:
        plt.show()

    return fig, ax


def find_xy_columns(df):
    candidates = [
        ("x_centroid", "y_centroid"),
        ("x", "y"),
        ("X", "Y"),
        ("centroid_x", "centroid_y"),
    ]
    for x, y in candidates:
        if x in df.columns and y in df.columns:
            return x, y

    x_like = [c for c in df.columns if "x" in c.lower()]
    y_like = [c for c in df.columns if "y" in c.lower()]
    for x in x_like:
        for y in y_like:
            if ("centroid" in x.lower()) or ("centroid" in y.lower()):
                return x, y

    raise ValueError(
        "Couldn't find x/y columns. "
        "Tried common names like x_centroid/y_centroid, x/y, X/Y, centroid_x/centroid_y. "
        f"Here are some likely candidates:\n"
        f"x-like: {x_like[:10]}\n"
        f"y-like: {y_like[:10]}"
    )
