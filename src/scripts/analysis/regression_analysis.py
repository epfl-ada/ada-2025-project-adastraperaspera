"""
Gene-wise regression and correlation analyses versus plaque distance.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from src.utils.logging_utils import logger
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR adjustment; returns adjusted p-values in original order."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)

    mask = np.isfinite(pvals)
    pv = pvals[mask]
    m = pv.size
    if m == 0:
        return out

    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m / (np.arange(1, m + 1))

    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    inv_order = np.empty_like(order)
    inv_order[order] = np.arange(m)
    out[mask] = q[inv_order]
    return out


def plot_nested_regression_adj_r2(
    nested_regression_results: pd.DataFrame,
    gene_col: str = "gene",
    model_col: str = "model",
    y_col: str = "adj_r_squared",
    p_col: str = "f_test_pval",
    alpha: float = 0.01,
    adjust_pvals: bool = True,
    model_order: list[str] | None = None,
    mode: str = "genes",
    print_summary: bool = True,
    annotate_minmax_genes: bool = True,
    annotate_fontsize: int = 8,
    add_model3_legend: bool = True,
    model3_legend_title: str = "Model 3(k) variants",
    model3_legend_text: str | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (12, 6),
    savepath: str | None = None,
    show: bool = True,
):
    """
    Visualize adjusted R² across nested regression models.

    mode="genes":
      - one line per gene
      - markers per point for models > 0:
          black square if (adjusted) p >= alpha
          red star if (adjusted) p < alpha
      - gene legend on the right

    mode="average":
      - one line showing mean adjusted R² across genes per model
      - error bars span min..max (across genes) per model
      - optional min/max gene labels at the error bar endpoints
      - optional Model 3(k) explanatory legend on the right

    Returns
    -------
    (fig, ax, summary_df)
      summary_df is only populated in mode="average"; otherwise None.
    """
    df = nested_regression_results.copy()

    required = {gene_col, model_col, y_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if p_col not in df.columns:
        raise ValueError(
            f"Missing p-value column '{p_col}'. Either add it or change p_col."
        )

    if model_order is None:
        models = df[model_col].astype(str).unique().tolist()
        base = [m for m in ["0", "1", "2"] if m in models]

        def parse_k(m: str):
            if isinstance(m, str) and m.startswith("3_"):
                try:
                    return int(m.split("_", 1)[1])
                except Exception:
                    return None
            return None

        m3 = [m for m in models if isinstance(m, str) and m.startswith("3_")]
        m3_sorted = sorted(
            m3, key=lambda m: (parse_k(m) is None, parse_k(m) or 10**9, m)
        )

        known = set(base) | set(m3_sorted)
        other = sorted([m for m in models if m not in known])
        model_order = base + m3_sorted + other

    def model_label(m: str) -> str:
        m = str(m)
        if m in {"0", "1", "2"}:
            return f"Model {m}"
        if m.startswith("3_"):
            try:
                k = int(m.split("_", 1)[1])
                return f"Model 3({k})"
            except Exception:
                return f"Model {m}"
        return f"Model {m}"

    m3_ks = []
    for m in model_order:
        sm = str(m)
        if sm.startswith("3_"):
            try:
                m3_ks.append(int(sm.split("_", 1)[1]))
            except Exception:
                pass
    has_model3_variants = len(m3_ks) > 0

    df["_p_for_sig_"] = pd.to_numeric(df[p_col], errors="coerce").astype(float)

    if adjust_pvals:
        df["_p_adj_"] = np.nan
        for m, g in df.groupby(model_col, sort=False):
            if str(m) == "0":
                continue
            df.loc[g.index, "_p_adj_"] = _bh_adjust(g["_p_for_sig_"].to_numpy())
        p_used = "_p_adj_"
        p_label = "BH-adjusted p"
    else:
        p_used = "_p_for_sig_"
        p_label = "p"

    x_labels = [model_label(m) for m in model_order]
    x = np.arange(len(model_order))

    fig, ax = plt.subplots(figsize=figsize)
    summary_df = None

    mode = mode.lower().strip()
    if mode not in {"genes", "average"}:
        raise ValueError("mode must be either 'genes' or 'average'")

    if mode == "genes":
        if add_model3_legend and has_model3_variants:
            fig.subplots_adjust(right=0.58)
        else:
            fig.subplots_adjust(right=0.78)
    else:
        if add_model3_legend and has_model3_variants:
            fig.subplots_adjust(right=0.72)

    if mode == "average":
        d2 = df.copy()
        d2[model_col] = d2[model_col].astype(str)
        d2[gene_col] = d2[gene_col].astype(str)
        d2[y_col] = pd.to_numeric(d2[y_col], errors="coerce")

        rows = []
        for m in model_order:
            sub = d2[d2[model_col] == str(m)][[gene_col, y_col]].copy()
            sub = sub[np.isfinite(sub[y_col])]

            if sub.empty:
                rows.append((str(m), np.nan, np.nan, np.nan, 0, None, None))
                continue

            sub_sorted = sub.sort_values([y_col, gene_col], ascending=[True, True])
            min_gene = sub_sorted.iloc[0][gene_col]

            sub_sorted_desc = sub.sort_values(
                [y_col, gene_col], ascending=[False, True]
            )
            max_gene = sub_sorted_desc.iloc[0][gene_col]

            rows.append(
                (
                    str(m),
                    float(sub[y_col].mean()),
                    float(sub[y_col].min()),
                    float(sub[y_col].max()),
                    int(len(sub)),
                    min_gene,
                    max_gene,
                )
            )

        summary_df = pd.DataFrame(
            rows,
            columns=[
                model_col,
                "mean_adj_r2",
                "min_adj_r2",
                "max_adj_r2",
                "n_genes",
                "min_gene",
                "max_gene",
            ],
        )

        means = summary_df["mean_adj_r2"].to_numpy(dtype=float)
        mins = summary_df["min_adj_r2"].to_numpy(dtype=float)
        maxs = summary_df["max_adj_r2"].to_numpy(dtype=float)

        yerr = np.vstack([means - mins, maxs - means])
        ax.errorbar(x, means, yerr=yerr, fmt="-o", linewidth=2, capsize=4)

        if annotate_minmax_genes:
            for i in range(len(summary_df)):
                if not (np.isfinite(mins[i]) and np.isfinite(maxs[i])):
                    continue
                min_gene = summary_df.loc[i, "min_gene"]
                max_gene = summary_df.loc[i, "max_gene"]

                if isinstance(max_gene, str):
                    ax.annotate(
                        max_gene,
                        xy=(x[i], maxs[i]),
                        xytext=(0, 4),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=annotate_fontsize,
                    )
                if isinstance(min_gene, str):
                    ax.annotate(
                        min_gene,
                        xy=(x[i], mins[i]),
                        xytext=(0, -4),
                        textcoords="offset points",
                        ha="center",
                        va="top",
                        fontsize=annotate_fontsize,
                    )

        ax.set_xticks(x)
        ax.set_xticklabels(x_labels)
        ax.set_ylabel("Adjusted R²")
        ax.set_xlabel("")
        ax.grid(
            True, which="major", axis="both", linestyle="--", linewidth=0.5, alpha=0.6
        )

        if title is None:
            title = "Mean adjusted R² across genes (error bars: min to max)"
        ax.set_title(title)

        if print_summary:
            to_print = summary_df.copy()
            to_print[model_col] = to_print[model_col].map(model_label)
            print(to_print.to_string(index=False))

    else:
        genes = sorted(df[gene_col].astype(str).unique().tolist())

        for gene in genes:
            gdf = df[df[gene_col].astype(str) == gene].copy()
            gdf[model_col] = gdf[model_col].astype(str)
            gdf = gdf.set_index(model_col)

            y = [gdf[y_col].get(str(m), np.nan) for m in model_order]
            ax.plot(x, y, linewidth=1.5, label=gene, zorder=1)

            xs_sig, ys_sig, xs_nsig, ys_nsig = [], [], [], []
            for i, m in enumerate(model_order):
                if str(m) == "0":
                    continue
                yy = gdf[y_col].get(str(m), np.nan)
                pp = gdf[p_used].get(str(m), np.nan)
                if not np.isfinite(yy) or not np.isfinite(pp):
                    continue
                if pp < alpha:
                    xs_sig.append(i)
                    ys_sig.append(yy)
                else:
                    xs_nsig.append(i)
                    ys_nsig.append(yy)

            if xs_nsig:
                ax.scatter(xs_nsig, ys_nsig, marker="s", s=40, c="black", zorder=3)
            if xs_sig:
                ax.scatter(xs_sig, ys_sig, marker="*", s=120, c="red", zorder=4)

        ax.set_xticks(x)
        ax.set_xticklabels(x_labels)
        ax.set_ylabel("Adjusted R²")
        ax.set_xlabel("Nested model")
        ax.grid(
            True, which="major", axis="both", linestyle="--", linewidth=0.5, alpha=0.6
        )

        if title is None:
            title = f"Adjusted R² across nested models (markers by {p_label} < {alpha})"
        ax.set_title(title)

        gene_leg = ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            frameon=False,
            title="Gene",
        )
        ax.add_artist(gene_leg)

    if add_model3_legend and has_model3_variants:
        handles = [Line2D([0], [0], linestyle="none", marker=None) for _ in m3_ks]

        labels = []
        for k in m3_ks:

            if k == 1:
                current_label = "Adds top-1 neighbor PIG mean-expression feature"
            else:
                current_label = f"Adds top-{k} neighbor PIG mean-expression features"
            labels.append(current_label)

        if mode == "genes":
            anchor = (1.02, 0.35)
        else:
            anchor = (1.02, 1.0)

        m3_leg = ax.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=anchor,
            frameon=False,
            title=model3_legend_title,
            handlelength=0,
            handletextpad=0.0,
            borderaxespad=0.0,
            fontsize=9,
        )
        ax.add_artist(m3_leg)

    if savepath is not None:
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
    if show:
        plt.show()

    return fig, ax, summary_df


def compute_gene_spatial_stats(
    df: pd.DataFrame,
    gene_cols: list[str],
    distance_col: str = "distance_to_plaque",
) -> pd.DataFrame:
    """
    For each gene, compute correlation and regression slope vs. plaque distance,
    plus FDR-corrected p-values.
    """
    if distance_col not in df.columns:
        raise ValueError(f"Missing distance column: {distance_col}")

    x = df[[distance_col]].to_numpy()
    results = []

    for g in gene_cols:
        y = df[g].to_numpy()
        if np.allclose(y, 0):
            continue

        r, p = spearmanr(x.ravel(), y)
        if np.isnan(r):
            continue

        model = LinearRegression().fit(x, y)
        slope = float(model.coef_[0])

        results.append((g, r, p, slope))

    res_df = pd.DataFrame(
        results, columns=["gene", "spearman_r", "spearman_p", "slope"]
    )

    if not res_df.empty:
        _, fdr, _, _ = multipletests(res_df["spearman_p"], method="fdr_bh")
        res_df["fdr_pval"] = fdr

    res_df.sort_values("spearman_r", ascending=False, inplace=True)
    logger.info(f"Computed spatial trends for {len(res_df)} genes (FDR-corrected).")
    return res_df


def rank_neighbor_pigs_by_correlation(
    cells_df: pd.DataFrame,
    pig_genes: list[str],
    use_spearman: bool = False,
) -> dict[str, list[tuple[str, float]]]:
    """
    For each target PIG g:
      - looks at columns neigh_mean_{h} for all h != g
      - computes correlation between g (cell expression) and each neigh_mean_h
      - ranks neighbor PIGs by absolute correlation |corr|

    Args:
        cells_df: DataFrame with cell data. Must contain:
            - Columns for each gene in pig_genes: gene expression values
            - Columns neigh_mean_{gene} for each gene in pig_genes: neighbor mean expression
        pig_genes: List of PIG gene names (e.g., ['Gfap', 'Apoe', 'Cst3', ...])
        use_spearman: If True, use Spearman correlation; otherwise use Pearson (default: False)

    Returns:
        rankings: dict mapping gene g -> list of (neighbor_col_name, |corr|)
                  sorted descending by absolute correlation.
    """
    rankings = {}

    present_genes = [g for g in pig_genes if g in cells_df.columns]
    if not present_genes:
        logger.warning("No PIG genes found in cells_df columns.")
        return rankings

    for target_gene in present_genes:

        y = cells_df[target_gene].to_numpy()

        if np.allclose(y, 0) or np.all(np.isnan(y)):
            continue

        neighbor_correlations = []
        for other_gene in present_genes:
            if other_gene == target_gene:
                continue

            neighbor_col = f"neigh_mean_{other_gene}"
            if neighbor_col not in cells_df.columns:
                logger.debug(f"Missing neighbor column: {neighbor_col}")
                continue

            x = cells_df[neighbor_col].to_numpy()

            if np.allclose(x, 0) or np.all(np.isnan(x)):
                continue

            if use_spearman:
                corr, _ = spearmanr(x, y, nan_policy="omit")
            else:

                mask = ~(np.isnan(x) | np.isnan(y))
                if np.sum(mask) < 2:
                    continue
                x_clean = x[mask]
                y_clean = y[mask]
                corr = np.corrcoef(x_clean, y_clean)[0, 1]

            if np.isnan(corr):
                continue

            neighbor_correlations.append((neighbor_col, abs(corr)))

        neighbor_correlations.sort(key=lambda x: x[1], reverse=True)
        rankings[target_gene] = neighbor_correlations

    logger.info(
        f"Ranked neighbor PIGs by correlation for {len(rankings)} target genes "
        f"(using {'Spearman' if use_spearman else 'Pearson'} correlation)."
    )
    return rankings


def regress_expression_with_spatial_features(
    cells_df: pd.DataFrame,
    pig_genes: list[str],
    neighbor_rankings: dict[str, list[tuple[str, float]]] | None = None,
    distance_col: str = "distance_to_plaque",
    k_values: list[int] | None = None,
) -> pd.DataFrame:
    """
    Fit nested regression models for each PIG gene to assess the contribution
    of different spatial feature groups.

    For each PIG gene g, fits a series of nested linear models:
    - Model 0 (baseline): expr_g ~ distance_to_plaque
    - Model 1 (+ plaque geometry): expr_g ~ distance + plaque_geometry_features
    - Model 2 (+ multi-plaque proximity): expr_g ~ distance + plaque_geom + multi_plaque_features
    - Model 3(k) (+ neighborhood PIGs): expr_g ~ distance + plaque_geom +
      multi_plaque + top_k_neighbor_PIGs

    Args:
        cells_df: DataFrame with cell data. Must contain:
            - Columns for each gene in pig_genes: gene expression values
            - distance_col: distance to nearest plaque
            - Plaque geometry features (optional): nearest_plaque_area, nearest_plaque_perimeter,
              nearest_plaque_major_axis, nearest_plaque_minor_axis, nearest_plaque_orientation
            - Multi-plaque proximity features (optional): n_plaques_within_*um,
              mean_dist_to_plaques_within_*um
            - Neighborhood PIG features (optional): neigh_mean_{gene} for each gene in pig_genes
        pig_genes: List of PIG gene names (e.g., ['Gfap', 'Apoe', 'Cst3', ...])
        neighbor_rankings: Dictionary mapping gene g -> list of (neighbor_col_name, |corr|)
            sorted descending by absolute correlation. If None, will be computed automatically.
        distance_col: Name of the column containing distance to nearest plaque
            (default: 'distance_to_plaque')
        k_values: List of k values for Model 3(k) (default: [1, 2, 4, 8, 15])

    Returns:
        DataFrame with columns:
            - gene: PIG gene name
            - model: Model identifier (0, 1, 2, or 3_k)
            - r_squared: R² value
            - adj_r_squared: Adjusted R² value
            - n_obs: Number of observations used
            - n_params: Number of parameters in the model
            - f_test_pval: P-value from nested F-test (comparing to previous model, NaN for Model 0)
            - f_test_fstat: F-statistic from nested F-test (NaN for Model 0)
    """
    if k_values is None:
        k_values = [1, 2, 4, 8, 15]

    if distance_col not in cells_df.columns:
        raise ValueError(
            f"Missing distance column: '{distance_col}'. "
            f"Available columns: {cells_df.columns.tolist()[:20]}..."
        )

    present_genes = [g for g in pig_genes if g in cells_df.columns]
    if not present_genes:
        logger.warning("No PIG genes found in cells_df columns.")
        return pd.DataFrame()

    if neighbor_rankings is None:
        logger.info("Computing neighbor PIG rankings...")
        neighbor_rankings = rank_neighbor_pigs_by_correlation(
            cells_df=cells_df,
            pig_genes=present_genes,
            use_spearman=False,
        )

    geometry_features = [
        "nearest_plaque_area",
        "nearest_plaque_perimeter",
        "nearest_plaque_major_axis",
        "nearest_plaque_minor_axis",
        "nearest_plaque_orientation",
    ]
    available_geometry = [f for f in geometry_features if f in cells_df.columns]

    multi_plaque_features = [
        col
        for col in cells_df.columns
        if col.startswith("n_plaques_within_")
        or col.startswith("mean_dist_to_plaques_within_")
    ]

    logger.info(
        f"Fitting nested regression models for {len(present_genes)} PIG genes. "
        f"Available features: {len(available_geometry)} geometry, "
        f"{len(multi_plaque_features)} multi-plaque proximity"
    )

    results = []

    for gene in present_genes:

        y = cells_df[gene].to_numpy()

        if np.allclose(y, 0) or np.all(np.isnan(y)):
            logger.debug(f"Skipping {gene}: constant or all NaN")
            continue

        valid_mask = ~(np.isnan(y) | np.isnan(cells_df[distance_col].values))
        if np.sum(valid_mask) < 10:
            logger.debug(
                f"Skipping {gene}: insufficient valid observations ({np.sum(valid_mask)})"
            )
            continue

        y_valid = y[valid_mask]
        n_obs = len(y_valid)

        fitted_models = {}

        try:
            x0 = cells_df[[distance_col]].values[valid_mask]
            x0 = sm.add_constant(x0)
            model0 = sm.OLS(y_valid, x0).fit()
            fitted_models[0] = model0

            results.append(
                {
                    "gene": gene,
                    "model": "0",
                    "r_squared": model0.rsquared,
                    "adj_r_squared": model0.rsquared_adj,
                    "n_obs": n_obs,
                    "n_params": model0.df_model + 1,
                    "f_test_pval": np.nan,
                    "f_test_fstat": np.nan,
                }
            )
        except Exception as e:
            logger.warning(f"Failed to fit Model 0 for {gene}: {e}")
            continue

        if available_geometry:
            try:
                x1_cols = [distance_col, *available_geometry]
                x1 = cells_df[x1_cols].values[valid_mask]

                geom_mask = ~np.isnan(x1).any(axis=1)
                if np.sum(geom_mask) < 10:
                    logger.debug(
                        f"Skipping Model 1 for {gene}: insufficient valid "
                        f"observations after geometry filtering"
                    )
                else:

                    x0_common = cells_df[[distance_col]].values[valid_mask][geom_mask]
                    y_common = y_valid[geom_mask]
                    x0_common = sm.add_constant(x0_common)
                    model0_common = sm.OLS(y_common, x0_common).fit()

                    x1_clean = x1[geom_mask]
                    y1_clean = y_valid[geom_mask]
                    x1_clean = sm.add_constant(x1_clean)
                    model1 = sm.OLS(y1_clean, x1_clean).fit()
                    fitted_models[1] = model1

                    f_test = model1.compare_f_test(model0_common)
                    f_stat = f_test[0]
                    f_pval = f_test[1]

                    results.append(
                        {
                            "gene": gene,
                            "model": "1",
                            "r_squared": model1.rsquared,
                            "adj_r_squared": model1.rsquared_adj,
                            "n_obs": np.sum(geom_mask),
                            "n_params": model1.df_model + 1,
                            "f_test_pval": f_pval,
                            "f_test_fstat": f_stat,
                        }
                    )
            except Exception as e:
                logger.debug(f"Failed to fit Model 1 for {gene}: {e}")

        if multi_plaque_features:
            try:

                if 1 in fitted_models and available_geometry:
                    x2_cols = [
                        distance_col,
                        *available_geometry,
                        *multi_plaque_features,
                    ]
                    base_cols = [distance_col, *available_geometry]
                else:
                    x2_cols = [distance_col, *multi_plaque_features]
                    base_cols = [distance_col]

                x2 = cells_df[x2_cols].values[valid_mask]

                multi_mask = ~np.isnan(x2).any(axis=1)
                if np.sum(multi_mask) < 10:
                    logger.debug(
                        f"Skipping Model 2 for {gene}: insufficient valid "
                        f"observations after multi-plaque filtering"
                    )
                else:

                    x_base_common = cells_df[base_cols].values[valid_mask][multi_mask]
                    y_common = y_valid[multi_mask]
                    x_base_common = sm.add_constant(x_base_common)
                    base_model_common = sm.OLS(y_common, x_base_common).fit()

                    x2_clean = x2[multi_mask]
                    y2_clean = y_valid[multi_mask]
                    x2_clean = sm.add_constant(x2_clean)
                    model2 = sm.OLS(y2_clean, x2_clean).fit()
                    fitted_models[2] = model2

                    f_test = model2.compare_f_test(base_model_common)
                    f_stat = f_test[0]
                    f_pval = f_test[1]

                    results.append(
                        {
                            "gene": gene,
                            "model": "2",
                            "r_squared": model2.rsquared,
                            "adj_r_squared": model2.rsquared_adj,
                            "n_obs": np.sum(multi_mask),
                            "n_params": model2.df_model + 1,
                            "f_test_pval": f_pval,
                            "f_test_fstat": f_stat,
                        }
                    )
            except Exception as e:
                logger.debug(f"Failed to fit Model 2 for {gene}: {e}")

        if neighbor_rankings.get(gene):

            if 2 in fitted_models:
                if available_geometry:
                    base_cols = [
                        distance_col,
                        *available_geometry,
                        *multi_plaque_features,
                    ]
                else:
                    base_cols = [distance_col, *multi_plaque_features]
            elif 1 in fitted_models and available_geometry:
                base_cols = [distance_col, *available_geometry]
            else:
                base_cols = [distance_col]

            for k in k_values:

                top_k_neighbors = [col for col, _ in neighbor_rankings[gene][:k]]

                top_k_neighbors = [
                    col for col in top_k_neighbors if col in cells_df.columns
                ]

                if not top_k_neighbors:
                    continue

                try:
                    x3_cols = [*base_cols, *top_k_neighbors]
                    x3 = cells_df[x3_cols].values[valid_mask]

                    neigh_mask = ~np.isnan(x3).any(axis=1)
                    if np.sum(neigh_mask) < 10:
                        logger.debug(
                            f"Skipping Model 3(k={k}) for {gene}: insufficient "
                            f"valid observations after neighbor filtering"
                        )
                        continue

                    x_base_common = cells_df[base_cols].values[valid_mask][neigh_mask]
                    y_common = y_valid[neigh_mask]
                    x_base_common = sm.add_constant(x_base_common)
                    base_model_common = sm.OLS(y_common, x_base_common).fit()

                    x3_clean = x3[neigh_mask]
                    y3_clean = y_valid[neigh_mask]
                    x3_clean = sm.add_constant(x3_clean)
                    model3k = sm.OLS(y3_clean, x3_clean).fit()
                    fitted_models[f"3_{k}"] = model3k

                    f_test = model3k.compare_f_test(base_model_common)
                    f_stat = f_test[0]
                    f_pval = f_test[1]

                    results.append(
                        {
                            "gene": gene,
                            "model": f"3_{k}",
                            "r_squared": model3k.rsquared,
                            "adj_r_squared": model3k.rsquared_adj,
                            "n_obs": np.sum(neigh_mask),
                            "n_params": model3k.df_model + 1,
                            "f_test_pval": f_pval,
                            "f_test_fstat": f_stat,
                        }
                    )
                except Exception as e:
                    logger.debug(f"Failed to fit Model 3(k={k}) for {gene}: {e}")

    results_df = pd.DataFrame(results)

    if not results_df.empty:
        logger.info(
            f" Completed nested regression analysis. "
            f"Fitted models for {results_df['gene'].nunique()} genes, "
            f"total {len(results_df)} model fits."
        )
    else:
        logger.warning("No regression models were successfully fitted.")

    return results_df


def separate_nested_regression_results(
    results_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Separate nested regression results into two groups:
    (A) Model summaries: Each model with its features, R², adj R², n_obs
    (B) Nested comparisons: Comparisons between models with ΔR² and p-values

    Args:
        results_df: DataFrame from regress_expression_with_spatial_features()
            Must contain columns: gene, model, r_squared, adj_r_squared, n_obs, n_params,
            f_test_pval, f_test_fstat

    Returns:
        Tuple of (model_summaries_df, nested_comparisons_df):
        - model_summaries_df: Columns: gene, model, features_used, r_squared, adj_r_squared,
          n_obs, n_params
        - nested_comparisons_df: Columns: gene, comparison, model_from, model_to, delta_r_squared,
          delta_adj_r_squared, f_test_pval, f_test_fstat
    """
    if results_df.empty:
        logger.warning("Empty results DataFrame provided.")
        return pd.DataFrame(), pd.DataFrame()

    feature_descriptions = {
        "0": "distance_to_plaque",
        "1": (
            "distance_to_plaque + plaque_geometry "
            "(area, perimeter, major_axis, minor_axis, orientation)"
        ),
        "2": (
            "distance_to_plaque + plaque_geometry + multi_plaque_proximity "
            "(n_plaques_within_radius, mean_dist_to_plaques)"
        ),
    }

    model_summaries = results_df[
        ["gene", "model", "r_squared", "adj_r_squared", "n_obs", "n_params"]
    ].copy()

    def get_features(model_name: str) -> str:
        """Get feature description for a model."""
        if model_name == "0":
            return feature_descriptions["0"]
        elif model_name == "1":
            return feature_descriptions["1"]
        elif model_name == "2":
            return feature_descriptions["2"]
        elif model_name.startswith("3_"):
            k = model_name.replace("3_", "")
            return (
                feature_descriptions["2"]
                + f" + top_{k}_neighbor_PIG_features (ranked by correlation)"
            )
        else:
            return "unknown"

    model_summaries["features_used"] = model_summaries["model"].apply(get_features)

    model_summaries = model_summaries[
        [
            "gene",
            "model",
            "features_used",
            "r_squared",
            "adj_r_squared",
            "n_obs",
            "n_params",
        ]
    ]

    comparisons = []

    for gene in results_df["gene"].unique():
        gene_results = results_df[results_df["gene"] == gene].copy()

        def model_sort_key(model_name: str) -> tuple[int, int]:
            """Sort key: (model_level, k_value)."""
            if model_name == "0":
                return (0, 0)
            elif model_name == "1":
                return (1, 0)
            elif model_name == "2":
                return (2, 0)
            elif model_name.startswith("3_"):
                k = int(model_name.replace("3_", ""))
                return (3, k)
            else:
                return (999, 0)

        gene_results["_sort_key"] = gene_results["model"].apply(model_sort_key)
        gene_results = gene_results.sort_values("_sort_key").drop(columns=["_sort_key"])

        for i in range(len(gene_results) - 1):
            model_from_row = gene_results.iloc[i]
            model_to_row = gene_results.iloc[i + 1]

            model_from = model_from_row["model"]
            model_to = model_to_row["model"]

            if pd.isna(model_to_row["f_test_pval"]):
                logger.debug(
                    f"Skipping comparison {model_from}→{model_to} for {gene}: "
                    "missing F-test p-value"
                )
                continue

            delta_r_squared = model_to_row["r_squared"] - model_from_row["r_squared"]
            delta_adj_r_squared = (
                model_to_row["adj_r_squared"] - model_from_row["adj_r_squared"]
            )

            comparison_label = f"M{model_from}→M{model_to}"

            comparisons.append(
                {
                    "gene": gene,
                    "comparison": comparison_label,
                    "model_from": model_from,
                    "model_to": model_to,
                    "delta_r_squared": delta_r_squared,
                    "delta_adj_r_squared": delta_adj_r_squared,
                    "f_test_pval": model_to_row["f_test_pval"],
                    "f_test_fstat": model_to_row["f_test_fstat"],
                }
            )

    nested_comparisons = pd.DataFrame(comparisons)

    if not nested_comparisons.empty:

        def comparison_sort_key(row: pd.Series) -> tuple[str, int, int]:
            """Sort key: (gene, model_level_from, model_level_to)."""
            model_from = row["model_from"]
            model_to = row["model_to"]

            def get_level(model_name: str) -> int:
                if model_name == "0":
                    return 0
                elif model_name == "1":
                    return 1
                elif model_name == "2":
                    return 2
                elif model_name.startswith("3_"):
                    return 3
                else:
                    return 999

            return (row["gene"], get_level(model_from), get_level(model_to))

        nested_comparisons["_sort_key"] = nested_comparisons.apply(
            comparison_sort_key, axis=1
        )
        nested_comparisons = nested_comparisons.sort_values("_sort_key").drop(
            columns=["_sort_key"]
        )

    logger.info(
        f" Separated results: {len(model_summaries)} model summaries, "
        f"{len(nested_comparisons)} nested comparisons"
    )

    return model_summaries, nested_comparisons


def apply_fdr_correction_to_comparisons(
    nested_comparisons: pd.DataFrame,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Apply Benjamini-Hochberg FDR correction to nested comparisons, grouped by comparison type.

    For each comparison type (e.g., M0→M1, M1→M2, M2→M3_1, etc.), collects all p-values
    and applies FDR correction separately. This accounts for multiple testing across genes.

    Args:
        nested_comparisons: DataFrame from separate_nested_regression_results()
            Must contain columns: comparison, f_test_pval
        alpha: FDR threshold (default: 0.05)

    Returns:
        DataFrame with added column: f_test_pval_adj (FDR-adjusted p-values)
    """
    if nested_comparisons.empty:
        logger.warning("Empty comparisons DataFrame provided.")
        return nested_comparisons.copy()

    comparisons_adj = nested_comparisons.copy()

    for comparison_type in comparisons_adj["comparison"].unique():
        mask = comparisons_adj["comparison"] == comparison_type
        pvals = comparisons_adj.loc[mask, "f_test_pval"].values

        valid_mask = ~pd.isna(pvals)
        if np.sum(valid_mask) == 0:
            continue

        pvals_valid = pvals[valid_mask]

        _, pvals_adj, _, _ = multipletests(pvals_valid, method="fdr_bh", alpha=alpha)

        pvals_adj_full = np.full(len(pvals), np.nan)
        pvals_adj_full[valid_mask] = pvals_adj

        comparisons_adj.loc[mask, "f_test_pval_adj"] = pvals_adj_full

    n_sig_before = (comparisons_adj["f_test_pval"] < alpha).sum()
    n_sig_after = (comparisons_adj["f_test_pval_adj"] < alpha).sum()

    logger.info(
        f" Applied FDR correction. "
        f"Significant before: {n_sig_before}, after: {n_sig_after} "
        f"(alpha = {alpha})"
    )

    return comparisons_adj


def create_feature_block_importance_table(
    nested_comparisons: pd.DataFrame,
    alpha: float = 0.05,
    use_fdr_corrected: bool = True,
) -> pd.DataFrame:
    """
    Create summary table showing feature block importance across all PIG genes.

    For each comparison type (M0→M1, M1→M2, M2→M3_k), computes:
    - Mean ΔR² across all genes
    - Median ΔR²
    - Number of significant genes (p < alpha)
    - Total genes tested

    Args:
        nested_comparisons: DataFrame from separate_nested_regression_results()
            or apply_fdr_correction_to_comparisons()
        alpha: Significance threshold (default: 0.05)
        use_fdr_corrected: If True, use f_test_pval_adj; otherwise use f_test_pval

    Returns:
        DataFrame with columns:
            - comparison: Comparison type (e.g., "M0→M1")
            - mean_delta_r_squared: Mean ΔR² across genes
            - median_delta_r_squared: Median ΔR² across genes
            - n_significant: Number of genes with significant improvement
            - n_tested: Total number of genes tested
            - pct_significant: Percentage of genes with significant improvement
    """
    if nested_comparisons.empty:
        logger.warning("Empty comparisons DataFrame provided.")
        return pd.DataFrame()

    pval_col = (
        "f_test_pval_adj"
        if (use_fdr_corrected and "f_test_pval_adj" in nested_comparisons.columns)
        else "f_test_pval"
    )

    if use_fdr_corrected and pval_col not in nested_comparisons.columns:
        logger.warning(
            "f_test_pval_adj not found. Using uncorrected p-values. "
            "Run apply_fdr_correction_to_comparisons() first."
        )
        pval_col = "f_test_pval"

    summary_rows = []

    for comparison_type in sorted(nested_comparisons["comparison"].unique()):
        comp_data = nested_comparisons[
            nested_comparisons["comparison"] == comparison_type
        ]

        comp_data_valid = comp_data[comp_data[pval_col].notna()]

        if len(comp_data_valid) == 0:
            continue

        mean_delta_r2 = comp_data_valid["delta_r_squared"].mean()
        median_delta_r2 = comp_data_valid["delta_r_squared"].median()
        mean_delta_adj_r2 = comp_data_valid["delta_adj_r_squared"].mean()
        median_delta_adj_r2 = comp_data_valid["delta_adj_r_squared"].median()
        n_tested = len(comp_data_valid)
        n_significant = (comp_data_valid[pval_col] < alpha).sum()
        pct_significant = 100 * n_significant / n_tested if n_tested > 0 else 0.0

        summary_rows.append(
            {
                "comparison": comparison_type,
                "mean_delta_r_squared": mean_delta_r2,
                "median_delta_r_squared": median_delta_r2,
                "n_significant": n_significant,
                "n_tested": n_tested,
                "pct_significant": pct_significant,
                "mean_delta_adj_r_squared": mean_delta_adj_r2,
                "median_delta_adj_r_squared": median_delta_adj_r2,
            }
        )

    summary_df = pd.DataFrame(summary_rows)

    logger.info(
        f" Created feature block importance table with {len(summary_df)} comparison types"
    )

    return summary_df


def create_best_neighbor_pigs_table(
    nested_comparisons: pd.DataFrame,
    neighbor_rankings: dict[str, list[tuple[str, float]]],
    alpha: float = 0.05,
    use_fdr_corrected: bool = True,
    score_col: str = "delta_adj_r_squared",
) -> pd.DataFrame:
    if nested_comparisons.empty:
        logger.warning("Empty comparisons DataFrame provided.")
        return pd.DataFrame()

    df = nested_comparisons.copy()

    pval_col = (
        "f_test_pval_adj"
        if (use_fdr_corrected and "f_test_pval_adj" in df.columns)
        else "f_test_pval"
    )
    if pval_col not in df.columns:
        raise ValueError(f"Missing p-value column '{pval_col}' in nested_comparisons.")
    if score_col not in df.columns:
        raise ValueError(f"Missing score column '{score_col}' in nested_comparisons.")

    df["model_from"] = df["model_from"].astype(str)
    df["model_to"] = df["model_to"].astype(str)

    def parse_k(m: str) -> int:
        m = str(m)
        if m.startswith("3_"):
            try:
                return int(m.split("_", 1)[1])
            except Exception:
                return 0
        return 0

    df = df[
        df["model_to"].str.startswith("3_")
        & (df["model_from"].eq("2") | df["model_from"].str.startswith("3_"))
    ].copy()
    if df.empty:
        logger.warning("No (2 or 3_k) → 3_k comparisons found.")
        return pd.DataFrame()

    df["k_from"] = df["model_from"].apply(lambda m: 0 if m == "2" else parse_k(m))
    df["k_to"] = df["model_to"].apply(parse_k)
    df = df[df["k_to"] > df["k_from"]].copy()
    if df.empty:
        logger.warning("No increasing-k transitions found.")
        return pd.DataFrame()

    df["delta_r_squared"] = pd.to_numeric(df.get("delta_r_squared"), errors="coerce")
    df["delta_adj_r_squared"] = pd.to_numeric(
        df.get("delta_adj_r_squared"), errors="coerce"
    )

    best_rows = []

    for gene, g in df.groupby("gene", sort=False):
        g = g.copy()

        g_score = pd.to_numeric(g[score_col], errors="coerce")
        g = g.assign(_score=g_score).dropna(subset=["_score"])
        if g.empty:
            continue

        sig = g[g[pval_col] < alpha]
        if not sig.empty:
            best = sig.loc[sig["_score"].idxmax()]
            is_sig = True
        else:
            best = g.loc[g["_score"].idxmax()]
            is_sig = False

        best_k = int(best["k_to"])
        best_model = str(best["model_to"])
        best_step = f"M{best['model_from']}→M{best['model_to']}"

        cum_r2 = {0: 0.0}
        cum_adj = {0: 0.0}

        g_sorted = g.sort_values(["k_to", "k_from"]).to_dict("records")
        for row in g_sorted:
            kf, kt = int(row["k_from"]), int(row["k_to"])
            if kf not in cum_r2:
                continue
            dr2 = row.get("delta_r_squared")
            da2 = row.get("delta_adj_r_squared")
            if pd.notna(dr2):
                cum_r2[kt] = float(cum_r2[kf] + float(dr2))
            if pd.notna(da2):
                cum_adj[kt] = float(cum_adj[kf] + float(da2))

        total_delta_r2 = cum_r2.get(best_k, np.nan)
        total_delta_adj = cum_adj.get(best_k, np.nan)

        added_cols, added_genes = [], []
        topk_cols, topk_genes = [], []
        if gene in neighbor_rankings and best_k > 0:
            ranked_cols = [c for c, _ in neighbor_rankings[gene]]
            k_from = int(best["k_from"])
            added_cols = ranked_cols[k_from:best_k]
            added_genes = [c.replace("neigh_mean_", "") for c in added_cols]
            topk_cols = ranked_cols[:best_k]
            topk_genes = [c.replace("neigh_mean_", "") for c in topk_cols]

        best_rows.append(
            {
                "gene": gene,
                "best_k": best_k,
                "best_model": best_model,
                "best_step": best_step,
                "score": float(best["_score"]),
                "score_metric": score_col,
                "delta_r_squared": total_delta_r2,
                "delta_adj_r_squared": total_delta_adj,
                "step_delta_r_squared": float(best.get("delta_r_squared", np.nan)),
                "step_delta_adj_r_squared": float(
                    best.get("delta_adj_r_squared", np.nan)
                ),
                "f_test_pval": float(best[pval_col]),
                "is_significant": bool(is_sig),
                "added_neighbor_genes": ", ".join(added_genes),
                "added_neighbor_cols": ", ".join(added_cols),
                "top_k_neighbor_genes": ", ".join(topk_genes),
                "top_k_neighbor_cols": ", ".join(topk_cols),
            }
        )

    out = pd.DataFrame(best_rows)
    if out.empty:
        return out

    out = out.sort_values("score", ascending=False, kind="stable").reset_index(
        drop=True
    )
    return out
