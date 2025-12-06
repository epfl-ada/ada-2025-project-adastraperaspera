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
            continue  # skip constant genes

        # Spearman correlation
        r, p = spearmanr(x.ravel(), y)
        if np.isnan(r):
            continue

        # Linear regression slope
        model = LinearRegression().fit(x, y)
        slope = float(model.coef_[0])

        results.append((g, r, p, slope))

    res_df = pd.DataFrame(results, columns=["gene", "spearman_r", "spearman_p", "slope"])

    # --- Apply FDR correction ---
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

    # Filter to genes that are actually present in the dataframe
    present_genes = [g for g in pig_genes if g in cells_df.columns]
    if not present_genes:
        logger.warning("No PIG genes found in cells_df columns.")
        return rankings

    for target_gene in present_genes:
        # Get target gene expression
        y = cells_df[target_gene].to_numpy()

        # Skip if constant or all NaN
        if np.allclose(y, 0) or np.all(np.isnan(y)):
            continue

        # Find all neighbor mean columns for other PIG genes
        neighbor_correlations = []
        for other_gene in present_genes:
            if other_gene == target_gene:
                continue  # Skip self

            neighbor_col = f"neigh_mean_{other_gene}"
            if neighbor_col not in cells_df.columns:
                logger.debug(f"Missing neighbor column: {neighbor_col}")
                continue

            # Get neighbor mean expression
            x = cells_df[neighbor_col].to_numpy()

            # Skip if constant or all NaN
            if np.allclose(x, 0) or np.all(np.isnan(x)):
                continue

            # Compute correlation
            if use_spearman:
                corr, _ = spearmanr(x, y, nan_policy="omit")
            else:
                # Pearson correlation
                # Remove NaN pairs
                mask = ~(np.isnan(x) | np.isnan(y))
                if np.sum(mask) < 2:
                    continue  # Need at least 2 valid pairs
                x_clean = x[mask]
                y_clean = y[mask]
                corr = np.corrcoef(x_clean, y_clean)[0, 1]

            if np.isnan(corr):
                continue

            # Store absolute correlation
            neighbor_correlations.append((neighbor_col, abs(corr)))

        # Sort by absolute correlation (descending)
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

    # Filter to genes that are actually present in the dataframe
    present_genes = [g for g in pig_genes if g in cells_df.columns]
    if not present_genes:
        logger.warning("No PIG genes found in cells_df columns.")
        return pd.DataFrame()

    # Compute neighbor rankings if not provided
    if neighbor_rankings is None:
        logger.info("Computing neighbor PIG rankings...")
        neighbor_rankings = rank_neighbor_pigs_by_correlation(
            cells_df=cells_df,
            pig_genes=present_genes,
            use_spearman=False,
        )

    # Identify available feature groups
    geometry_features = [
        "nearest_plaque_area",
        "nearest_plaque_perimeter",
        "nearest_plaque_major_axis",
        "nearest_plaque_minor_axis",
        "nearest_plaque_orientation",
    ]
    available_geometry = [f for f in geometry_features if f in cells_df.columns]

    # Find multi-plaque proximity features
    # (columns starting with n_plaques_within_ or mean_dist_to_plaques_within_)
    multi_plaque_features = [
        col
        for col in cells_df.columns
        if col.startswith("n_plaques_within_") or col.startswith("mean_dist_to_plaques_within_")
    ]

    logger.info(
        f"Fitting nested regression models for {len(present_genes)} PIG genes. "
        f"Available features: {len(available_geometry)} geometry, "
        f"{len(multi_plaque_features)} multi-plaque proximity"
    )

    results = []

    for gene in present_genes:
        # Get target gene expression
        y = cells_df[gene].to_numpy()

        # Skip if constant or all NaN
        if np.allclose(y, 0) or np.all(np.isnan(y)):
            logger.debug(f"Skipping {gene}: constant or all NaN")
            continue

        # Create mask for valid observations (non-NaN in y and distance)
        valid_mask = ~(np.isnan(y) | np.isnan(cells_df[distance_col].values))
        if np.sum(valid_mask) < 10:  # Need at least 10 observations
            logger.debug(f"Skipping {gene}: insufficient valid observations ({np.sum(valid_mask)})")
            continue

        y_valid = y[valid_mask]
        n_obs = len(y_valid)

        # Store fitted models for nested F-tests
        fitted_models = {}

        # Model 0: Baseline (distance only)
        try:
            x0 = cells_df[[distance_col]].values[valid_mask]
            x0 = sm.add_constant(x0)  # Add intercept
            model0 = sm.OLS(y_valid, x0).fit()
            fitted_models[0] = model0

            results.append(
                {
                    "gene": gene,
                    "model": "0",
                    "r_squared": model0.rsquared,
                    "adj_r_squared": model0.rsquared_adj,
                    "n_obs": n_obs,
                    "n_params": model0.df_model + 1,  # +1 for intercept
                    "f_test_pval": np.nan,
                    "f_test_fstat": np.nan,
                }
            )
        except Exception as e:
            logger.warning(f"Failed to fit Model 0 for {gene}: {e}")
            continue

        # Model 1: Add plaque geometry features
        if available_geometry:
            try:
                x1_cols = [distance_col, *available_geometry]
                x1 = cells_df[x1_cols].values[valid_mask]
                # Check for NaN in geometry features
                geom_mask = ~np.isnan(x1).any(axis=1)
                if np.sum(geom_mask) < 10:
                    logger.debug(
                        f"Skipping Model 1 for {gene}: insufficient valid "
                        f"observations after geometry filtering"
                    )
                else:
                    # Use the same observations for Model 0 (for nested F-test)
                    x0_common = cells_df[[distance_col]].values[valid_mask][geom_mask]
                    y_common = y_valid[geom_mask]
                    x0_common = sm.add_constant(x0_common)
                    model0_common = sm.OLS(y_common, x0_common).fit()

                    x1_clean = x1[geom_mask]
                    y1_clean = y_valid[geom_mask]
                    x1_clean = sm.add_constant(x1_clean)
                    model1 = sm.OLS(y1_clean, x1_clean).fit()
                    fitted_models[1] = model1

                    # Nested F-test: Model 1 vs Model 0 (using same observations)
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

        # Model 2: Add multi-plaque proximity features
        if multi_plaque_features:
            try:
                # Start with Model 1 features if available, otherwise Model 0
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
                # Check for NaN in multi-plaque features
                multi_mask = ~np.isnan(x2).any(axis=1)
                if np.sum(multi_mask) < 10:
                    logger.debug(
                        f"Skipping Model 2 for {gene}: insufficient valid "
                        f"observations after multi-plaque filtering"
                    )
                else:
                    # Use the same observations for base model (for nested F-test)
                    x_base_common = cells_df[base_cols].values[valid_mask][multi_mask]
                    y_common = y_valid[multi_mask]
                    x_base_common = sm.add_constant(x_base_common)
                    base_model_common = sm.OLS(y_common, x_base_common).fit()

                    x2_clean = x2[multi_mask]
                    y2_clean = y_valid[multi_mask]
                    x2_clean = sm.add_constant(x2_clean)
                    model2 = sm.OLS(y2_clean, x2_clean).fit()
                    fitted_models[2] = model2

                    # Nested F-test: Model 2 vs base model (using same observations)
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

        # Model 3(k): Add top-k neighborhood PIG features
        if neighbor_rankings.get(gene):
            # Determine base model columns for Model 3
            # (prefer Model 2, then Model 1, then Model 0)
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
                # Get top-k neighbor features
                top_k_neighbors = [col for col, _ in neighbor_rankings[gene][:k]]
                # Filter to features that actually exist in the dataframe
                top_k_neighbors = [col for col in top_k_neighbors if col in cells_df.columns]

                if not top_k_neighbors:
                    continue

                try:
                    x3_cols = [*base_cols, *top_k_neighbors]
                    x3 = cells_df[x3_cols].values[valid_mask]
                    # Check for NaN in neighbor features
                    neigh_mask = ~np.isnan(x3).any(axis=1)
                    if np.sum(neigh_mask) < 10:
                        logger.debug(
                            f"Skipping Model 3(k={k}) for {gene}: insufficient "
                            f"valid observations after neighbor filtering"
                        )
                        continue

                    # Use the same observations for base model (for nested F-test)
                    x_base_common = cells_df[base_cols].values[valid_mask][neigh_mask]
                    y_common = y_valid[neigh_mask]
                    x_base_common = sm.add_constant(x_base_common)
                    base_model_common = sm.OLS(y_common, x_base_common).fit()

                    x3_clean = x3[neigh_mask]
                    y3_clean = y_valid[neigh_mask]
                    x3_clean = sm.add_constant(x3_clean)
                    model3k = sm.OLS(y3_clean, x3_clean).fit()
                    fitted_models[f"3_{k}"] = model3k

                    # Nested F-test: Model 3(k) vs base model (using same observations)
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
            f"✅ Completed nested regression analysis. "
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

    # Define feature descriptions for each model
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

    # (A) Model Summaries: Extract all models, remove F-test columns
    model_summaries = results_df[
        ["gene", "model", "r_squared", "adj_r_squared", "n_obs", "n_params"]
    ].copy()

    # Add feature descriptions
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

    # Reorder columns
    model_summaries = model_summaries[
        ["gene", "model", "features_used", "r_squared", "adj_r_squared", "n_obs", "n_params"]
    ]

    # (B) Nested Comparisons: Create comparison rows
    comparisons = []

    for gene in results_df["gene"].unique():
        gene_results = results_df[results_df["gene"] == gene].copy()

        # Sort by model number for proper ordering
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

        # Create comparisons between consecutive models
        for i in range(len(gene_results) - 1):
            model_from_row = gene_results.iloc[i]
            model_to_row = gene_results.iloc[i + 1]

            model_from = model_from_row["model"]
            model_to = model_to_row["model"]

            # Skip if F-test p-value is NaN (shouldn't happen, but handle gracefully)
            if pd.isna(model_to_row["f_test_pval"]):
                logger.debug(
                    f"Skipping comparison {model_from}→{model_to} for {gene}: "
                    "missing F-test p-value"
                )
                continue

            # Calculate ΔR²
            delta_r_squared = model_to_row["r_squared"] - model_from_row["r_squared"]
            delta_adj_r_squared = model_to_row["adj_r_squared"] - model_from_row["adj_r_squared"]

            # Create comparison label
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

    # Sort comparisons for readability
    if not nested_comparisons.empty:
        # Sort by gene, then by comparison order
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

        nested_comparisons["_sort_key"] = nested_comparisons.apply(comparison_sort_key, axis=1)
        nested_comparisons = nested_comparisons.sort_values("_sort_key").drop(columns=["_sort_key"])

    logger.info(
        f"✅ Separated results: {len(model_summaries)} model summaries, "
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

    # Group by comparison type and apply FDR correction separately
    for comparison_type in comparisons_adj["comparison"].unique():
        mask = comparisons_adj["comparison"] == comparison_type
        pvals = comparisons_adj.loc[mask, "f_test_pval"].values

        # Filter out NaN values for correction
        valid_mask = ~pd.isna(pvals)
        if np.sum(valid_mask) == 0:
            continue

        pvals_valid = pvals[valid_mask]

        # Apply Benjamini-Hochberg FDR correction
        _, pvals_adj, _, _ = multipletests(pvals_valid, method="fdr_bh", alpha=alpha)

        # Create adjusted p-values array (NaN for invalid, adjusted for valid)
        pvals_adj_full = np.full(len(pvals), np.nan)
        pvals_adj_full[valid_mask] = pvals_adj

        # Store adjusted p-values
        comparisons_adj.loc[mask, "f_test_pval_adj"] = pvals_adj_full

    # Count significant after correction
    n_sig_before = (comparisons_adj["f_test_pval"] < alpha).sum()
    n_sig_after = (comparisons_adj["f_test_pval_adj"] < alpha).sum()

    logger.info(
        f"✅ Applied FDR correction. "
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

    # Determine which p-value column to use
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

    # Group by comparison type
    summary_rows = []

    for comparison_type in sorted(nested_comparisons["comparison"].unique()):
        comp_data = nested_comparisons[nested_comparisons["comparison"] == comparison_type]

        # Filter out NaN p-values
        comp_data_valid = comp_data[comp_data[pval_col].notna()]

        if len(comp_data_valid) == 0:
            continue

        # Compute statistics
        mean_delta_r2 = comp_data_valid["delta_r_squared"].mean()
        median_delta_r2 = comp_data_valid["delta_r_squared"].median()
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
            }
        )

    summary_df = pd.DataFrame(summary_rows)

    logger.info(
        f"✅ Created feature block importance table with {len(summary_df)} comparison types"
    )

    return summary_df


def create_best_neighbor_pigs_table(
    nested_comparisons: pd.DataFrame,
    neighbor_rankings: dict[str, list[tuple[str, float]]],
    alpha: float = 0.05,
    use_fdr_corrected: bool = True,
) -> pd.DataFrame:
    """
    Create table showing best neighbor PIG features for each target gene.

    For each target PIG gene, looks at all M2→M3_k comparisons and identifies:
    - The model where ΔR² is largest AND significant
    - The value of k (number of neighbor features)
    - Which neighbor PIG genes were included in the top-k
    - The amount of ΔR² they explained

    Args:
        nested_comparisons: DataFrame from separate_nested_regression_results()
            or apply_fdr_correction_to_comparisons()
        neighbor_rankings: Dictionary from rank_neighbor_pigs_by_correlation()
            Maps gene -> list of (neighbor_col_name, |corr|) sorted descending
        alpha: Significance threshold (default: 0.05)
        use_fdr_corrected: If True, use f_test_pval_adj; otherwise use f_test_pval

    Returns:
        DataFrame with columns:
            - gene: Target PIG gene name
            - best_k: Optimal k value (number of neighbor features)
            - best_model: Model identifier (e.g., "3_4")
            - delta_r_squared: ΔR² explained by neighbor features
            - delta_adj_r_squared: Δadj R² explained
            - f_test_pval: P-value (corrected if available)
            - is_significant: Whether improvement is significant
            - top_k_neighbor_genes: List of neighbor gene names in top-k
            - top_k_neighbor_cols: List of neighbor column names in top-k
    """
    if nested_comparisons.empty:
        logger.warning("Empty comparisons DataFrame provided.")
        return pd.DataFrame()

    # Determine which p-value column to use
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

    # Filter to M2→M3_k comparisons only
    m2_to_m3_comparisons = nested_comparisons[
        nested_comparisons["comparison"].str.startswith("M2→M3_")
    ].copy()

    if m2_to_m3_comparisons.empty:
        logger.warning("No M2→M3_k comparisons found.")
        return pd.DataFrame()

    # Extract k value from comparison
    def extract_k(comparison: str) -> int:
        """Extract k value from comparison string like 'M2→M3_4'."""
        if "→M3_" in comparison:
            k_str = comparison.split("→M3_")[1]
            try:
                return int(k_str)
            except ValueError:
                return 0
        return 0

    m2_to_m3_comparisons["k"] = m2_to_m3_comparisons["comparison"].apply(extract_k)

    # For each gene, find best model (largest significant ΔR²)
    best_rows = []

    for gene in m2_to_m3_comparisons["gene"].unique():
        gene_comparisons = m2_to_m3_comparisons[m2_to_m3_comparisons["gene"] == gene].copy()

        # Filter to significant improvements
        significant = gene_comparisons[gene_comparisons[pval_col] < alpha].copy()

        if significant.empty:
            # No significant improvements, use largest ΔR² regardless of significance
            best_row = gene_comparisons.loc[gene_comparisons["delta_r_squared"].idxmax()]
            is_sig = False
        else:
            # Use largest significant ΔR²
            best_row = significant.loc[significant["delta_r_squared"].idxmax()]
            is_sig = True

        k = int(best_row["k"])
        model_to = best_row["model_to"]

        # Get top-k neighbor genes from rankings
        top_k_neighbor_cols = []
        top_k_neighbor_genes = []

        if neighbor_rankings.get(gene):
            # Get top-k neighbor columns
            top_k_neighbor_cols = [col for col, _ in neighbor_rankings[gene][:k]]

            # Extract gene names from column names (neigh_mean_GeneName -> GeneName)
            top_k_neighbor_genes = [col.replace("neigh_mean_", "") for col in top_k_neighbor_cols]

        best_rows.append(
            {
                "gene": gene,
                "best_k": k,
                "best_model": model_to,
                "delta_r_squared": best_row["delta_r_squared"],
                "delta_adj_r_squared": best_row["delta_adj_r_squared"],
                "f_test_pval": best_row[pval_col],
                "is_significant": is_sig,
                "top_k_neighbor_genes": ", ".join(top_k_neighbor_genes),
                "top_k_neighbor_cols": ", ".join(top_k_neighbor_cols),
            }
        )

    best_neighbors_df = pd.DataFrame(best_rows)

    # Sort by delta_r_squared (descending)
    best_neighbors_df = best_neighbors_df.sort_values("delta_r_squared", ascending=False)

    logger.info(
        f"✅ Created best neighbor PIGs table for {len(best_neighbors_df)} genes. "
        f"{best_neighbors_df['is_significant'].sum()} genes have significant improvements."
    )

    return best_neighbors_df
