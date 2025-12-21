from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import logging
from typing import Any, Protocol

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


class ModelRunner(Protocol):
    """
    Protocol for model runner callables.

    A conforming function should accept (X_train, X_test, y_train, y_test, gene_cols)
    and return (fitted_model, results_dict_or_series, importance_df).
    """

    def __call__(
        self,
        X_train: np.ndarray,
        X_test: np.ndarray,
        y_train: np.ndarray,
        y_test: np.ndarray,
        gene_cols: Sequence[str],
    ) -> tuple[Any, Mapping[str, Any] | pd.Series | dict[str, Any], pd.DataFrame]: ...


@dataclass
class GeneModelingOutput:
    """Container for key outputs of the gene-distance modeling pipeline."""

    trained_models: dict[str, Any]
    results_df: pd.DataFrame
    importance_df: pd.DataFrame
    top_genes: list[str]
    scaler: StandardScaler
    y_true: pd.Series
    y_pred: np.ndarray
    residuals: pd.Series
    pearson_corr: pd.DataFrame
    spearman_corr: pd.DataFrame
    delta_expression: pd.Series


def run_gene_distance_modeling(
    cells_with_distances: pd.DataFrame,
    gene_cols: Sequence[str],
    model_runners: Mapping[str, ModelRunner],
    *,
    target_col: str = "nearest_plaque_center_dist",
    best_model_name: str = "xgb",
    test_size: float = 0.20,
    random_state: int = 42,
    quantile_threshold: float = 0.90,
    n_top_genes_for_overlay: int = 5,
    top_n_importances_for_plot: int = 20,
    plot_model_performance: Callable[[pd.DataFrame], None] | None = None,
    plot_model_performance_interactive: Callable[[pd.DataFrame], None] | None = None,
    plot_top_gene_importances: Callable[[pd.DataFrame, int], None] | None = None,
    plot_top_gene_importances_interactive: (
        Callable[[pd.DataFrame, int], None] | None
    ) = None,
    plot_spatial_overlay: Callable[[pd.DataFrame, str], None] | None = None,
    plot_multi_gene_signature: (
        Callable[[pd.DataFrame, Sequence[str]], None] | None
    ) = None,
    plot_pred_vs_true: Callable[[pd.Series, np.ndarray, str], None] | None = None,
    plot_residual_hist: Callable[[pd.Series, np.ndarray, str], None] | None = None,
    plot_spatial_residual_map: (
        Callable[[pd.DataFrame, pd.Series, np.ndarray, str], None] | None
    ) = None,
    plot_residual_figure: (
        Callable[[pd.DataFrame, pd.Series, np.ndarray, Sequence[str], Any], None] | None
    ) = None,
    plot_bivariate_resid_dist_spatial: (
        Callable[[pd.DataFrame, pd.Series, np.ndarray, str], None] | None
    ) = None,
    plot_bivariate_resid_distance_spatial_interactive: (
        Callable[[pd.DataFrame, pd.Series, np.ndarray, str], None] | None
    ) = None,
    plot_radius_color_spatial: (
        Callable[[pd.DataFrame, pd.Series, np.ndarray, str], None] | None
    ) = None,
    plot_true_pred_kde: Callable[[pd.Series, np.ndarray, str], None] | None = None,
    plot_true_pred_kde_interactive: (
        Callable[[pd.Series, np.ndarray, str], None] | None
    ) = None,
    plot_residual_vs_distance: (
        Callable[[pd.Series, np.ndarray, str], None] | None
    ) = None,
    plot_residual_vs_distance_interactive: (
        Callable[[pd.Series, np.ndarray, str], None] | None
    ) = None,
    plaques_poly: Any | None = None,
    logger: logging.Logger | None = None,
    dropna: bool = True,
) -> GeneModelingOutput:
    """
    Train several models to predict plaque distance from gene expression, collect results,
    compute residual-based diagnostics, and (optionally) run plotting callbacks.

    Parameters
    ----------
    cells_with_distances : pd.DataFrame
        Full cell-level dataframe containing gene columns and the target column.
    gene_cols : Sequence[str]
        Names of gene expression columns to use as features.
    model_runners : Mapping[str, ModelRunner]
        Dict mapping model names to callables with signature:
            (X_train, X_test, y_train, y_test, gene_cols)
            -> (fitted_model, results_mapping_or_series, importance_df)
        where importance_df contains at minimum columns ["gene", "importance"].
    target_col : str, default "nearest_plaque_center_dist"
        Name of the distance target column.
    best_model_name : str, default "xgb"
        The model key used for full-dataset predictions & residual diagnostics.
    test_size : float, default 0.20
        Proportion of the dataset to include in the test split.
    random_state : int, default 42
        Random seed used for reproducibility in the train/test split.
    quantile_threshold : float, default 0.90
        Quantile cutoff on absolute residuals for "poorly predicted" cells.
    n_top_genes_for_overlay : int, default 5
        Number of top genes (from the selected best model) for overlay/signature plots.
    top_n_importances_for_plot : int, default 20
        Number of top features to show in the importance plot callback, if provided.
    plot_model_performance : Optional[Callable[[pd.DataFrame], None]], optional
        Callback to visualize model performance summary. Receives `results_df`.
    plot_top_gene_importances : Optional[Callable[[pd.DataFrame, int], None]], optional
        Callback to plot importances. Receives (`importance_df`, `top_n`).
    plot_spatial_overlay : Optional[Callable[[pd.DataFrame, str], None]], optional
        Callback to plot a spatial overlay for a single gene. Receives (df, gene_name).
    plot_multi_gene_signature : Optional[Callable[[pd.DataFrame, Sequence[str]], None]], optional
        Callback to plot a multi-gene signature. Receives (df, genes).
    plot_pred_vs_true : Optional[Callable[[pd.Series, np.ndarray, str], None]], optional
        Callback to plot predictions vs. truth. Receives (y_true, y_pred, model_name).
    plot_residual_hist : Optional[Callable[[pd.Series, np.ndarray, str], None]], optional
        Callback to plot residual histogram. Receives (y_true, y_pred, model_name).
    plot_spatial_residual_map : Optional[Callable[[pd.DataFrame, pd.Series, np.ndarray, str], None]], optional
        Callback to plot a spatial residual map. Receives (df, y_true, y_pred, model_name).
    plot_residual_figure : Optional[Callable[[pd.DataFrame, pd.Series, np.ndarray, Sequence[str], Any], None]], optional
        Callback for a composite residual figure. Receives (df, y_true, y_pred, gene_cols, plaques_poly).
    plaques_poly : Any, optional
        Optional geometry or other object passed through to `plot_residual_figure`.
    logger : Optional[logging.Logger], optional
        Logger for informational messages. If None, a module-level logger is used.
    dropna : bool, default True
        If True, drop rows with NaNs across `gene_cols + [target_col]` prior to training.

    Returns
    -------
    GeneModelingOutput
        Dataclass containing trained models, result/importance tables, predictions,
        residual diagnostics, and the fitted scaler.
    """
    log = logger or logging.getLogger(__name__)

    required_cols = list(gene_cols) + [target_col]
    df = cells_with_distances.copy()

    if dropna:
        before = len(df)
        df = df.dropna(subset=required_cols)
        dropped = before - len(df)
        if dropped:
            log.info(f"Dropped {dropped} rows with NaNs in {required_cols}.")

    X = df[gene_cols].to_numpy(copy=True)
    y = df[target_col].to_numpy()

    log.info(f"X shape: {X.shape}")
    log.info(f"y shape: {y.shape}")
    log.info(f"Target NaN count: {int(np.isnan(y).sum())}")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    log.info(f"Train shape: {X_train_scaled.shape}, Test shape: {X_test_scaled.shape}")
    log.info(
        f"Mean of first 5 genes (train, after scaling): "
        f"{np.mean(X_train_scaled[:, :5], axis=0)}"
    )
    log.info(
        f"Std  of first 5 genes (train, after scaling): "
        f"{np.std(X_train_scaled[:, :5], axis=0)}"
    )

    trained_models: dict[str, Any] = {}
    results_list: list[dict[str, Any]] = []
    importance_frames: list[pd.DataFrame] = []

    for name, fn in model_runners.items():
        log.info(f"Running {name.upper()}")
        model, results, imp = fn(
            X_train_scaled, X_test_scaled, y_train, y_test, gene_cols
        )

        if isinstance(results, pd.Series):
            results = results.to_dict()
        else:
            results = dict(results)
        results["model"] = name

        trained_models[name] = model
        results_list.append(results)

        imp = imp.copy()
        imp["model"] = name
        importance_frames.append(imp)

    results_df = pd.DataFrame(results_list)
    importance_df = (
        pd.concat(importance_frames, ignore_index=True)
        if importance_frames
        else pd.DataFrame(columns=["gene", "importance", "model"])
    )

    log.info("Model performance summary created.")
    log.info(results_df.head())

    if plot_model_performance is not None:
        plot_model_performance(results_df)
    if plot_model_performance_interactive is not None:
        plot_model_performance_interactive(results_df)
    if plot_top_gene_importances is not None and not importance_df.empty:
        plot_top_gene_importances(importance_df, top_n_importances_for_plot)
    if plot_top_gene_importances_interactive is not None and not importance_df.empty:
        plot_top_gene_importances_interactive(
            importance_df, top_n=top_n_importances_for_plot
        )

    top_genes: list[str] = []
    if (
        not importance_df.empty
        and "importance" in importance_df.columns
        and "gene" in importance_df.columns
    ):
        mask_best = importance_df["model"].str.lower() == best_model_name.lower()
        imp_best = importance_df.loc[mask_best]
        if not imp_best.empty:
            top_genes = imp_best.nlargest(n_top_genes_for_overlay, "importance")[
                "gene"
            ].tolist()

    if top_genes:
        log.info(f"Top genes for {best_model_name.upper()}: {top_genes}")
        if plot_spatial_overlay is not None:
            plot_spatial_overlay(cells_with_distances, gene=top_genes[0])
        if plot_multi_gene_signature is not None:
            plot_multi_gene_signature(cells_with_distances, top_genes)
    else:
        log.info(
            f"No top genes found for model '{best_model_name}'. "
            "Check importance_df contains ['gene','importance','model']."
        )

    if best_model_name not in trained_models:
        raise KeyError(
            f"best_model_name='{best_model_name}' not found in trained_models: "
            f"{list(trained_models.keys())}"
        )

    X_scaled_all = scaler.transform(cells_with_distances[gene_cols].to_numpy())
    best_model = trained_models[best_model_name]
    y_true = cells_with_distances[target_col]
    y_pred = best_model.predict(X_scaled_all)

    residuals = y_true - y_pred
    abs_resid = residuals.abs()

    if plot_pred_vs_true is not None:
        plot_pred_vs_true(y_true, y_pred, model_name=best_model_name)
    if plot_residual_hist is not None:
        plot_residual_hist(y_true, y_pred, model_name=best_model_name)
    if plot_spatial_residual_map is not None:
        plot_spatial_residual_map(
            cells_with_distances, y_true, y_pred, model_name=best_model_name
        )
    if plot_residual_figure is not None:
        plot_residual_figure(
            cells_with_distances, y_true, y_pred, gene_cols, plaques_poly
        )

    import matplotlib.pyplot as plt

    if plot_bivariate_resid_dist_spatial is not None:
        fig, ax = plt.subplots(1, 1, figsize=(7.5, 6))

        fig.subplots_adjust(left=0.26, right=0.98, top=0.90, bottom=0.08)

        plot_bivariate_resid_dist_spatial(
            cells_with_distances,
            y_true,
            y_pred,
            model_name=best_model_name,
            ax=ax,
        )

        fig.suptitle(
            f"{best_model_name.upper()}: residual structure around plaques", fontsize=14
        )

        if plot_bivariate_resid_distance_spatial_interactive is not None:
            plot_bivariate_resid_distance_spatial_interactive(
                cells_with_distances, y_true, y_pred, model_name=best_model_name
            )

        if plot_true_pred_kde is not None:
            plot_true_pred_kde(y_true, y_pred, model_name=best_model_name)
        if plot_true_pred_kde_interactive is not None:
            plot_true_pred_kde_interactive(y_true, y_pred, model_name=best_model_name)

        if plot_residual_vs_distance is not None:
            plot_residual_vs_distance(y_true, y_pred, model_name=best_model_name)
        if plot_residual_vs_distance_interactive is not None:
            plot_residual_vs_distance_interactive(
                y_true, y_pred, model_name=best_model_name
            )

    pearson_corrs: list[tuple[str, float]] = []
    for g in gene_cols:
        g_series = cells_with_distances[g]
        if np.nanvar(g_series.values) > 0:
            r = np.corrcoef(g_series.values, residuals.values)[0, 1]
            if not np.isnan(r):
                pearson_corrs.append((g, float(r)))

    pearson_df = pd.DataFrame(pearson_corrs, columns=["gene", "pearson_r"])
    if not pearson_df.empty:
        pearson_df["abs_r"] = pearson_df["pearson_r"].abs()
        pearson_df = pearson_df.sort_values("abs_r", ascending=False).reset_index(
            drop=True
        )

    spearman_rows: list[tuple[str, float]] = []
    for g in gene_cols:
        rho, _ = spearmanr(cells_with_distances[g], residuals, nan_policy="omit")
        if not np.isnan(rho):
            spearman_rows.append((g, float(rho)))

    spearman_df = pd.DataFrame(spearman_rows, columns=["gene", "spearman_rho"])
    if not spearman_df.empty:
        spearman_df["abs_rho"] = spearman_df["spearman_rho"].abs()
        spearman_df = spearman_df.sort_values("abs_rho", ascending=False).reset_index(
            drop=True
        )

    thr = abs_resid.quantile(quantile_threshold)
    poor_mask = abs_resid >= thr
    good_mask = abs_resid < thr

    poor_means = cells_with_distances.loc[poor_mask, gene_cols].mean(numeric_only=True)
    good_means = cells_with_distances.loc[good_mask, gene_cols].mean(numeric_only=True)
    delta_expression = (poor_means - good_means).sort_values(ascending=False)
    log.info(delta_expression.head(15))

    log.info("Gene modeling pipeline completed.")

    return GeneModelingOutput(
        trained_models=trained_models,
        results_df=results_df,
        importance_df=importance_df,
        top_genes=top_genes,
        scaler=scaler,
        y_true=y_true,
        y_pred=y_pred,
        residuals=residuals,
        pearson_corr=pearson_df,
        spearman_corr=spearman_df,
        delta_expression=delta_expression,
    )
