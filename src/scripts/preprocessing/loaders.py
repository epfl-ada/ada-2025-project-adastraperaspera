"""
Functions for loading Xenium data tables and expression matrices,
and merging expression with per-cell metadata.
"""

from __future__ import annotations

import os

import pandas as pd
import scanpy as sc

from src.utils.logging_utils import logger


def load_cells_table(path: str) -> pd.DataFrame:
    """
    Load the Xenium cells table (.parquet).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cells table not found: {path}")

    logger.info(f"Loading cells table from {path}")
    df = pd.read_parquet(path)
    logger.info(f"Loaded {df.shape[0]} cells with {df.shape[1]} columns")
    return df


def load_expression_matrix(path: str, normalize: bool = True) -> sc.AnnData:
    """
    Load the Xenium expression matrix (.h5) into an AnnData object,
    optionally normalized with log1p.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Expression matrix not found: {path}")

    logger.info(f"Loading expression matrix from {path}")
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()

    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    logger.info(f"Loaded expression matrix: {adata.n_obs} cells * {adata.n_vars} genes")
    return adata


def expression_to_dataframe(adata: sc.AnnData) -> pd.DataFrame:
    """
    Convert AnnData expression matrix to a pandas DataFrame with clean cell IDs.
    Supports both sparse and dense matrices.
    """
    logger.info("Converting AnnData to DataFrame...")

    x = adata.X
    if hasattr(x, "tocsc"):  # sparse matrix
        expr_df = pd.DataFrame.sparse.from_spmatrix(
            x, index=adata.obs_names.astype(str), columns=adata.var_names
        )
    else:  # dense numpy array
        expr_df = pd.DataFrame(x, index=adata.obs_names.astype(str), columns=adata.var_names)

    expr_df.index = expr_df.index.astype(str).str.strip("b'").str.replace("'", "")
    expr_df.index.name = "cell_id"

    logger.info(f"Expression DataFrame: {expr_df.shape[0]} * {expr_df.shape[1]}")
    return expr_df


def merge_expression_with_cells(cells_df: pd.DataFrame, adata: sc.AnnData) -> pd.DataFrame:
    """
    Merge filtered cell metadata (with distances) and gene expression values.

    Args:
        cells_df (pd.DataFrame): e.g. cells_with_dist
        adata (sc.AnnData): expression matrix

    Returns:
        pd.DataFrame: combined morphology + expression dataframe.
    """
    expr_df = expression_to_dataframe(adata)
    cells_df = cells_df.copy()
    cells_df["cell_id"] = cells_df["cell_id"].astype(str).str.strip("b'").str.replace("'", "")

    common_ids = set(cells_df["cell_id"]) & set(expr_df.index)
    logger.info(f"Merging {len(common_ids)} common cell IDs")

    merged = cells_df[cells_df["cell_id"].isin(common_ids)].merge(
        expr_df.loc[list(common_ids)], on="cell_id", how="inner"
    )

    logger.info(f"Merged dataframe shape: {merged.shape}")
    return merged
