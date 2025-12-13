from __future__ import annotations

import os
from pathlib import Path

# Standard Library Imports
import warnings

import numpy as np

# Third-Party Library Imports
import pandas as pd
from shapely.geometry import Polygon
import yaml

# Local Module Imports
from src.data.download import download_xenium_dataset
from src.scripts.preprocessing.exploration import (
    summarize_missing_by_column,
    summarize_missing_by_row,
)
from src.scripts.preprocessing.geometry import (
    build_polygons,
    close_polygon_if_needed,
    filter_plaques,
    load_plaques,
    rename_centroid_columns,
)
from src.scripts.preprocessing.loaders import (
    combine_cells_and_expression,
    load_cells_table,
    load_expression_matrix,
)
from src.scripts.preprocessing.normalization import pydeseq2_normalize_global
from src.scripts.preprocessing.qc import (
    filter_cells,
)
from src.utils.logging_utils import logger

warnings.filterwarnings("ignore")

# Get absolute path to this script’s directory
preprocessing_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.dirname(preprocessing_dir)
src_dir = os.path.dirname(scripts_dir)
root_dir = os.path.dirname(src_dir)
config_path = os.path.join(root_dir, "configs/plaque_preprocessing.yaml")

# Load YAML safely
with open(os.path.realpath(config_path)) as f:
    plaque_preprocessing_cfg = yaml.safe_load(f)


def preprocess_cells_and_expression(
    url: str, output_dir: str, PATH_TO_DATA_FOLDER: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Download, preprocess, and normalize Xenium single-cell dataset by combining
    cell metadata with expression data, filtering invalid entries, and generating
    both raw and normalized combined datasets.

    Workflow:
        1. Download the dataset if not already available.
        2. Load the cell metadata (`cells.parquet`) and expression matrix (`cell_feature_matrix.h5`).
        3. Combine the two datasets into a single DataFrame.
        4. Analyze missing values per column and per row.
        5. Identify gene columns and filter cells below the 5th percentile of
           transcript count, cell/nucleus area, or total nonzero genes.
        6. Normalize gene expression using PyDESeq2 and recombine with cell metadata.

    Args:
        url (str):
            The URL pointing to the Xenium dataset Zip archive.
        output_dir (str):
            Directory where the downloaded and extracted dataset will be stored.
        PATH_TO_DATA_FOLDER (str):
            Path to the folder where processed CSV outputs (combined and normalized) will be saved.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, list[str]]:
            - filtered_combined_df: Filtered DataFrame containing valid cells and their raw transcript counts.
            - combined_df_normalized: Normalized combined DataFrame (transcript counts normalized via PyDESeq2).
            - gene_cols: List of gene columns in the expression matrix (i.e. gene names).
    """
    # 1. Download the dataset if not already downloaded
    logger.info("\nDownloading the dataset...\n")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(PATH_TO_DATA_FOLDER, exist_ok=True)
    xenium_path = download_xenium_dataset(url, output_dir)

    # 2. Load the data and combine the cells with the expression matrix
    logger.info("\nLoading the data...\n")
    cells_path = f"{xenium_path}/cells.parquet"
    expr_path = f"{xenium_path}/cell_feature_matrix.h5"
    cells_df = load_cells_table(cells_path)
    logger.info(cells_df.head())
    adata = load_expression_matrix(expr_path, False)

    combined_df = combine_cells_and_expression(
        cells_df, adata, f"{PATH_TO_DATA_FOLDER}/combined.csv"
    )
    logger.info(f"Shape of combined_df: {combined_df.shape}")

    # 3. Investigate the missing values
    logger.info("\nInvestigating the missing values...\n")
    col_stats_df = summarize_missing_by_column(combined_df)
    _, row_aggregates = summarize_missing_by_row(combined_df)
    logger.info(
        f"Found {sum(col_stats_df['pct_missing'] > 0)} columns with missing values. "
        "5 columns with the most missing values:\n"
    )
    logger.info(col_stats_df.head())
    logger.info("\nRow-based statistics on missing values:\n")
    logger.info(row_aggregates)

    # 4. Extract metadata columns and gene columns
    gene_cols = adata.var_names.tolist()
    X = adata.X
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

    # Make sure cell IDs are normalized
    cell_ids = adata.obs_names.astype(str).str.strip("b'").str.replace("'", "")
    n_genes_map = pd.DataFrame({"cell_id": cell_ids, "n_genes": n_genes})

    # Attach n_genes to the combined table
    combined_df = combined_df.merge(n_genes_map, on="cell_id", how="inner")

    # 5. Filter cells below quality thresholds
    logger.info("\nFiltering cells...\n")
    filtered_combined_df = filter_cells(
        combined_df.copy(),
        transcript_q=0.05,
        area_q=0.05,
        nucleus_q=0.05,
        genes_q=0.05,
        gene_cols=gene_cols,
    )
    filtered_combined_df = rename_centroid_columns(filtered_combined_df)
    logger.info(
        f"Original: {cells_df.shape[0]} cells → Filtered: {filtered_combined_df.shape[0]} cells"
    )

    # 6. Normalize data using PyDESeq2
    logger.info("\nNormalizing the data...\n")
    adata_normalized, _ = pydeseq2_normalize_global(adata, min_genes=20, apply_log1p=True)
    combined_df_normalized = combine_cells_and_expression(
        filtered_combined_df[cells_df.columns],
        adata_normalized,
        f"{PATH_TO_DATA_FOLDER}/combined_normalized.csv",
    )

    return filtered_combined_df, combined_df_normalized, gene_cols, combined_df


def preprocess_plaques(PATH_TO_DATA_FOLDER: str | Path) -> tuple[pd.DataFrame, Polygon]:
    """
    Preprocess plaque polygons against a brain boundary and remove small plaques.

    This routine:
      1) Loads plaque vertices and builds valid plaque polygons.
      2) Loads the brain polygon, coerces numeric coords, drops NaNs, and closes it if needed.
      3) Filters plaques using geometric rules (holes/convexity/nesting/overlaps).
      4) Removes small plaques below the 5th percentile of area (computed on the
         post-geometry-filtered set).

    Args:
        PATH_TO_DATA_FOLDER: Directory containing:
            - ``plaque_polygons.csv`` (long-form vertices)
            - ``brain_polygon.csv``   (name,x,y; first two lines skipped)

    Returns:
        pd.DataFrame: The filtered plaques dataframe. Expected to include at least:
            - ``geometry``: polygon geometry for the plaque
            - ``area``: area of the polygon
            - ``removed_small``: bool flag (all False in the returned frame since small ones are dropped)

    Raises:
        FileNotFoundError: If an expected CSV does not exist.
        KeyError: If required columns (e.g., ``area``) are missing after polygon construction.

    Notes:
        - Relies on project-scoped helpers: ``load_plaques``, ``build_polygons``,
          ``close_polygon_if_needed``, and ``filter_plaques``.
        - Uses config constants: ``FILL_HOLES``, ``ENFORCE_CONVEX``, ``REMOVE_NESTED``,
          and ``MERGE_OVERLAPS``.
        - Emits progress via the module-level ``logger``.
    """
    base_path = Path(PATH_TO_DATA_FOLDER)

    # 1) Load the plaque data, build valid polygons
    logger.info("\nLoading the plaque data...\n")
    plaques_csv = base_path / "plaque_polygons.csv"
    if not plaques_csv.exists():
        raise FileNotFoundError(f"Plaque CSV not found: {plaques_csv}")

    plaques_long = load_plaques(plaques_csv)
    plaques_poly = build_polygons(plaques_long)

    # 2) Load the brain polygon and ensure it is valid
    logger.info("\nLoading the brain polygon...\n")
    brain_csv = base_path / "brain_polygon.csv"
    if not brain_csv.exists():
        raise FileNotFoundError(f"Brain CSV not found: {brain_csv}")

    brain = (
        pd.read_csv(
            brain_csv,
            sep=",",
            header=None,
            skiprows=2,
            names=["name", "x", "y"],
            engine="python",
        )
        .assign(
            x=lambda d: pd.to_numeric(d["x"], errors="coerce"),
            y=lambda d: pd.to_numeric(d["y"], errors="coerce"),
        )
        .dropna(subset=["x", "y"])
        .reset_index(drop=True)
    )

    brain_xy = brain[["x", "y"]].to_numpy()
    brain_geom = close_polygon_if_needed(brain_xy)

    logger.info(
        "plaques_long: %s vertices | plaques_poly: %s polygons | brain polygon valid: %s",
        f"{len(plaques_long):,}",
        f"{len(plaques_poly):,}",
        getattr(brain_geom, "is_valid", None),
    )

    # 3) Filter the plaques with geometric rules
    logger.info("\nFiltering the plaques...\n")
    plaques_poly = filter_plaques(
        plaques_poly,
        brain_geom,
        FILL_HOLES=plaque_preprocessing_cfg["fill_holes"],
        ENFORCE_CONVEX=plaque_preprocessing_cfg["enforce_convex"],
        REMOVE_NESTED=plaque_preprocessing_cfg["remove_nested"],
        MERGE_OVERLAPS=plaque_preprocessing_cfg["merge_overlaps"],
    )

    # 4) Filter out small plaques by area percentile (5th)
    if "area" not in plaques_poly.columns:
        raise KeyError("Expected 'area' column in plaques_poly after geometry construction.")

    areas = plaques_poly["area"].to_numpy(dtype=float)
    p5 = float(np.percentile(areas, 5)) if len(areas) > 0 else 0.0

    n_before = len(plaques_poly)
    plaques_poly["removed_small"] = plaques_poly["area"] < p5
    plaques_poly = plaques_poly[~plaques_poly["removed_small"]].copy()
    n_after = len(plaques_poly)

    logger.info("Kept after area cutoff: %s out of %s", f"{n_after:,}", f"{n_before:,}")

    return plaques_poly, brain_geom
