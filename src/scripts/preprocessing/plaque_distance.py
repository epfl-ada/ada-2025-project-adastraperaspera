"""
Compute per-cell distances to the nearest amyloid-beta plaque polygon.
"""

from __future__ import annotations

import geopandas as gpd
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from shapely.strtree import STRtree

from src.utils.logging_utils import logger


def load_plaque_polygons(csv_path: str) -> gpd.GeoDataFrame:
    """
    Load plaque polygons from a QuPath-style CSV file with columns Selection, X, Y.

    Handles files that include header comments and groups coordinates by 'Selection'.

    Args:
        csv_path (str): Path to plaque_polygons.csv

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with Polygon geometries and selection IDs.
    """
    import io

    logger.info(f"Parsing plaque polygons from {csv_path}")

    with open(csv_path, encoding="utf-8") as f:
        lines = [ln for ln in f.readlines() if not ln.startswith("#") and ln.strip()]

    df = pd.read_csv(io.StringIO("".join(lines)))

    if not {"Selection", "X", "Y"}.issubset(df.columns):
        raise ValueError("CSV must have 'Selection', 'X', 'Y' columns after headers.")

    polygons = []
    names = []
    for sel, group in df.groupby("Selection"):
        coords = list(zip(group["X"], group["Y"], strict=False))

        if len(coords) >= 3:
            try:
                poly = Polygon(coords)
                if poly.is_valid:
                    polygons.append(poly)
                    names.append(sel)
            except Exception as e:
                logger.warning(f"Failed to create polygon for {sel}: {e}")
                continue

    gdf = gpd.GeoDataFrame({"selection": names, "geometry": polygons})
    gdf.set_crs(epsg=4326, inplace=True, allow_override=True)
    logger.info(f" Loaded {len(gdf)} plaque polygons.")
    return gdf


def compute_cell_to_plaque_distance(
    cells_df: pd.DataFrame,
    plaques_gdf: gpd.GeoDataFrame,
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
) -> pd.DataFrame:
    """
    Compute distance from each cell centroid to the nearest plaque polygon.
    """
    if plaques_gdf.empty:
        logger.warning("Plaque GeoDataFrame is empty; distances set to NaN.")
        cells_df["distance_to_plaque"] = float("nan")
        return cells_df

    logger.info(
        f"Computing distances for {len(cells_df)} cells → {len(plaques_gdf)} plaques"
    )

    cell_points = gpd.GeoSeries(
        gpd.points_from_xy(cells_df[x_col], cells_df[y_col]), crs=plaques_gdf.crs
    )

    tree = STRtree(plaques_gdf.geometry.values)

    distances = []
    for pt in cell_points:
        nearest_geom = tree.geometries[tree.nearest(pt)]
        distances.append(pt.distance(nearest_geom))

    cells_out = cells_df.copy()
    cells_out["distance_to_plaque"] = distances
    logger.info(" Distance computation complete.")
    return cells_out
