from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.ops import nearest_points, unary_union
from shapely.strtree import STRtree
from tqdm import tqdm

from src.utils.logging_utils import logger


def is_poly(g):
    return isinstance(g, (Polygon, MultiPolygon))


def fill_holes(g):
    if isinstance(g, Polygon):
        return Polygon(g.exterior)
    if isinstance(g, MultiPolygon):
        return MultiPolygon([Polygon(p.exterior) for p in g.geoms if not p.exterior.is_empty])
    return g


def convexify(g):
    return g.convex_hull


# Filter plaques
def filter_plaques(
    plaques: pd.DataFrame,
    brain_geom: Polygon,
    FILL_HOLES: bool = True,
    ENFORCE_CONVEX: bool = True,
    REMOVE_NESTED: bool = True,
    MERGE_OVERLAPS: bool = True,
):
    plaques_poly = plaques.copy()
    plaques_poly["geometry"] = plaques_poly["geometry"].apply(repair)
    if FILL_HOLES:
        plaques_poly["geometry"] = plaques_poly["geometry"].map(fill_holes)
    if ENFORCE_CONVEX:
        plaques_poly["geometry"] = plaques_poly["geometry"].map(convexify)
    if REMOVE_NESTED and len(plaques_poly):
        geoms = plaques_poly["geometry"].tolist()
        tree = STRtree(geoms)
        nested = set()
        try:
            src_idx, tree_idx = tree.query(geoms, predicate="contains", return_indices=True)
            for i, j in zip(src_idx, tree_idx, strict=False):
                if i == j:
                    continue  # skip self-pairs
                nested.add(j)  # j is the contained geometry (should be dropped)
        except TypeError:
            # Older Shapely versions
            for i, g in enumerate(geoms):
                for cand in tree.query(g):  # returns geometries
                    if cand is g:
                        continue
                    try:
                        if cand.contains(g):
                            nested.add(i)  # i is the contained geometry (should be dropped)
                            break
                    except Exception:
                        pass
        if nested:
            keep_idx = [k for k in range(len(geoms)) if k not in nested]
            plaques_poly = plaques_poly.iloc[keep_idx].copy()
    parts = []
    if MERGE_OVERLAPS and len(plaques_poly):
        # Remove overlapping plaques (if any)
        # Overlapping or touching polygons get fused. Separate islands remain separate but are returned together.
        u = unary_union(plaques_poly["geometry"].tolist())
        if isinstance(u, Polygon):
            parts = [u]
        elif isinstance(u, MultiPolygon):
            parts = list(u.geoms)
        else:
            parts = []
    else:
        parts = plaques_poly["geometry"].tolist()

    plaques_poly = pd.DataFrame({"geometry": parts})
    if len(plaques_poly) == 0:
        plaques_poly["plaque_id"] = []
    else:
        plaques_poly["plaque_id"] = range(1, len(plaques_poly) + 1)

    plaques_poly = plaques_poly[
        plaques_poly["geometry"].map(is_poly)
        & plaques_poly["geometry"].map(lambda g: g.is_valid and g.area > 0)
    ].copy()

    plaques_poly["isPolygon"] = plaques_poly["geometry"].map(lambda g: isinstance(g, Polygon))
    plaques_poly["isMultiPolygon"] = plaques_poly["geometry"].map(
        lambda g: isinstance(g, MultiPolygon)
    )
    plaques_poly["is_valid"] = plaques_poly["geometry"].map(lambda g: getattr(g, "is_valid", False))
    plaques_poly["has_holes"] = plaques_poly["geometry"].map(has_holes)
    plaques_poly["is_convex"] = plaques_poly["geometry"].map(
        lambda g: hasattr(g, "convex_hull") and g.equals(g.convex_hull)
    )
    plaques_poly["area"] = plaques_poly["geometry"].map(lambda g: getattr(g, "area", 0.0))
    plaques_poly["centroid_x"] = plaques_poly["geometry"].map(
        lambda g: getattr(getattr(g, "centroid", None), "x", np.nan)
    )
    plaques_poly["centroid_y"] = plaques_poly["geometry"].map(
        lambda g: getattr(getattr(g, "centroid", None), "y", np.nan)
    )
    plaques_poly["in_brain"] = plaques_poly["geometry"].map(
        lambda g, bg=brain_geom: (g is not None) and g.intersects(bg)
    )

    plaques_poly = plaques_poly[
        plaques_poly["is_valid"] & (plaques_poly["area"] > 0) & plaques_poly["in_brain"]
    ].copy()
    plaques_poly["plaque_id"] = range(1, len(plaques_poly) + 1)

    logger.info(
        {
            "n_final": int(len(plaques_poly)),
            "polygons": int(plaques_poly["isPolygon"].sum()),
            "multipolygons": int(plaques_poly["isMultiPolygon"].sum()),
            "holes_present": int(plaques_poly["has_holes"].sum()),
            "convex": int(plaques_poly["is_convex"].sum()),
        }
    )
    return plaques_poly


def compute_cell_to_plaque_distances(
    cells_df: pd.DataFrame,
    plaques_df: pd.DataFrame,
) -> pd.DataFrame:
    pg = plaques_df["geometry"].tolist()
    pid = plaques_df["plaque_id"].to_numpy()
    parea = plaques_df["area"].to_numpy()
    pcent = np.array([[g.centroid.x, g.centroid.y] for g in pg])
    tree = STRtree(pg)
    results = []
    for _, row in tqdm(
        cells_df[["x_centroid", "y_centroid"]].iterrows(),
        total=len(cells_df),
        desc="Computing distances",
        ncols=100,
    ):
        x, y = float(row.x_centroid), float(row.y_centroid)
        p = Point(x, y)
        # Get nearest geometry index directly
        j = tree.nearest(p)  # returns integer index of nearest geometry
        nearest_geom = pg[j]
        # Is this cell occurring inside aplaque?
        inside = nearest_geom.contains(p)
        # Closest boundary point and distances
        np1, np2 = nearest_points(p, nearest_geom)
        dist_boundary = p.distance(nearest_geom)
        # Distance to plaque centroid
        cx, cy = pcent[j]
        dist_center = np.hypot(x - cx, y - cy)
        results.append(
            (
                inside,
                pid[j],
                dist_boundary,
                dist_center,
                np2.x,
                np2.y,
                parea[j],
            )
        )

    cols = [
        "inside_any_plaque",
        "nearest_plaque_id",
        "nearest_plaque_dist",
        "nearest_plaque_center_dist",
        "closest_point_x",
        "closest_point_y",
        "nearest_plaque_area",
    ]
    dist_df = pd.DataFrame(results, columns=cols, index=cells_df.index)
    cells_aug = pd.concat([cells_df, dist_df], axis=1)

    logger.info({"cells": int(len(cells_aug)), "plaques_used": int(len(plaques_df))})
    return cells_aug


def rename_centroid_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Helper function to rename the centroid columns.

    Args:
        df: The dataframe to rename the centroid columns of.

    Returns:
        The dataframe with the renamed centroid columns.
    """
    cand_x = [c for c in ["x_centroid", "x", "x_center", "centroid_x"] if c in df.columns]
    cand_y = [c for c in ["y_centroid", "y", "y_center", "centroid_y"] if c in df.columns]
    if not cand_x or not cand_y:
        raise ValueError(
            "Missing centroid columns. Provide x_centroid/y_centroid (or synonyms: x,y / x_center,y_center / centroid_x,centroid_y)."
        )
    if cand_x[0] != "x_centroid":
        df = df.rename(columns={cand_x[0]: "x_centroid"})
    if cand_y[0] != "y_centroid":
        df = df.rename(columns={cand_y[0]: "y_centroid"})
    return df


def repair(g: Polygon) -> Polygon:
    """
    If needed, reconstructs the polygon from its boundary at distance 0
    During reconstruction, GEOS cleans up self-touching rings, bowties, duplicate vertices, and similar artifacts

    Args:
        g: A Polygon object

    Returns:
        A Polygon object
    """
    try:
        if g.is_valid:
            return g
        else:
            return g.buffer(0)
    except Exception:
        return g


def has_holes(g: Polygon | MultiPolygon) -> bool:
    """
    Checks if a polygon or multi-polygon has holes

    Args:
        g: A Polygon or MultiPolygon object

    Returns:
        A boolean value
    """
    if isinstance(g, Polygon):
        return len(g.interiors) > 0
    if isinstance(g, MultiPolygon):
        return any(len(p.interiors) > 0 for p in g.geoms)
    return False


def load_plaques(
    plaque_path: Path,
) -> pd.DataFrame:
    """
    Loads the plaques from a CSV file.

    Args:
        plaque_path: Path to the CSV file containing the plaques.

    Returns:
        A DataFrame containing the plaques.
    """
    raw = plaque_path.read_text(encoding="utf-8").splitlines()
    areas_meta = []
    # Load the areas metadata from the first line if it exists
    if raw and raw[0].lstrip().startswith("#Areas"):
        parts = [p.strip() for p in raw[0].split(",")]
        areas_meta = [float(x) for x in parts[1:] if x.replace(".", "", 1).isdigit()]
    # Load the plaques from the remaining lines
    rows = []
    for line in raw[2:]:
        if not line.strip():
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 3:
            continue
        name, xs, ys = parts[0], parts[1], parts[2]
        try:
            x = float(xs)
            y = float(ys)
        except Exception:
            continue
        rows.append((name, x, y))

    plaques_long = pd.DataFrame(rows, columns=["plaque_name", "x", "y"])
    plaques_long["plaque_id"] = plaques_long["plaque_name"].str.extract(r"(\d+)").astype(int)
    if areas_meta:
        area_map = {
            i + 1: areas_meta[i]
            for i in range(min(len(areas_meta), plaques_long["plaque_id"].max()))
        }
        plaques_long["area_meta"] = plaques_long["plaque_id"].map(area_map)
    return plaques_long


def close_polygon_if_needed(
    xy: np.ndarray,
) -> Polygon | None:
    """
    Closes the polygon if needed.

    Args:
        xy: numpy array of shape (n, 2) containing the polygon vertices.

    Returns:
        A Polygon object.
    """
    if len(xy) < 3:
        return None
    # If the first and last points are not exactly identical, append the first point to the end
    if not (xy[0] == xy[-1]).all():
        # Ensure that the polygon is explicitly closed
        xy = np.vstack([xy, xy[0]])
    return Polygon(xy)


def build_polygons(
    df: pd.DataFrame,
    id_col: str = "plaque_id",
) -> pd.DataFrame:
    """
    Builds the polygons from a DataFrame.

    Args:
        df: DataFrame containing the polygons.
        id_col: Column name containing the polygon IDs.

    Returns:
        A DataFrame containing the polygons.
    """
    geoms, ids = [], []
    for pid, grp in df.groupby(id_col, sort=True):
        xy = grp[["x", "y"]].to_numpy()
        geom = close_polygon_if_needed(xy)
        if geom is not None:
            geoms.append(geom)
            ids.append(pid)
    return pd.DataFrame({"plaque_id": ids, "geometry": geoms})
