from __future__ import annotations

import csv
import pathlib

import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon, shape
from skimage.measure import ransac
from skimage.transform import AffineTransform, SimilarityTransform

from src.utils.logging_utils import logger


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    """Root-mean-square error between two (N,2) arrays."""
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def points_from_keypoints_df(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract source (IF/alignment) and destination (morphology/fixed) arrays from keypoints.

    Returns
    -------
    src_if : (N,2) float array  -- alignment (IF) pixel coordinates
    dst_morph : (N,2) float array -- fixed (morphology) pixel coordinates
    """
    cols_needed = {"fixedX", "fixedY", "alignmentX", "alignmentY"}
    missing = cols_needed.difference(df.columns)
    if missing:
        logger.error("Keypoints CSV missing columns: %s", missing)
        raise ValueError(f"Keypoints CSV missing columns: {missing}")

    src_if = df[["alignmentX", "alignmentY"]].to_numpy(float)
    dst_morph = df[["fixedX", "fixedY"]].to_numpy(float)
    logger.debug("Loaded %d keypoints", len(df))
    return src_if, dst_morph


def fit_transform_with_ransac(
    src_if_px: np.ndarray,
    dst_morph_px: np.ndarray,
    model_class: type[SimilarityTransform | AffineTransform],
    min_samples: int,
    residual_threshold: float,
    max_trials: int,
):
    """
    Robustly estimate a geometric transform (Similarity or Affine) from IF -> Morphology pixels.

    Parameters
    ----------
    src_if_px : (N,2) float
        IF (alignment) pixel coordinates.
        Here N is the number of keypoints used for alignment.
        2 columns correspond to x/y pixel coordinates.
        IF pixels are source (src) in this case: we want to begin with the IF image (source)
        and then map it into the morphology image grid (destination, dst).
    dst_morph_px : (N,2) float
        Morphology (fixed) pixel coordinates.
    model_class : class
        SimilarityTransform or AffineTransform (from skimage.transform).
    Returns
    -------
    model : fitted transform with .params (3x3)
    inliers :
        boolean mask of inlier matches
        It is computed in the following fashion:
        After fitting a tentative model from a random subset,
        RANSAC computes each pair's residual (how far the transformed source point lands
        from the target). Pairs with residual below the threshold are marked True (inliers);
        the rest are False (outliers)
    rmse_all : float, root mean squared error (RMSE) over all pairs using model
    """
    logger.debug(
        "Fitting %s with RANSAC (min_samples=%d, residual_threshold=%.3f, max_trials=%d)",
        model_class.__name__,
        min_samples,
        residual_threshold,
        max_trials,
    )
    model_robust, inliers = ransac(
        (src_if_px, dst_morph_px),
        model_class,
        min_samples=min_samples,
        residual_threshold=residual_threshold,
        max_trials=max_trials,
    )
    # Evaluate RMSE over *all* matches (not just inliers), to decide upgrade
    pred_all = model_robust(src_if_px)
    rmse_all = _rmse(pred_all, dst_morph_px)
    logger.info(
        "RANSAC %s → inliers: %d/%d | RMSE(px): %.3f",
        model_class.__name__,
        int(inliers.sum()),
        len(inliers),
        rmse_all,
    )
    return model_robust, inliers, rmse_all


def build_shapely_xy_transform(matrix_3x3: np.ndarray):
    """
    Return a function usable with shapely.ops.transform that applies the 3x3 transform.

    The function takes (x, y) and returns transformed (x', y').

    Parameters
    ----------
    matrix_3x3 : np.ndarray
        3x3 homogeneous transform mapping IF pixels → morphology pixels.

    Returns
    -------
    f : callable
        f(x, y, z=None) -> (x', y')  (ignores z)
    """

    def _f(x, y, z=None):
        xy1 = np.vstack([np.asarray(x), np.asarray(y), np.ones_like(x)])
        out = matrix_3x3 @ xy1
        return out[0], out[1]  # drop homogeneous coordinate

    return _f


def microns_to_pixels(xy_um: np.ndarray, pixel_size_um: float) -> np.ndarray:
    """Convert microns → pixels for Xenium morphology image."""
    return xy_um / pixel_size_um


def pixels_to_microns(xy_px: np.ndarray, pixel_size_um: float) -> np.ndarray:
    """Convert pixels → microns for Xenium morphology image."""
    return xy_px * pixel_size_um


def polygon_name(base: str, idx: int) -> str:
    """Generate a human-friendly, unique selection name."""
    return f"{base} {idx+1}"


def iter_exterior_polygons(geom_obj) -> list[Polygon]:
    """Yield exterior polygons from Polygon/MultiPolygon, ignore holes."""
    if isinstance(geom_obj, Polygon):
        yield geom_obj
    elif isinstance(geom_obj, MultiPolygon):
        yield from geom_obj.geoms
    else:
        # unsupported geometry types are skipped
        return


def load_polygons_from_geojson(gj: dict) -> list[tuple[str, Polygon]]:
    """
    Load exterior polygons and derive a selection name for each.

    Returns
    -------
    list of (name, polygon)
    """
    selections: list[tuple[str, Polygon]] = []
    feat_list = gj.get("features", [])
    counter = 0

    for feat in feat_list:
        base_name = "Selection"  # this format is required for import into Xenium Explorer

        geom = shape(feat.get("geometry"))
        for poly in iter_exterior_polygons(geom):
            # Ensure polygon is valid and has enough points
            if not isinstance(poly, Polygon) or len(poly.exterior.coords) < 3:
                continue
            sel_name = polygon_name(base_name, counter)
            selections.append((sel_name, poly))
            counter += 1

    return selections


def write_selections_csv(
    selections: list[tuple[str, Polygon]],
    out_csv: pathlib.Path,
    decimals: int = 3,
    top_n: int | None = None,
):
    """
    Write selections to CSV matching Xenium Explorer's importable format.

    Format (multiple):
      #Selection names: Sel 1, Sel 2, ...
      #Areas (µm^2): 123.45, 678.90
      Selection,X,Y
      Sel 1, x, y
      Sel 1, x, y
      ...
      Sel 2, x, y
      ...

    Format (single):
      #Selection name: Sel 1
      #Area (µm^2): 123.45
      Selection,X,Y
      Sel 1, x, y
      ...
    """
    if not selections:
        raise ValueError("No valid polygons found in GeoJSON.")

    names = [name for name, _ in selections]
    areas = [poly.area for _, poly in selections]  # microns squared (coords are microns)
    if top_n is not None:
        # Let us sort and pick top_n plaques with the largest area
        sorted_indices = np.argsort(areas)[::-1]
        top_indices = sorted_indices[:top_n]
        names = [names[i].split(" ")[0] + f" {idx+1}" for idx, i in enumerate(top_indices)]
        areas = [areas[i] for i in top_indices]
        selections = [(names[idx], selections[i][1]) for idx, i in enumerate(top_indices)]

    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        # --- Write the header lines as plain text (not via csv.writer) ---
        if len(selections) == 1:
            f.write(f"#Selection name: {names[0]}\n")
            f.write(f"#Area (µm^2): {areas[0]:.2f}\n")
        else:
            f.write("#Selection names: " + ", ".join(names) + "\n")
            f.write("#Areas (µm^2): " + ", ".join(f"{a:.2f}" for a in areas) + "\n")

        # --- Now write the table using csv.writer ---
        writer = csv.writer(f, lineterminator="\n")  # prevents extra CRLF on Windows
        writer.writerow(["Selection", "X", "Y"])

        fmt = f"{{:.{decimals}f}}"
        for name, poly in selections:
            # exterior coords are typically closed (first point repeated last) - keep as-is
            for x, y in poly.exterior.coords:
                writer.writerow([name, fmt.format(float(x)), fmt.format(float(y))])
