"""
aligner.py

Transform plaque polygons drawn on the ImmunoFluorescence (IF) image
(derived from QuPath GeoJSON, with pixel coordinates) into the Xenium
morphology coordinate frame (microns), using keypoints exported from
Xenium Explorer.

Inputs
------
- keypoints CSV: columns [fixedX,fixedY,alignmentX,alignmentY]
  * fixed   = morphology image pixels (DAPI; "base")
  * alignment = IF image pixels (to be transformed)
- plaque GeoJSON:
    QuPath-exported polygons around the amyloid
    beta plaques on the IF image (pixel coords)
- cells.csv.gz:
    Xenium cell centroids, in microns
    (x_centroid,y_centroid)

Outputs
-------
- plaque_polygons_in_morphology_microns.geojson
    (geometry in microns)
- cell_plaque_distances.csv.gz
    (cell_id, x/y in µm, distance_to_nearest_plaque_µm,
    inside_any_plaque)

Notes
-----
- Xenium morphology pixel size: 0.2125 µm / pixel
  (cf. 10x docs: cells/transcripts in microns; origin = top-left).
  [Source](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html)
- We estimate a SimilarityTransform via RANSAC; if RMS error is high,
  we retry with an AffineTransform.
"""

from __future__ import annotations

import logging

import geopandas as gpd
import pandas as pd
from shapely.ops import transform as shp_transform
from skimage.transform import AffineTransform, SimilarityTransform

from src.scripts.plaque_alignment.config import AppCfg
from src.scripts.plaque_alignment.utils import (
    build_shapely_xy_transform,
    fit_transform_with_ransac,
    load_polygons_from_geojson,
    points_from_keypoints_df,
    write_selections_csv,
)
from src.utils.logging_utils import setup_logging

logger = logging.getLogger(__name__)


def run_align_and_export(cfg: AppCfg) -> int:
    """
    Run the end-to-end alignment and export workflow.

    Parameters
    ----------
    cfg : AppCfg
        Fully-loaded application configuration (absolute paths, validated params).

    Returns
    -------
    int
        0 on success, non-zero on failure.
    """
    setup_logging(
        cfg.paths,  # uses cfg.paths.logs_dir internally
        level_str=str(cfg.logging.get("level", "INFO")),
        tzname=str(cfg.logging.get("timezone", "Europe/Zurich")),
    )
    logger.info("=== Xenium plaque alignment run started ===")
    logger.info(
        "Inputs: keypoints=%s | plaques=%s",
        cfg.paths.keypoints_csv,
        cfg.paths.plaque_geojson,
    )
    logger.info("Output (Selections CSV): %s", cfg.paths.out_xenium_format)
    # 1) Read keypoints and fit transform IF→Morphology
    kp = pd.read_csv(cfg.paths.keypoints_csv)
    src_if_px, dst_morph_px = points_from_keypoints_df(kp)

    # First try Similarity + RANSAC
    sim_model, sim_inliers, sim_rmse = fit_transform_with_ransac(
        src_if_px,
        dst_morph_px,
        model_class=SimilarityTransform,
        min_samples=cfg.params.ransac.min_samples,
        residual_threshold=cfg.params.ransac.residual_threshold_px,
        max_trials=cfg.params.ransac.max_trials,
    )

    # If residual too large, upgrade to Affine
    if sim_rmse > cfg.params.upgrade_to_affine_rmse_px:
        logger.warning(
            "Similarity RMSE (%.3f px) > threshold (%.3f px). Trying Affine.",
            sim_rmse,
            cfg.params.upgrade_to_affine_rmse_px,
        )
        aff_model, aff_inliers, aff_rmse = fit_transform_with_ransac(
            src_if_px,
            dst_morph_px,
            model_class=AffineTransform,
            min_samples=cfg.params.ransac.min_samples,
            residual_threshold=cfg.params.ransac.residual_threshold_px,
            max_trials=cfg.params.ransac.max_trials,
        )
        if aff_rmse < sim_rmse:
            model, inliers, rmse = aff_model, aff_inliers, aff_rmse
            model_name = "AffineTransform (RANSAC)"
        else:
            model, inliers, rmse = sim_model, sim_inliers, sim_rmse
            model_name = "SimilarityTransform (RANSAC)"
    else:
        model, inliers, rmse = sim_model, sim_inliers, sim_rmse
        model_name = "SimilarityTransform (RANSAC)"

    logger.info(
        "Selected model: %s | RMSE(px)=%.3f | Inliers=%d/%d",
        model_name,
        rmse,
        int(inliers.sum()),
        len(inliers),
    )
    logger.debug("Transform matrix (3x3):\n%s", model.params)

    # 2) Transform plaque polygons (IF pixels) → morphology pixels → microns
    gdf = gpd.read_file(cfg.paths.plaque_geojson)
    logger.info("Loaded %d plaque geometries from GeoJSON", len(gdf))

    # (a) IF px → morphology px
    f_px_to_px = build_shapely_xy_transform(model.params)
    gdf_morph_px = gdf.copy()
    gdf_morph_px["geometry"] = gdf_morph_px["geometry"].apply(
        lambda geom: shp_transform(f_px_to_px, geom)
    )

    # (b) morphology px → microns (Xenium units)
    def _px_to_um_geom(geom):
        # Return arrays in the signature expected by shapely.ops.transform: (x', y')
        scale = cfg.params.xenium_pixel_size_um
        return shp_transform(lambda x, y, z=None: (x * scale, y * scale), geom)

    gdf_morph_um = gdf_morph_px.copy()
    gdf_morph_um["geometry"] = gdf_morph_um["geometry"].apply(_px_to_um_geom)
    # Add metadata to properties
    gdf_morph_um["coord_units"] = "micron"
    gdf_morph_um["transform_model"] = model_name
    filtered_selections_around_plaques = load_polygons_from_geojson(gdf_morph_um.__geo_interface__)
    selections_cfg = cfg.params.selections
    write_selections_csv(
        selections=filtered_selections_around_plaques,
        out_csv=cfg.paths.out_xenium_format,
        decimals=selections_cfg.get("decimals"),
        top_n=selections_cfg.get("top_n"),
    )
    logger.info(
        "[OK] Wrote %d selection(s) to: %s",
        len(filtered_selections_around_plaques),
        cfg.paths.out_xenium_format,
    )
    return 0
