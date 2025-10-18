from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("yaml")  # config/logging may import yaml transitively
skimage = pytest.importorskip("skimage")
shapely = pytest.importorskip("shapely")

from shapely.geometry import Polygon  # noqa: E402
from skimage.transform import SimilarityTransform  # noqa: E402

from src.scripts.plaque_alignment import utils as plaque_utils  # noqa: E402


def test_points_from_keypoints_df_ok() -> None:
    df = pd.DataFrame(
        {
            "fixedX": [0.0, 10.0, 20.0],
            "fixedY": [0.0, 10.0, 20.0],
            "alignmentX": [1.0, 2.0, 3.0],
            "alignmentY": [4.0, 5.0, 6.0],
        }
    )
    src_if, dst_morph = plaque_utils.points_from_keypoints_df(df)
    assert src_if.shape == (3, 2)
    assert dst_morph.shape == (3, 2)
    assert np.allclose(src_if, np.array([[1, 4], [2, 5], [3, 6]], dtype=float))
    assert np.allclose(dst_morph, np.array([[0, 0], [10, 10], [20, 20]], dtype=float))


def test_points_from_keypoints_df_missing_col_raises() -> None:
    df = pd.DataFrame({"fixedX": [0], "fixedY": [0], "alignmentX": [0]})
    with pytest.raises(ValueError, match="Keypoints CSV missing columns"):
        plaque_utils.points_from_keypoints_df(df)


def test_fit_transform_with_ransac_similarity_recovers_transform() -> None:
    rng = np.random.default_rng(42)
    n = 50
    src = rng.uniform(low=-50, high=50, size=(n, 2))

    true_tform = SimilarityTransform(scale=1.25, rotation=0.15, translation=(5.0, -3.0))
    dst = true_tform(src)
    noise = rng.normal(scale=0.05, size=dst.shape)
    dst_noisy = dst + noise

    model, inliers, rmse = plaque_utils.fit_transform_with_ransac(
        src_if_px=src,
        dst_morph_px=dst_noisy,
        model_class=SimilarityTransform,
        min_samples=3,
        residual_threshold=0.2,
        max_trials=200,
    )
    assert inliers.sum() >= n * 0.8
    # Check that model maps a sample point close to true
    test_pt = np.array([[10.0, -7.0]])
    pred = model(test_pt)
    true = true_tform(test_pt)
    assert np.allclose(pred, true, atol=0.2)
    assert rmse < 0.2


def test_build_shapely_xy_transform_applies_matrix() -> None:
    # 2x scale and translation (tx, ty) = (1, -2)
    mat = np.array([[2.0, 0.0, 1.0], [0.0, 2.0, -2.0], [0.0, 0.0, 1.0]])
    f = plaque_utils.build_shapely_xy_transform(mat)

    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    from shapely.ops import transform as shp_transform

    out = shp_transform(f, poly)
    xs, ys = zip(*list(out.exterior.coords)[:4], strict=False)
    assert np.allclose(xs, [1.0, 3.0, 3.0, 1.0], atol=1e-6)
    assert np.allclose(ys, [-2.0, -2.0, 0.0, 0.0], atol=1e-6)


def test_unit_conversions() -> None:
    arr = np.array([[2.0, 4.0], [10.0, 20.0]], dtype=float)
    pix = plaque_utils.microns_to_pixels(arr, pixel_size_um=0.5)
    assert np.allclose(pix, np.array([[4.0, 8.0], [20.0, 40.0]]))
    um = plaque_utils.pixels_to_microns(pix, pixel_size_um=0.5)
    assert np.allclose(um, arr)


def test_load_and_write_selections_csv(tmp_path: Path) -> None:
    # Build GeoJSON dict with two polygons of differing areas
    gj = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[(0, 0), (2, 0), (2, 1), (0, 1), (0, 0)]],
                },
                "properties": {},
            },
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]],
                },
                "properties": {},
            },
        ],
    }

    sels = plaque_utils.load_polygons_from_geojson(gj)
    assert len(sels) == 2
    assert all(isinstance(name, str) for name, _ in sels)

    out_csv = tmp_path / "selections.csv"
    plaque_utils.write_selections_csv(sels, out_csv, decimals=2, top_n=1)

    content = out_csv.read_text(encoding="utf-8").splitlines()
    # Expect single selection header when top_n=1
    assert any(line.startswith("#Selection name:") for line in content)
    assert any(line == "Selection,X,Y" for line in content)
    # Ensure selection names appear in rows
    assert any(line.startswith("Selection 1,") for line in content)
