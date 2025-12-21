import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Polygon

from src.scripts.preprocessing.plaque_distance import (
    compute_cell_to_plaque_distance,
    load_plaque_polygons,
)


@pytest.fixture
def tmp_plaque_csv(tmp_path):
    """Create a temporary CSV mimicking plaque polygon export."""
    csv_path = tmp_path / "plaques.csv"
    with open(csv_path, "w") as f:
        f.write("Selection,X,Y\n")
        f.write("Selection 1,0,0\n")
        f.write("Selection 1,0,1\n")
        f.write("Selection 1,1,1\n")
        f.write("Selection 1,1,0\n")
        f.write("Selection 2,5,5\n")
        f.write("Selection 2,6,5\n")
        f.write("Selection 2,6,6\n")
        f.write("Selection 2,5,6\n")
    return csv_path


def test_load_plaque_polygons(tmp_plaque_csv):
    """Ensure plaque polygons CSV loads correctly into a GeoDataFrame."""
    gdf = load_plaque_polygons(tmp_plaque_csv)
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert "geometry" in gdf.columns
    assert len(gdf) == 2
    assert gdf.geometry.iloc[0].geom_type == "Polygon"


def test_load_plaque_invalid(tmp_path):
    """Ensure informative error if geometry columns are missing."""
    bad_csv = tmp_path / "bad.csv"
    pd.DataFrame({"foo": [1, 2], "bar": [3, 4]}).to_csv(bad_csv, index=False)
    with pytest.raises(ValueError, match="CSV must have 'Selection', 'X', 'Y' columns"):
        load_plaque_polygons(bad_csv)


@pytest.fixture
def dummy_cells():
    """Create a fake cell DataFrame."""
    return pd.DataFrame(
        {
            "cell_id": ["A", "B", "C"],
            "x_centroid": [0.5, 5.5, 10.0],
            "y_centroid": [0.5, 5.5, 10.0],
        }
    )


@pytest.fixture
def dummy_plaques():
    """Create a simple GeoDataFrame of plaques."""
    polys = [
        Polygon([(0, 0), (0, 1), (1, 1), (1, 0)]),
        Polygon([(5, 5), (6, 5), (6, 6), (5, 6)]),
    ]
    return gpd.GeoDataFrame({"geometry": polys}, crs="EPSG:4326")


def test_compute_distance(dummy_cells, dummy_plaques):
    """Check correct distance calculation."""
    df = compute_cell_to_plaque_distance(dummy_cells, dummy_plaques)
    assert "distance_to_plaque" in df.columns
    assert pytest.approx(df.loc[0, "distance_to_plaque"], rel=1e-3) == 0.0
    assert pytest.approx(df.loc[1, "distance_to_plaque"], rel=1e-3) == 0.0
    assert df.loc[2, "distance_to_plaque"] > 3.0


def test_compute_distance_missing_cols(dummy_plaques):
    """Ensure function raises if missing x/y columns."""
    df = pd.DataFrame({"cell_id": ["A"], "x": [0.1]})
    with pytest.raises(KeyError):
        compute_cell_to_plaque_distance(df, dummy_plaques)
