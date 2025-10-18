from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

geopandas = pytest.importorskip("geopandas")
gpd = geopandas
from shapely.geometry import Polygon  # noqa: E402

pytest.importorskip("yaml")
from src.scripts.plaque_alignment.app import run_align_and_export  # noqa: E402
from src.scripts.plaque_alignment.config import load_config  # noqa: E402


def _make_minimal_dataset(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    base = tmp_path / "base"
    base.mkdir(parents=True, exist_ok=True)
    # Keypoints: identity-like mapping (src==dst) with small noise
    kp = pd.DataFrame(
        {
            "fixedX": [0.0, 10.0, 50.0, 100.0],
            "fixedY": [0.0, 10.0, 50.0, 100.0],
            "alignmentX": [0.0, 10.0, 50.0, 100.0],
            "alignmentY": [0.0, 10.0, 50.0, 100.0],
        }
    )
    kp_path = base / "keypoints.csv"
    kp.to_csv(kp_path, index=False)

    # Simple square polygon in IF pixels
    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)])
    gdf = gpd.GeoDataFrame({"geometry": [poly]}, crs="EPSG:4326")
    plaques_path = base / "plaques.geojson"
    gdf.to_file(plaques_path, driver="GeoJSON")

    out_sel = base / "out.csv"
    logs_dir = base / "logs"
    return kp_path, plaques_path, out_sel, logs_dir


def _write_cfg(tmp_path: Path, kp: Path, plaques: Path, out_sel: Path, logs_dir: Path) -> Path:
    cfg = {
        "paths": {
            "base_dir": str(tmp_path),
            "keypoints_csv": str(kp.relative_to(tmp_path)),
            "plaque_geojson": str(plaques.relative_to(tmp_path)),
            "out_xenium_format": str(out_sel.relative_to(tmp_path)),
            "logs_dir": str(logs_dir.relative_to(tmp_path)),
        },
        "params": {
            "xenium_pixel_size_um": 0.2125,
            "ransac": {"min_samples": 3, "residual_threshold_px": 1.0, "max_trials": 50},
            "upgrade_to_affine_rmse_px": 5.0,
            "selections": {"decimals": 3, "top_n": None},
        },
        "logging": {"level": "INFO", "timezone": "UTC"},
    }
    cfg_path = tmp_path / "cfg.yaml"
    import yaml

    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return cfg_path


def test_run_align_and_export_integration(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    kp_path, plaques_path, out_sel, logs_dir = _make_minimal_dataset(tmp_path)
    cfg_path = _write_cfg(tmp_path, kp_path, plaques_path, out_sel, logs_dir)

    cfg = load_config(cfg_path)
    rc = run_align_and_export(cfg)
    assert rc == 0
    assert out_sel.exists()
    # Check header and at least one row
    lines = out_sel.read_text(encoding="utf-8").splitlines()
    assert any(line.startswith("#Selection") for line in lines)
    assert any(line == "Selection,X,Y" for line in lines)


def test_cli_end_to_end(tmp_path: Path) -> None:
    kp_path, plaques_path, out_sel, logs_dir = _make_minimal_dataset(tmp_path)
    cfg_path = _write_cfg(tmp_path, kp_path, plaques_path, out_sel, logs_dir)

    env = dict(**os.environ)
    env["PYTHONPATH"] = str((Path.cwd() / "src").resolve())

    cmd = [
        sys.executable,
        "-m",
        "scripts.plaque_alignment.cli",
        "--config",
        str(cfg_path),
        "--log-level",
        "INFO",
    ]
    res = subprocess.run(cmd, cwd=Path.cwd(), capture_output=True, text=True, env=env)
    assert res.returncode == 0, res.stderr
    assert out_sel.exists()
