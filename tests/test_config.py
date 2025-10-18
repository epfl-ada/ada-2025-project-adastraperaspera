from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("yaml")
import yaml  # noqa: E402

from src.scripts.plaque_alignment.config import (  # noqa: E402
    AppCfg,
    apply_cli_overrides,
    load_config,
)


def _write_yaml(p: Path, data: dict) -> None:
    p.write_text(yaml.safe_dump(data), encoding="utf-8")


def _base_cfg_dict(tmp_path: Path) -> dict:
    base_dir = tmp_path / "base"
    base_dir.mkdir(parents=True, exist_ok=True)
    return {
        "paths": {
            "base_dir": str(base_dir),
            "keypoints_csv": "keypoints.csv",
            "plaque_geojson": "plaques.geojson",
            "out_xenium_format": "out.csv",
            "logs_dir": "logs",
        },
        "params": {
            "xenium_pixel_size_um": 0.5,
            "ransac": {
                "min_samples": 3,
                "residual_threshold_px": 0.5,
                "max_trials": 50,
            },
            "upgrade_to_affine_rmse_px": 2.0,
            "selections": {"decimals": 3, "top_n": None},
        },
        "logging": {"level": "INFO", "timezone": "UTC"},
    }


def test_load_config_resolves_paths_and_validates(tmp_path: Path) -> None:
    cfg_dict = _base_cfg_dict(tmp_path)
    cfg_path = tmp_path / "cfg.yaml"
    _write_yaml(cfg_path, cfg_dict)

    cfg = load_config(cfg_path)
    assert isinstance(cfg, AppCfg)
    assert cfg.paths.base_dir.is_absolute()
    assert cfg.paths.keypoints_csv.is_absolute()
    assert cfg.paths.plaque_geojson.is_absolute()
    assert cfg.paths.out_xenium_format.is_absolute()
    assert cfg.paths.logs_dir.is_absolute()
    assert cfg.params.xenium_pixel_size_um == 0.5
    assert cfg.params.ransac.min_samples == 3
    assert cfg.logging["level"] == "INFO"


def test_apply_cli_overrides_updates_fields(tmp_path: Path) -> None:
    cfg_dict = _base_cfg_dict(tmp_path)
    cfg_path = tmp_path / "cfg.yaml"
    _write_yaml(cfg_path, cfg_dict)
    cfg = load_config(cfg_path)

    args = SimpleNamespace(
        keypoints="new_kp.csv",
        plaques="new_plaques.geojson",
        out_selections="results/out.csv",
        pixel_size_um=0.25,
        ransac_min_samples=5,
        ransac_residual_threshold=0.25,
        ransac_max_trials=100,
        upgrade_affine_rmse=1.0,
        log_level="DEBUG",
        log_tz="Europe/Zurich",
    )
    cfg2 = apply_cli_overrides(cfg, args)
    # same object mutated and returned
    assert cfg2 is cfg
    assert cfg.paths.keypoints_csv.name == "new_kp.csv"
    assert cfg.paths.plaque_geojson.name == "new_plaques.geojson"
    assert cfg.paths.out_xenium_format.name == "out.csv"
    assert cfg.params.xenium_pixel_size_um == 0.25
    assert cfg.params.ransac.min_samples == 5
    assert cfg.params.ransac.max_trials == 100
    assert cfg.params.upgrade_to_affine_rmse_px == 1.0
    assert cfg.logging["level"] == "DEBUG"
    assert cfg.logging["timezone"] == "Europe/Zurich"


def test_load_config_invalid_params_raise(tmp_path: Path) -> None:
    bad = _base_cfg_dict(tmp_path)
    bad["params"]["xenium_pixel_size_um"] = 0.0
    cfg_path = tmp_path / "cfg_bad.yaml"
    _write_yaml(cfg_path, bad)
    with pytest.raises(ValueError, match="Invalid configuration format"):
        load_config(cfg_path)
