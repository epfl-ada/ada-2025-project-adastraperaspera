"""
Configuration loading & validation.
Expects a `../../../configs/plaque_alignment.yaml` file for user-editable settings; CLI flags can override.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import os
from pathlib import Path
import typing as t

import yaml

# ===========================================================
# Dataclass models for validating the plaque_alignment.yaml
# ===========================================================


@dataclass
class RansacCfg:
    """RANSAC hyperparameters."""

    min_samples: int
    residual_threshold_px: float
    max_trials: int


@dataclass
class ParamsCfg:
    """Algorithm parameters and tuning knobs."""

    xenium_pixel_size_um: float
    ransac: RansacCfg
    upgrade_to_affine_rmse_px: float
    selections: dict[str, t.Any]


@dataclass
class PathsCfg:
    """
    All file system locations for inputs/outputs.
    All members are absolute Paths after loading & resolution.
    """

    base_dir: Path
    keypoints_csv: Path
    plaque_geojson: Path
    out_xenium_format: Path
    logs_dir: Path


@dataclass
class AppCfg:
    """Top-level application configuration."""

    paths: PathsCfg
    params: ParamsCfg
    logging: dict[str, t.Any]


# =========================
# Helpers
# =========================


def _expand(p: str | Path) -> Path:
    """
    Expand ~ and $VARS then return a Path (without resolving).
    """
    s = str(p)
    return Path(os.path.expandvars(os.path.expanduser(s)))


def _resolve_under(base: Path, p: str | Path) -> Path:
    """
    Resolve a path relative to base if it is not absolute.
    """
    pp = _expand(p)
    return (pp if pp.is_absolute() else (base / pp)).resolve()


def _validate_params(params: ParamsCfg) -> None:
    """
    Validate parameter ranges and raise ValueError on invalid inputs.
    """
    if params.xenium_pixel_size_um <= 0:
        raise ValueError("xenium_pixel_size_um must be > 0")
    if params.ransac.min_samples < 2:
        raise ValueError("ransac.min_samples must be >= 2")
    if params.ransac.residual_threshold_px <= 0:
        raise ValueError("ransac.residual_threshold_px must be > 0")
    if params.ransac.max_trials <= 0:
        raise ValueError("ransac.max_trials must be > 0")
    if params.upgrade_to_affine_rmse_px <= 0:
        raise ValueError("upgrade_to_affine_rmse_px must be > 0")


def _coerce_paths(raw_paths: dict[str, t.Any]) -> PathsCfg:
    """
    Coerce raw dict (from YAML) to a fully-resolved PathsCfg.

    All fields become absolute Paths; relative ones are resolved under base_dir.
    """
    # Base dir first
    base_dir_raw = raw_paths.get("base_dir")
    base_dir = _expand(base_dir_raw).resolve()

    def resolve_key(key: str) -> Path:
        if key not in raw_paths:
            raise KeyError(f"Missing required path '{key}' in config.paths")
        return _resolve_under(base_dir, raw_paths[key])

    return PathsCfg(
        base_dir=base_dir,
        keypoints_csv=resolve_key("keypoints_csv"),
        plaque_geojson=resolve_key("plaque_geojson"),
        out_xenium_format=resolve_key("out_xenium_format"),
        logs_dir=resolve_key("logs_dir"),
    )


def _coerce_params(raw_params: dict[str, t.Any]) -> ParamsCfg:
    """
    Coerce raw dict (from YAML) to ParamsCfg and validate it.
    """
    r_raw = raw_params.get("ransac") or {}
    r = RansacCfg(
        min_samples=int(r_raw.get("min_samples")),
        residual_threshold_px=float(r_raw.get("residual_threshold_px")),
        max_trials=int(r_raw.get("max_trials")),
    )
    params = ParamsCfg(
        xenium_pixel_size_um=float(raw_params.get("xenium_pixel_size_um")),
        ransac=r,
        upgrade_to_affine_rmse_px=float(raw_params.get("upgrade_to_affine_rmse_px")),
        selections=t.cast(dict[str, t.Any], raw_params.get("selections") or {}),
    )
    _validate_params(params)
    return params


# =========================
# Public API
# =========================


def load_config(path: str | Path) -> AppCfg:
    """
    Load an AppCfg from a YAML file.

    Parameters
    ----------
    path : str | Path
        Path to a YAML config file.

    Returns
    -------
    AppCfg
        Fully-coerced configuration with absolute paths and validated params.
    """
    cfg_path = _expand(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}

    raw_paths = raw.get("paths")
    raw_params = raw.get("params")
    raw_logging = raw.get("logging")

    paths = _coerce_paths(raw_paths)
    params = _coerce_params(raw_params)

    logging_cfg: dict[str, t.Any] = {
        "level": raw_logging.get("level", "INFO"),
        "timezone": raw_logging.get("timezone", "Europe/Zurich"),
        **{k: v for k, v in raw_logging.items() if k not in {"level", "timezone"}},
    }

    return AppCfg(paths=paths, params=params, logging=logging_cfg)


def apply_cli_overrides(cfg: AppCfg, args: t.Any) -> AppCfg:
    """
    Apply argparse-style overrides on top of the loaded config.

    Only attributes present and non-None in `args` are applied. The function
    returns a (shallow) updated AppCfg for convenience and also mutates `cfg`.

    Parameters
    ----------
    cfg : AppCfg
        The loaded configuration object to modify.
    args : Any
        An argparse.Namespace (or any object with similarly named attributes).

    Returns
    -------
    AppCfg
        The same instance with fields updated from CLI flags.
    """

    # Helper for attribute presence AND non-None value.
    def has(args: t.Any, name: str) -> bool:
        return hasattr(args, name) and getattr(args, name) is not None

    # --- Paths ---
    b = cfg.paths.base_dir  # absolute
    if has(args, "keypoints"):
        cfg.paths = replace(cfg.paths, keypoints_csv=_resolve_under(b, args.keypoints))
    if has(args, "plaques"):
        cfg.paths = replace(cfg.paths, plaque_geojson=_resolve_under(b, args.plaques))
    if has(args, "out_selections"):
        cfg.paths = replace(cfg.paths, out_xenium_format=_resolve_under(b, args.out_selections))

    # --- Params: pixel size & RANSAC ---
    if has(args, "pixel_size_um"):
        cfg.params = replace(cfg.params, xenium_pixel_size_um=float(args.pixel_size_um))

    r = cfg.params.ransac
    changed = False
    if has(args, "ransac_min_samples"):
        r = replace(r, min_samples=int(args.ransac_min_samples))
        changed = True
    if has(args, "ransac_residual_threshold"):
        r = replace(r, residual_threshold_px=float(args.ransac_residual_threshold))
        changed = True
    if has(args, "ransac_max_trials"):
        r = replace(r, max_trials=int(args.ransac_max_trials))
        changed = True
    if changed:
        cfg.params = replace(cfg.params, ransac=r)
    if has(args, "upgrade_affine_rmse"):
        cfg.params = replace(cfg.params, upgrade_to_affine_rmse_px=float(args.upgrade_affine_rmse))

    # Re-validate after overrides
    _validate_params(cfg.params)

    # --- Logging ---
    if has(args, "log_level"):
        cfg.logging["level"] = str(args.log_level)
    if has(args, "log_tz"):
        cfg.logging["timezone"] = str(args.log_tz)

    return cfg
