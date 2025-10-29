from __future__ import annotations

import argparse

from src.scripts.plaque_alignment.app import run_align_and_export
from src.scripts.plaque_alignment.config import apply_cli_overrides, load_config


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="xenium-plaque",
        description="Align IF plaque polygons to Xenium morphology space and export selections.",
    )
    p.add_argument("--config", default="configs/plaque_alignment.yaml", help="Path to YAML config.")
    p.add_argument("--keypoints")
    p.add_argument("--plaques")
    p.add_argument("--cells")
    p.add_argument("--out-geojson")
    p.add_argument("--out-selections")
    p.add_argument("--pixel-size-um", type=float)
    p.add_argument("--ransac-min-samples", type=int)
    p.add_argument("--ransac-residual-threshold", type=float)
    p.add_argument("--ransac-max-trials", type=int)
    p.add_argument("--upgrade-affine-rmse", type=float)
    p.add_argument("--log-level")
    p.add_argument("--log-tz")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    cfg = apply_cli_overrides(cfg, args)
    return run_align_and_export(cfg)


if __name__ == "__main__":
    raise SystemExit(main())
