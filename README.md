
## ADA 2025 Project Template

Clean, consistent code is a major grading factor in ADA. This repository is configured with pre-commit hooks and Ruff so that style, quality, and formatting are enforced automatically on every commit.

### Quickstart

```bash
# clone project
git clone git@github.com:epfl-ada/ada-2025-project-adastraperaspera.git
cd ada-2025-project-adastraperaspera

# Create a venv and activate it
python3.11 -m venv ada-venv
source ada-venv/bin/activate
export PATH="/my/path/to/ada-venv/bin:$PATH" # <- modify this to your venv path

# install dependencies (includes Ruff and pre-commit)
pip install -r pip_requirements.txt

# install git hooks
pre-commit install

# run hooks on all files once (recommended)
pre-commit run --all-files
```

### What are pre-commit hooks?

Pre-commit hooks are scripts that run automatically before each commit. They prevent common mistakes (style violations, stray whitespace, large files, and more) from entering the repository. This keeps the codebase clean and consistent for collaborators and graders.

### Hooks in this repo

- Ruff: lint, auto-fix, and format Python according to `ruff.toml`.
- Base safety checks: large-file checks, merge-conflict detection, YAML validation, end-of-file newline, and trailing whitespace.
- Codespell: catches common spelling mistakes in text and code.
- nbstripout: strips notebook outputs so `.ipynb` files remain lightweight and diff-friendly.
- GeoJSON note: commits will print a reminder to consider Git LFS for large `data/*.geojson` files.

You can run Ruff directly as well:

```bash
ruff check --fix .
ruff format .
```

### How to work with the hooks

- Normal development: just commit as usual; fixes will run automatically.
- Run on demand: `pre-commit run --all-files` to validate the whole repo.
- Rarely bypass hooks: `git commit -m "msg" --no-verify` (not recommended; marks will reflect code quality).

### GeoJSON data context

This repo contains `data/Xenium_V1_FFPE_TgCRND8_17_9_months_plaque_polygons.geojson`, added from QuPath, based on slice:1 channel only and a Random Trees classifier (see your recent commit). This can be a large artifact. The hooks will remind you to consider Git LFS for such files. If the file must live in the repo for grading, keep it in `data/` and document its provenance and size here.

### Alignment keypoints CSV (OME-TIFF registration)

`data/Xenium_V1_FFPE_TgCRND8_17_9_months_if_image_keypoints_point.csv` provides matched 2D control points relating two images:

- `Xenium_V1_FFPE_TgCRND8_17_9_months_if_image.ome.tif`
- `morphology_focus.ome.tif`

The CSV columns are:

- `fixedX`, `fixedY`: coordinates in the fixed/reference image
- `alignmentX`, `alignmentY`: coordinates in the other image to be aligned

Convention used here: the morphology-focus image is treated as the fixed/reference image, and the Xenium IF image is the one being aligned. You can use these keypoints to estimate a transform (e.g., similarity/affine) with your preferred library.

Example (Python):

```python
import pandas as pd

keypoints = pd.read_csv(
    "data/Xenium_V1_FFPE_TgCRND8_17_9_months_if_image_keypoints_point.csv"
)
fixed_points = keypoints[["fixedX", "fixedY"]].to_numpy()
moving_points = keypoints[["alignmentX", "alignmentY"]].to_numpy()
# Fit a transform with skimage, OpenCV, or similar
```

### Project structure (orientation)

```
├── data/                       # Project data files
├── src/                        # Source code
│   ├── data/
│   ├── models/
│   ├── utils/
│   └── scripts/
├── tests/                      # Tests
├── results.ipynb               # Results notebook
├── ruff.toml                   # Ruff configuration (lint + format)
├── .pre-commit-config.yaml     # Git hooks configuration
├── pip_requirements.txt        # Python dependencies (incl. dev tools)
└── README.md
```

### Notes on Ruff configuration

- Targets Python 3.11, line length 100, import sorting enabled.
- Enforces naming, bugbear, pyupgrade, comprehensions, pytest style, and annotations.
- Docstring style: Google; missing-docstring rules are relaxed for pragmatism.
- Tests and small scripts are less strict (see `ruff.toml`).

Keeping code clean and consistent is part of your grade—use the hooks locally and push code that passes them.
