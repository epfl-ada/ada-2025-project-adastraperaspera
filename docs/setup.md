## Quickstart

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

# Enable Jupyter Notebook Table of Contents (TOC)
jupyter contrib nbextension install --sys-prefix
jupyter nbextensions_configurator enable --sys-prefix
jupyter nbextension enable toc2/main --sys-prefix

# install git hooks
pre-commit install

# run hooks on all files once (recommended)
pre-commit run --all-files
```

# Environment Variables

Example of an environment variable file:
```
BASE_XENIUM_DIR=your/path/to/xenium/data
```


## Testing

We provide unit and integration tests for our codebase.
To run tests locally, use the following command:

```bash
export PYTHONPATH="$PWD/src"
# activate venv (see Quickstart), then:
pytest -q

# with coverage report
pytest --cov=src --cov-report=term-missing
```

Notes:
- Tests automatically add `src` to `PYTHONPATH` via `tests/conftest.py`.
- GeoJSON I/O uses the `pyogrio` engine by default via GeoPandas.
- Ensure dependencies from `pip_requirements.txt` are installed (includes PyYAML for config parsing).

## Project structure (orientation)

```
.
├── configs
│   └── plaque_alignment.yaml                                # Config used by the plaque-alignment CLI.
├── src
│   ├── data
│   │   └── Xenium_V1_FFPE_TgCRND8_17_9_months
│   │       ├── Xenium_V1_FFPE_TgCRND8_17_9_months_if_image.qpdata   # QuPath project: 11 negative rects + 9 positive polygons + classifier outputs (IF space).
│   │       ├── image_keypoints.csv                # 26 matched control points (morphology↔IF).
│   │       ├── plaque_polygons.csv                # 1,938 plaque polygons transformed into morphology coords (post-alignment).
│   │       ├── brain_polygon.csv                  # The coordinates delimiting the brain region.
│   │       └── qupath_plaque_polygons.geojson     # Predicted plaque polygons exported from QuPath in IF-image coords (pre-alignment).
│   ├── scripts
│   │   └── plaque_alignment
│   │       ├── app.py                             # Fits a similarity or an affine transform, computes RMSE, transforms polygons.
│   │       ├── cli.py                             # Entry point for CLI.
│   │       ├── config.py                          # YAML config loader.
│   │       └── utils.py                           # I/O, geometry helpers.
│   └── utils
│       └── logging_utils.py                       # Shared logger setup.
├── tests                                           # Unit and integration tests.
├── results.ipynb                                   # Notebook to showcase the results.
├── ruff.toml                                       # Ruff config.
├── .pre-commit-config.yaml                         # Git hooks to auto-run Ruff/formatting on commits.
├── pip_requirements.txt                            # Python dependencies.
└── README.md                                       # Project overview and detailed notes.

```

### Reproducibility

To rerun the preprocessing and analysis refer to results.ipynb

### Ruff Configuration

- Targets Python 3.11, line length 100, import sorting enabled.
- Enforces naming, bugbear, pyupgrade, comprehensions, pytest style, and annotations.
- Docstring style: Google; missing-docstring rules are relaxed for pragmatism.
- Tests and small scripts are less strict (see `ruff.toml`).
