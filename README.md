# Team ADAstraPerAspera: Milestone P2

## Introduction

The goal of this project is to build upon the “proximal vs. distal” analysis in Xenium’s application note and inspect the continuous change in gene expression levels and cell type proportions as we move away from the nearest amyloid-beta plaque.

## Developer Best Practices
- [ ] I remember that code quality in ADA *IS GRADED!*
- [ ] I don't use `print()`; instead I always use the logger from `src/utils/logging_utils.py`, for example:

    ```python
    from src.utils.logging_utils import logger
    logger.info("Hello, world!")
    ```

- [ ] I always install (`pre-commit install`) pre-commit hooks and never skip them when committing.
- [ ] I always leave plenty of comments in the code.
- [ ] Every public class, method, file, and function has a docstring.
- [ ] I always use type annotations in function signatures, for example:

    ```python
    def add(a: int, b: int) -> int:
        return a + b
    ```

- [ ] Whenever I commit code, I also add unit and integration tests for any new functionality.
- [ ] I never run `git add .`; instead, I always run `git status` to see what files have changed and then add them one by one.
- [ ] I never commit secrets, API keys, or irrelevant files such as duplicates, caches, boilerplate, etc.


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

# install git hooks
pre-commit install

# run hooks on all files once (recommended)
pre-commit run --all-files
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
│   └── config.yaml                                # Config used by the plaque-alignment CLI.
├── src
│   ├── data
│   │   └── Xenium_V1_FFPE_TgCRND8_17_9_months
│   │       ├── Xenium_V1_FFPE_TgCRND8_17_9_months_if_image.qpdata   # QuPath project: 11 negative rects + 9 positive polygons + classifier outputs (IF space).
│   │       ├── image_keypoints.csv                # 26 matched control points (morphology↔IF).
│   │       ├── plaque_polygons.csv                # 1,938 plaque polygons transformed into morphology coords (post-alignment).
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

## Data Preprocessing

### Plaque classification from Xenium IF image

- In Xenium Alzheimer Disease experiment, the amyloid beta plaques were stained with antibodies and later visualized with immunofluorescence.
- As a result, 10X Genomics company provides `Xenium_V1_FFPE_TgCRND8_17_9_months_if_image.ome.tif` file with the IF data.
- This file contains two channels called `slice:1` and `slice:3`.
- `slice:1` contains coarse, unevenly spaced groups of bright pixels that represent the amyloid beta plaques.
- `slice:3` contains small and much more evenly spaced bright pixels that represent the nuclei of the cells.
- To identify the plaques, we restricted the analysis only to the `slice:1` channel.
- Next, we labeled the data.
- We produced 11 data points carrying rectangles around plaque-free regions (negative samples).
- Similarly, we labeled 9 data points with hand-drawn polygons around plaques (positive samples).
- We then used the object classifier from QuPath to use the provided training data for training.
- We selected the default implementation of the Random trees based classifier due to the limited amount of training data.
- We then used the classifier to identify the plaques in the `slice:1` channel.
- The resulting objects (along with the ground truth labels used to train the classifier) are stored in `src/data/Xenium_V1_FFPE_TgCRND8_17_9_months/Xenium_V1_FFPE_TgCRND8_17_9_months_if_image.qpdata`
- Only the predicted plaque polygons are exported to `src/data/Xenium_V1_FFPE_TgCRND8_17_9_months/qupath_plaque_polygons.geojson`
- However, these coordinates are mismatched from the morphology coordinate space.
- And morphology coordinate space is also used to identify the cell locations.
- Thus, as is, we cannot use the obtained plaque polygons for further analysis.
- In the following section, we describe how we aligned the IF image to the morphology image.

### Aligning the IF image to the morphology image

- To align IF image to the morphology space, we produced 26 keypoint pairs around important anatomical landmarks that constitute a map between the two images.
- In other words, each pair contains closely visually aligned points on the IF image and the morphology image.
- These keypoints are stored in `src/data/Xenium_V1_FFPE_TgCRND8_17_9_months/image_keypoints.csv`
- To reiterate, this file provides matched 2D control points relating two images:
- `Xenium_V1_FFPE_TgCRND8_17_9_months_if_image.ome.tif`
- `morphology_focus.ome.tif`
The CSV columns of this file are:
- `fixedX`, `fixedY`: coordinates in the fixed/reference image, which is the morphology image in our case.
- `alignmentX`, `alignmentY`: coordinates in the other image to be aligned, which is the IF image in our case.
- We would now like to learn a transformation function that would map the keypoints from the IF image as closely as possible on average to their corresponding counterparts in the morphology image.
- Namely, these key points were used to train a SimilarityTransform (RANSAC) model (which proved to perform better than an AffineTransform), resulting in a root mean squared error of 15.213 pixels or 15.213 * 0.2125 µm = 3.231 µm, which is considered acceptable for our purposes.
- To repeat this analysis, run the following simple command below:

```bash
export PYTHONPATH="$PWD/src"
python -m scripts.plaque_alignment.cli --config configs/config.yaml
```

- Now that we have a mapping from the IF coordinate space into the morphology coordinate space, we can transform the plaque polygons into the morphology coordinate space.
- As a result, we obtained 1938 exterior polygons around the plaques.
- The transformed plaque polygons which are stored in `src/data/Xenium_V1_FFPE_TgCRND8_17_9_months/plaque_polygons.csv`
- Importing these transformed polygons into the Xenium Explorer reveals close visual alignment with the morphogy image, as expected.

### Ruff Configuration

- Targets Python 3.11, line length 100, import sorting enabled.
- Enforces naming, bugbear, pyupgrade, comprehensions, pytest style, and annotations.
- Docstring style: Google; missing-docstring rules are relaxed for pragmatism.
- Tests and small scripts are less strict (see `ruff.toml`).
