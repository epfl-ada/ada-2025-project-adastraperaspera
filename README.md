# [Quantifying How Amyloid-beta Plaques Reshape the Tissue Microenvironment](https://epfl-ada.github.io/ada-2025-project-adastraperaspera/)

## Abstract

Amyloid beta plaques (Aβ) are a known hallmark of Alzheimer's disease (AD) with known effects including changes in gene expression, glial activation and neuronal death. However, the bulk of the existing research only examines this influence across rough distance bins.  Building a finer model of plaque-induced microenvironment has important downstream applications. With the knowledge of *which* cells and *which* genes respond *where* around the plaques, drug developers can pre-filter therapeutic targets accessible from the vasculature.

In this project, we utilize Xenium murine transcriptomics dataset to develop and interpret quantitative models. We investigate (1) temporal and spatial trends in gene expression, and (2) time and distance-dependent changes in cell composition. We describe the interplay between the two and provide a plausible causal relationship. Finally, we flip the direction and benchmark predictive models to infer plaque distance from multiple modalities. 

## Full Analysis

The full report with interactive visuals is available on [our website](https://epfl-ada.github.io/ada-2025-project-adastraperaspera/). Derivations and supplementary analysis can be found in [`results.ipynb`](https://github.com/epfl-ada/ada-2025-project-adastraperaspera/blob/main/results.ipynb).

## Research Questions

In this project, we explore the following research questions (RQs).

#### **How does cell-type composition change with plaque proximity?**

- What clusters are under- and over-represented near plaques?
- How quickly does the cluster proportion change with distance?

#### **How are cell-type composition, PIG expression, and plaque distance related?**

- Are the changes in Plaque-Induced Genes (PIGs) mostly due to the cell-type composition changes?
- Does distance to plaque still have an effect after fixing the cell type?

#### **How does the Plaque Induced Gene (PIG) expression change in plaque proximity?**

- Which PIGs have the highest differential expression around the plaques?
- How quickly does the PIG expression decrease with distance?

#### **When modeling plaque distance, which feature modalities are most important?**

- When comparing morphology, spatial coordinates, and cell type, which has the strongest predictive power?
- Do the predictive models generalize in a spatial cross-validation setup?

#### **How does the gene expression change with age for each cell type and mouse group?**

- Which cell types show the largest differences in gene expression between transgenic (Tg) and wild type (Wt) mice?
- What are the age-specific changes in gene expression in Tg mice?

## Datasets

- [**10X Genomics Xenium AD mouse brain dataset**](https://www.10xgenomics.com/datasets/xenium-in-situ-analysis-of-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard)
  - 6 mice (Wt at 2.5, 5.7, 13.4 months; Tg at 2.5, 5.7, 17.9 months)
  - 347 genes, 351,714 cells, 78,885,074 transcripts
  - 6 morphology images and 1 immunofluorescence (IF) image
  - 1,736 Aβ plaques after quality control
- **Plaque annotations and alignment**
  - A [QuPath](https://qupath.github.io/)-based annotation file with labeled plaque and plaque-free regions
  - Landmark pairs of points for morphology/IF image alignment
  - Plaque polygons after alignment, convexity-enforcement, and quality control

## Methods

### Preprocessing and spatial features

- **Plaque detection**: random forest-based segmentation trained on 20 hand-labeled polygons (11 plaque-free regions and 9 plaques).
- **Coordinate alignment**: RANSAC-based transformation from IF coordinate space to morphology coordinate space using 26 hand-labeledlandmark pairs. Resulted in the alignment RMSE of 3.2 µm.
- **Cell-to-plaque distance**: Euclidean distance from Xenium-provided cell centroid to the nearest boundary of plaque polygon.
- **Extra features**
  - Plaque geometry: area, perimeter, major axis length, orientation.
  - Multi-plaque proximity within $R = 61$ µm: plaque count and average distance.
  - Neighborhood expression: mean expression of 15 other PIGs in $100$ nearest neighbors.

### Statistics and modeling

- **Joint clustering**: Dimensionality reduction with PCA -> k-Nearest Neighbors graph (k = 15) -> Leiden clustering -> UMAP for visualization.
- **Cell type composition vs distance**: logistic regression per Leiden cluster with multiple testing correction.
- **PIG spatial pattern analysis**: distance binning, Analysis of Variance (ANOVA) across bins, and per-gene regressions with the Benjamini-Hochberg FDR correction.
- **Distance prediction**: LASSO, ElasticNet, Random Forest, XGBoost, Partial Least Squares (PLS), and a linear regression on gene expression. Ablation analysis on additional features includingmorphology, spatial coordinates, and cluster labels.
- **Spatial cross-validation**: rectangular tile-based segmentation for generalizability analysis. Regression model ablations (within-tile cell permutation, replacement of nearest neighbors with farthest neighbors) to probe the impact of local environment.

## Timeline and contributions

```
Milestone P3 - Analyses, Ownership, and Due-Point 2025 Timeline 
  [P2]                                           [P3]
    │                                              │
  ├─────────┬─────────┬─────────┬─────────┬─────────┤►
  Nov10   Nov18     Nov26     Dec4      Dec12     Dec20
  start    ▲1        ▲2        ▲3        ▲4        ▲5
           
  Nov18: Preprocessing + clustering                        [Alexander]      Progress: [██░░░░░░░░]
  Nov26: Multi-mouse analysis + individualized modeling    [Sogand]         Progress: [████░░░░░░]
  Dec4: Feature engineering + regression modeling          [Zayed]          Progress: [██████░░░░]
  Dec12: Spatial cross-validation setup + ablations        [Walid]          Progress: [████████░░]
  Dec20: Frontend development + CI/CD                      [Rosa]           Progress: [██████████]

  Key:
    ▲1..▲5  deliverable due-points
    Progress blocks: ░ not started | █ done
```

- **Alexander**: 
  - Transcriptomics data acquisition and preprocessing; 
  - Plaque classification, keypoint labelling and alignment, quality control; 
  - Joint cell clustering and cell type annotation;
  - Results interpretation, repository linting, and reproducibility setup.
- **Sogand**: 
  - Independent cell clustering;
  - Distance to plaque inference with regression and tree-based models;
  - Generalization capability analysis of regression and tree-based models; 
  - Model extension to multiple mice (WT and Tg, 2.5 to 17.9 months) with subject-specific pattern detection.
- **Rosa**: 
  - Writing and editing the report;
  - Interactive plot development;
  - Frontend development;
  - GitHub Pages CI/CD.
- **Zayed**:
  - Distributional modeling, normalization,and anomaly detection;
  - Feature engineering (plaque geometry, multi-plaque proximity, neighborhood gene expression);
  - Nested regression modeling;
  - False Discovery Rate (FDR) corrections.
- **Walid**: 
  - Decision tree interpretation;
  - Regression model diagnostics and ablations;
  - Spatial cross-validation setup;
  - Code testing.

## Environment Setup

Follow the commands below for a quick start:

```shell
# Clone the repository
git clone git@github.com:epfl-ada/ada-2025-project-adastraperaspera.git
cd ada-2025-project-adastraperaspera

# Create and activate a Python 3.11 virtual environment
python3.11 -m venv ada-venv
source ada-venv/bin/activate

# Ensure the venv bin is on PATH
export PATH="$PWD/ada-venv/bin:$PATH"

# Install Python dependencies
pip install -r pip_requirements.txt

# Enable Jupyter Table of Contents (TOC)
jupyter contrib nbextension install --sys-prefix
jupyter nbextensions_configurator enable --sys-prefix
jupyter nbextension enable toc2/main --sys-prefix

# Set common environment variables
source activate_env.sh

# Install pre-commit hooks
pre-commit install

# Run hooks once on all files
pre-commit run --all-files
```

### Environment variables

Set the Xenium data path:

```shell
export BASE_XENIUM_DIR=/path/to/xenium/data
```

### Reproducibility

To reproduce our results, follow the steps in [`results.ipynb`](https://github.com/epfl-ada/ada-2025-project-adastraperaspera/blob/main/results.ipynb).

## Frontend setup

```shell
cd frontend
npm ci
npm run build
```

### View the website locally

Set up the development server:

```shell
cd frontend
npm run dev
```

Open the local URL available from terminal (usually `http://localhost:5173`).


## Testing

```shell
export PYTHONPATH="$PWD/src"
pytest -q

# With coverage report
pytest --cov=src --cov-report=term-missing
```

Notes:

- Tests automatically add `src` to `PYTHONPATH` via `tests/conftest.py`.
- GeoJSON I/O uses the `pyogrio` engine by default via GeoPandas.
- Ensure dependencies from `pip_requirements.txt` are installed (includes PyYAML for config parsing).

## Repository structure 

```shell
.
├── .github/
│   └── workflows/
│       └── static.yml                 # Deploy frontend to GitHub Pages
├── frontend/                          # Vite + React frontend
│   ├── src/                           # Frontend source code
│   └── public/
│       └── plots/                     # Pre-rendered interactive HTML plots
├── configs/                           # YAML configs
│   ├── column_annotations.yaml
│   ├── plaque_alignment.yaml
│   ├── plaque_preprocessing.yaml
│   └── statistical_constants.yaml
├── src/                               # Python code
│   ├── data/                          # Data tables, figures, and raw microscopyimages
│   │   ├── download.py
│   │   ├── Alzeimer_xenium_predicted_cell_types_with_ids.csv
│   │   ├── Cell_label_reference.xlsx
│   │   ├── pig_by_distance_interactive.html
│   │   ├── raw_images/
│   │   └── figures/
│   ├── scripts/
│   │   ├── analysis/
│   │   ├── preprocessing/
│   │   ├── plaque_alignment/
│   │   ├── plaque_signature/
│   │   ├── research_questions/
│   │   └── visualization/
│   └── utils/
│       └── logging_utils.py
├── tests/                             # Unit and integration tests
├── docs/                              # Setup notes and developer tutorial
├── results.ipynb                      # Main notebook
├── activate_env.sh                    # Virtual environment activation
├── pip_requirements.txt               # Python dependencies
├── ruff.toml                          # Ruff linter configuration
├── .pre-commit-config.yaml            # Pre-commit hooks (code styling, and linting)
├── .gitignore
└── README.md
```
