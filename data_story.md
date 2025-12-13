# Introduction
This project investigates how Aβ plaques impact the surrounding microenvironment. We aim to develop a quantitative model to precisely describe the influence of the plaques on the surrounding tissue. Knowing which genes and cells are impacted at different plaque distances can refine our understanding of Alzheimer's development. Further, the developers of new drugs can use our model to select realistic targets within the plaque regions accessible from the vasculature.

Our story explores how plaque proximity impacts the cellular, molecular, and tissue environments. We aim to analyze the cell-to-plaque distances and apply rigorous statistical tests to describe the spatial trends in gene expression and cell composition. Further, we look at the interplay between gene expression and cell composition, asking ourselves which of the two phenomena is the underlying cause. Finally, we flip the perspective and benchmark predictive models to infer plaque distance from multigene expression.

# Dataset Analysis

In total, Xenium dataset contains:
- 6 mice
- 6 morphology images
- 1 Immunofluorescence image
- 347 genes
- 351,714 cells
- 78,885,074 transcripts
- 34.7 GB of data
- 0 missing values
- >= 1 zero-transcript gene in each cell (sparse)

## Microscopy data

In this project, we analyze the Xenium dataset from 10X Genomics which contains trascriptomic data accompanied by morphology images. The data comes from sagittal brain slices of 6 mice with 4′,6-diamidino-2-phenylindole (DAPI) fluorescent DNA-binding nucleus dye. Three mice constitute healthy controls (wild type, no induced mutations) at 2.5, 5.7, and 13.4 months of age. The remaining mice are mutated (i.e., transgenic) at 2.5, 5.7, and 17.9 months of age.

The induced mutation forces the murine cells to express the amyloid precursor protein (App) carrying known Alzheimer's disease familial mutations. As a result of mutations, the transgenic mice express up to 5 times more of the endogenous App. This leads to early and aggressive cerebral amyloid beta (Aβ) plaque deposition as soon as 3 months of age. 

The Aβ plaques are then revealed with immunofluorescence (IF) staining, but only in transgenic mice at 17.9 months of age; in the figure below, the plaques appear in red. i

// ToDo: 
// 1. Add interactive zoom-in
// 2. Make background transparent so that it blends into the dark theme of the website
// Current plot was obtained with: src.scripts.visualization.plots.make_image_grid; exact inputs can be found in results.ipynb
<p align="center">
  <img src="figures/microscopy_6_mice.png" width="480">
  <br><em>Microscopy images of 6 mice</em>
</p>

Since each brain slice comes from a different mouse, the inter-mouse variation in brain morphology is very significant. Our best attempt to align a pair of most similar mice (transgenic at 17.9 and 5.7 months) resulted in a Root Mean Square Error (RMSE) of 3,390 µm, which is over 50 times larger than the median cell to plaque distance (61 µm).

// ToDo: 
// 1. Add interactive zoom-in
// 2. Make background transparent so that it blends into the dark theme of the website
// Current plot is a simple cropped screenshot in a png format
<p align="center">
  <img src="src/data/figures/Tg_17_Tg_5_alignment.png" width="160">
  <br><em>Attempted alignment of Tg 5.7 months old mouse onto Tg 17.9 months old mouse</em>
</p>

// ToDo: Cell-to-plaque distance map
// Add plot based on src.scripts.visualization.plots.plot_cell_to_plaque_map_visible; exact inputs can be found in results.ipynb
// Make it interactive by addding zoom
// Make background transparent so that it blends into the dark theme of the website

// ToDo: Add distribution plot of cell-to-plaque distance
// Based on src.scripts.visualization.plots.analyze_plaque_distance; exact inputs can be found in results.ipynb
// Use only linear scale
// Zoom in, transparent background

## Gene expression

We are dealing with a spatial transcriptomics dataset which contains single cell gene expression measurements of 347 genes, out of which 248 represent markers for 8 main cell types, canonical neuronal cortical layer markers, and non-neuronal markers; 83 genes related to activated microglia and astrocytes; and 16 PIGs curated from primary literature. The figure below reveals the individual cells as filled circles with clustering component.

// ToDo: 
// 1. Add interactive zoom-in
// 2. Make background transparent so that it blends into the dark theme of the website
// Current plot was obtained with: src.scripts.visualization.plots.make_image_grid; exact inputs can be found in results.ipynb
<p align="center">
  <img src="figures/microscopy_cells_unified.png" width="480">
  <br><em>Single cell images of 6 mice</em>
</p>

First, let us explore the data by visualizing the distribution of log1p-transformed transcript counts for the 16 PIGs against the average distribution of all 347 genes. We will focus on the mouse with the most advanced stage of the Alzheimer's disease (transgenic at 17.9 months of age).

// ToDo: 
// 1. Add interactive zoom-in
// 2. Make background transparent so that it blends into the dark theme of the website
// Current plot was obtained with: src.scripts.visualization.plots.plot_gene_distributions; exact inputs can be found in results.ipynb
<p align="center">
  <img src="figures/Expression_Distribution.png" width="480">
  <br><em>Probability distribution of log1p-transformed PIG transcript counts</em>
</p>

We can see that in 12 out of 16 PIGs, the proportion of cells with zero transcript count is higher than average, which is expected as they are indicators of Alzheimer's disease. 

This motivates us to focus on the spatial effects of plaque distance on the expression level of PIGs.


# Research Questions

## Predicting Cell Composition

### Cell clustering

## Predicting Gene Expression

// ToDo: 
// Add a plot src.scripts.visualization.plots.plot_gene_trends; all 16 PIGs should be plotted
// Zoom in, transparent background

// Todo:
// Add a plot src.scripts.visualization.plots.plot_expression_heatmap; exact inputs can be found in results.ipynb
// Zoom in, transparent background

## Predicting Plaque distance