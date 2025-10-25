"""
DESeq2-style median-of-ratios normalization for single-cell gene expression data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from pydeseq2.preprocessing import deseq2_norm

from src.utils.logging_utils import logger


def geometric_mean_per_gene(adata: sc.AnnData) -> np.ndarray:
    """
    Compute geometric mean per gene using only positive counts (ignoring zeros).
    
    For each gene, compute the geometric mean across cells using only positive counts.
    If a gene is zero in every cell, its geometric mean is undefined (0).
    
    Args:
        adata (sc.AnnData): Expression data with genes as columns, cells as rows.
        
    Returns:
        np.ndarray: Geometric mean for each gene (genes with all zeros = 0).
    """
    logger.info("Computing geometric mean per gene (ignoring zeros)...")
    
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    
    # Convert to dense if needed and ensure non-negative
    X = np.asarray(X, dtype=np.float64)
    
    # For each gene, compute geometric mean of positive values only
    geometric_means = np.zeros(X.shape[1])  # shape[1] = n_genes
    
    for gene_idx in range(X.shape[1]):
        gene_counts = X[:, gene_idx]
        positive_counts = gene_counts[gene_counts > 0]
        
        if len(positive_counts) > 0:
            # Geometric mean = exp(mean(log(positive_values)))
            geometric_means[gene_idx] = np.exp(np.mean(np.log(positive_counts)))
        else:
            # All zeros for this gene
            geometric_means[gene_idx] = 0.0
    
    n_defined = np.sum(geometric_means > 0)
    logger.info(f"Computed geometric means for {n_defined}/{X.shape[1]} genes")
    
    return geometric_means



def compute_size_factors(adata: sc.AnnData, min_genes: int = 20) -> np.ndarray:
    """
    Compute DESeq2-style size factors using median-of-ratios normalization.
    
    Steps:
    1. Compute geometric mean per gene (ignoring zeros)
    2. Calculate per-cell ratios for genes with defined geometric means
    3. Take median of ratios as cell size factor
    4. Handle cells with too few genes using fallback
    
    Args:
        adata (sc.AnnData): Expression data.
        min_genes (int): Minimum number of genes required for size factor calculation.
        
    Returns:
        np.ndarray: Size factors for each cell.
    """
    logger.info("Computing DESeq2-style size factors...")
    
    # Step 1: Geometric mean per gene
    geometric_means = geometric_mean_per_gene(adata)
    
    # Step 2: Per-cell ratios and size factors
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    
    n_cells, n_genes = X.shape
    size_factors = np.zeros(n_cells)
    
    # Find genes with defined geometric means
    valid_genes = geometric_means > 0
    n_valid_genes = np.sum(valid_genes)
    
    logger.info(f"Using {n_valid_genes} genes with defined geometric means")
    
    if n_valid_genes == 0:
        logger.warning("No genes with defined geometric means found!")
        return np.ones(n_cells)
    
    # Calculate size factors for each cell
    for cell_idx in range(n_cells):
        cell_counts = X[cell_idx, :]
        valid_counts = cell_counts[valid_genes]
        valid_geometric_means = geometric_means[valid_genes]
        
        # Calculate ratios for this cell
        ratios = valid_counts / valid_geometric_means
        
        # Take median of ratios as size factor
        if len(ratios) >= min_genes:
            size_factors[cell_idx] = np.median(ratios)
        else:
            # Too few genes - will be handled later with fallback
            size_factors[cell_idx] = np.nan
    
    # Step 3: Handle cells with too few genes
    nan_mask = np.isnan(size_factors)
    n_nan = np.sum(nan_mask)
    
    if n_nan > 0:
        logger.warning(f"{n_nan} cells with too few genes, using fallback size factors")
        
        # Use median of valid size factors as fallback
        valid_size_factors = size_factors[~nan_mask]
        if len(valid_size_factors) > 0:
            fallback_size_factor = np.median(valid_size_factors)
        else:
            fallback_size_factor = 1.0
            
        size_factors[nan_mask] = fallback_size_factor
    
    # Step 4: Center the factors
    global_median = np.median(size_factors)
    if global_median > 0:
        size_factors = size_factors / global_median
    else:
        logger.warning("Global median is 0, setting all size factors to 1.0")
        size_factors = np.ones_like(size_factors)
    
    logger.info(f"Size factors computed: median={np.median(size_factors):.3f}, "
                f"range=[{np.min(size_factors):.3f}, {np.max(size_factors):.3f}]")
    
    return size_factors


def deseq2_normalize(adata: sc.AnnData, min_genes: int = 20) -> sc.AnnData:
    """
    Apply DESeq2-style median-of-ratios normalization to gene expression data.
    
    Steps:
    1. Compute geometric mean per gene (ignoring zeros)
    2. Calculate per-cell size factors using median-of-ratios
    3. Center size factors by global median
    4. Apply normalization (divide counts by size factors)
    5. Apply log1p transformation for variance stabilization
    
    Args:
        adata (sc.AnnData): Raw expression data.
        min_genes (int): Minimum number of genes required for size factor calculation.
        
    Returns:
        sc.AnnData: Normalized AnnData object with log1p-transformed counts.
    """
    logger.info("🔄 Applying DESeq2-style median-of-ratios normalization...")
    
    # Create a copy to avoid modifying original data
    adata_norm = adata.copy()
    
    # Compute size factors
    size_factors = compute_size_factors(adata_norm, min_genes=min_genes)
    
    # Store size factors in obs
    adata_norm.obs['size_factors'] = size_factors
    
    # Apply normalization: divide each cell's counts by its size factor
    logger.info("Applying size factor normalization...")
    
    X = adata_norm.X
    if sparse.issparse(X):
        # For sparse matrices, normalize row by row
        X_norm = X.copy()
        for i in range(X.shape[0]):
            if size_factors[i] > 0:
                X_norm[i, :] = X[i, :] / size_factors[i]
    else:
        # For dense matrices, vectorized operation
        X_norm = X / size_factors[:, np.newaxis]
    
    adata_norm.X = X_norm
    
    # Apply log1p transformation for variance stabilization
    logger.info("Applying log1p transformation...")
    sc.pp.log1p(adata_norm)
    
    logger.info("✅ DESeq2 normalization complete")
    return adata_norm


def plot_size_factor_distribution(adata: sc.AnnData, figsize: tuple[int, int] = (8, 6)) -> None:
    """
    Plot distribution of size factors after DESeq2 normalization.
    
    Args:
        adata (sc.AnnData): Normalized AnnData object with size factors in obs.
        figsize (tuple[int, int]): Figure size.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if 'size_factors' not in adata.obs.columns:
        raise KeyError("Size factors not found in adata.obs. Run deseq2_normalize first.")
    
    logger.info("Plotting size factor distribution...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Histogram
    sns.histplot(adata.obs['size_factors'], bins=50, ax=ax1, color='skyblue', alpha=0.7)
    ax1.axvline(1.0, color='red', linestyle='--', label='Median (1.0)')
    ax1.set_xlabel('Size Factor')
    ax1.set_ylabel('Number of Cells')
    ax1.set_title('Distribution of Size Factors')
    ax1.legend()
    
    # Box plot
    sns.boxplot(y=adata.obs['size_factors'], ax=ax2, color='lightgreen')
    ax2.set_ylabel('Size Factor')
    ax2.set_title('Size Factor Distribution')
    
    plt.tight_layout()
    plt.show()
    
    # Print summary statistics
    sf = adata.obs['size_factors']
    print(f"Size Factor Summary:")
    print(f"  Mean: {sf.mean():.3f}")
    print(f"  Median: {sf.median():.3f}")
    print(f"  Range: [{sf.min():.3f}, {sf.max():.3f}]")
    print(f"  Std: {sf.std():.3f}")
def pydeseq2_normalize_global(
    adata: sc.AnnData, 
    min_genes: int = 20,
    apply_log1p: bool = True
) -> tuple[sc.AnnData, np.ndarray]:
    """
    Apply DESeq2-style median-of-ratios normalization using PyDESeq2.
    
    This method:
    1. Uses PyDESeq2 to compute median-of-ratios size factors
    2. Applies size factors globally across all slides/time points
    3. Applies log1p transformation for exploration
    4. Preserves raw counts + size factors for DE analysis
    
    Args:
        adata (sc.AnnData): Raw expression data with cells as samples.
        min_genes (int): Minimum number of genes required for size factor calculation.
        apply_log1p (bool): Whether to apply log1p transformation for exploration.
        
    Returns:
        tuple[sc.AnnData, np.ndarray]: 
            - Normalized AnnData object (log1p transformed if apply_log1p=True)
            - Size factors array for use in DE analysis
    """
    logger.info("🔄 Applying PyDESeq2 median-of-ratios normalization...")
    
    # Create a copy to avoid modifying original data
    adata_norm = adata.copy()
    
    # Convert to pandas DataFrame for PyDESeq2 (cells as rows, genes as columns)
    logger.info("Converting to pandas DataFrame for PyDESeq2...")
    X = adata_norm.X
    if sparse.issparse(X):
        X = X.toarray()
    
    # Create DataFrame with cells as rows, genes as columns
    counts_df = pd.DataFrame(
        X, 
        index=adata_norm.obs_names, 
        columns=adata_norm.var_names
    )
    
    logger.info(f"Input data: {counts_df.shape[0]} cells × {counts_df.shape[1]} genes")
    
    # Apply PyDESeq2 median-of-ratios normalization
    logger.info("Computing size factors using PyDESeq2...")
    try:
        normalized_counts, size_factors = deseq2_norm(counts_df)
        
        # Check if size factors are valid (not all NaN)
        if np.isnan(size_factors).all():
            logger.warning("PyDESeq2 returned all NaN size factors, falling back to custom implementation")
            # Fall back to custom implementation
            size_factors = compute_size_factors(adata_norm, min_genes=min_genes)
            # Apply size factors manually
            X_norm = X / size_factors[:, np.newaxis]
            normalized_counts = pd.DataFrame(
                X_norm, 
                index=adata_norm.obs_names, 
                columns=adata_norm.var_names
            )
        else:
            logger.info("✅ PyDESeq2 normalization completed successfully")
            
    except Exception as e:
        logger.warning(f"PyDESeq2 normalization failed: {e}, falling back to custom implementation")
        # Fall back to custom implementation
        size_factors = compute_size_factors(adata_norm, min_genes=min_genes)
        # Apply size factors manually
        X_norm = X / size_factors[:, np.newaxis]
        normalized_counts = pd.DataFrame(
            X_norm, 
            index=adata_norm.obs_names, 
            columns=adata_norm.var_names
        )
    
    # Store size factors in obs
    adata_norm.obs['size_factors'] = size_factors
    adata_norm.obs['size_factors_global'] = size_factors  # Global application
    
    # Apply normalized counts to AnnData
    adata_norm.X = normalized_counts
    
    # Apply log1p transformation for exploration if requested
    if apply_log1p:
        logger.info("Applying log1p transformation for exploration...")
        sc.pp.log1p(adata_norm)
        logger.info("✅ Log1p transformation applied")
    
    # Log size factor statistics
    sf = adata_norm.obs['size_factors']
    logger.info(f"Size factors computed: median={sf.median():.3f}, "
                f"range=[{sf.min():.3f}, {sf.max():.3f}]")
    
    logger.info("✅ PyDESeq2 normalization complete")
    return adata_norm, size_factors


def apply_global_size_factors(
    adata: sc.AnnData, 
    global_size_factors: np.ndarray,
    apply_log1p: bool = True
) -> sc.AnnData:
    """
    Apply globally computed size factors to new data.
    
    This is useful when you want to apply size factors computed from one dataset
    to another dataset (e.g., different slides or time points).
    
    Args:
        adata (sc.AnnData): Raw expression data.
        global_size_factors (np.ndarray): Size factors computed from reference dataset.
        apply_log1p (bool): Whether to apply log1p transformation.
        
    Returns:
        sc.AnnData: Normalized AnnData object.
    """
    logger.info("🔄 Applying global size factors...")
    
    # Create a copy
    adata_norm = adata.copy()
    
    # Store global size factors
    adata_norm.obs['size_factors_global'] = global_size_factors
    
    # Apply normalization
    X = adata_norm.X
    if sparse.issparse(X):
        X = X.toarray()
    
    # Apply size factors
    X_norm = X / global_size_factors[:, np.newaxis]
    adata_norm.X = X_norm
    
    # Apply log1p if requested
    if apply_log1p:
        logger.info("Applying log1p transformation...")
        sc.pp.log1p(adata_norm)
    
    logger.info("✅ Global size factors applied")
    return adata_norm


def prepare_for_de_analysis(
    adata: sc.AnnData, 
    size_factors: np.ndarray
) -> sc.AnnData:
    """
    Prepare data for differential expression analysis.
    
    This preserves raw counts and provides size factors as offsets for DE analysis.
    The size factors can be used as offsets in statistical models.
    
    Args:
        adata (sc.AnnData): Raw expression data.
        size_factors (np.ndarray): Size factors for offset in DE analysis.
        
    Returns:
        sc.AnnData: AnnData object with raw counts and size factors for DE.
    """
    logger.info("🔄 Preparing data for DE analysis...")
    
    # Create a copy with raw counts
    adata_de = adata.copy()
    
    # Store size factors as offsets
    adata_de.obs['size_factors_offset'] = size_factors
    adata_de.obs['log_size_factors'] = np.log(size_factors)
    
    # Ensure we have raw counts (no normalization applied)
    logger.info("Preserving raw counts for DE analysis")
    logger.info(f"Data shape: {adata_de.shape}")
    logger.info(f"Size factors range: [{size_factors.min():.3f}, {size_factors.max():.3f}]")
    
    logger.info("✅ Data prepared for DE analysis")
    return adata_de
