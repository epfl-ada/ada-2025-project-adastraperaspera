"""
Tests for normalization module functions.
"""

import numpy as np
import pytest
import scanpy as sc
from scipy import sparse

from src.scripts.preprocessing.normalization import (
    apply_global_size_factors,
    compute_size_factors,
    deseq2_normalize,
    geometric_mean_per_gene,
    plot_size_factor_distribution,
    prepare_for_de_analysis,
    pydeseq2_normalize_global,
)


@pytest.fixture
def dummy_adata():
    """Create a dummy AnnData object for testing."""
    # Create a small expression matrix (5 cells, 10 genes)
    np.random.seed(42)
    X = np.random.poisson(5, size=(5, 10)).astype(float)

    # Add some zeros to test geometric mean calculation
    X[0, 0] = 0  # First gene in first cell is zero
    X[1, 1] = 0  # Second gene in second cell is zero

    # Create AnnData object
    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(10)]
    adata.obs_names = [f"Cell_{i}" for i in range(5)]

    return adata


@pytest.fixture
def sparse_adata():
    """Create a dummy sparse AnnData object for testing."""
    np.random.seed(42)
    X = np.random.poisson(3, size=(3, 5)).astype(float)
    X[X < 2] = 0  # Make it sparse

    # Convert to sparse matrix
    X_sparse = sparse.csr_matrix(X)

    adata = sc.AnnData(X=X_sparse)
    adata.var_names = [f"Gene_{i}" for i in range(5)]
    adata.obs_names = [f"Cell_{i}" for i in range(3)]

    return adata


@pytest.fixture
def all_zeros_adata():
    """Create an AnnData object with all zeros to test edge cases."""
    X = np.zeros((3, 5))
    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(5)]
    adata.obs_names = [f"Cell_{i}" for i in range(3)]

    return adata


def test_geometric_mean_per_gene(dummy_adata):
    """Test geometric mean calculation per gene."""
    geometric_means = geometric_mean_per_gene(dummy_adata)

    # Check output shape
    assert len(geometric_means) == dummy_adata.n_vars
    assert geometric_means.shape == (10,)

    # Check that geometric means are non-negative
    assert np.all(geometric_means >= 0)

    # Check that we have some defined geometric means
    assert np.sum(geometric_means > 0) > 0


def test_geometric_mean_per_gene_sparse(sparse_adata):
    """Test geometric mean calculation with sparse data."""
    geometric_means = geometric_mean_per_gene(sparse_adata)

    # Check output shape
    assert len(geometric_means) == sparse_adata.n_vars
    assert geometric_means.shape == (5,)

    # Check that geometric means are non-negative
    assert np.all(geometric_means >= 0)


def test_geometric_mean_per_gene_all_zeros(all_zeros_adata):
    """Test geometric mean calculation with all zeros."""
    geometric_means = geometric_mean_per_gene(all_zeros_adata)

    # All geometric means should be zero
    assert np.all(geometric_means == 0)


def test_compute_size_factors(dummy_adata):
    """Test size factor computation."""
    size_factors = compute_size_factors(dummy_adata, min_genes=5)

    # Check output shape
    assert len(size_factors) == dummy_adata.n_obs
    assert size_factors.shape == (5,)

    # Check that size factors are positive
    assert np.all(size_factors > 0)

    # Check that size factors are reasonable (not too extreme)
    assert np.all(size_factors > 0.1)
    assert np.all(size_factors < 10.0)


def test_compute_size_factors_sparse(sparse_adata):
    """Test size factor computation with sparse data."""
    size_factors = compute_size_factors(sparse_adata, min_genes=2)

    # Check output shape
    assert len(size_factors) == sparse_adata.n_obs
    assert size_factors.shape == (3,)

    # Check that size factors are positive
    assert np.all(size_factors > 0)


def test_compute_size_factors_min_genes(dummy_adata):
    """Test size factor computation with different min_genes parameter."""
    # Test with very high min_genes (should trigger fallback)
    size_factors = compute_size_factors(dummy_adata, min_genes=100)

    # Should still return valid size factors
    assert len(size_factors) == dummy_adata.n_obs
    assert np.all(size_factors > 0)


def test_deseq2_normalize(dummy_adata):
    """Test DESeq2 normalization function."""
    adata_norm = deseq2_normalize(dummy_adata, min_genes=5)

    # Check that original data is not modified
    assert dummy_adata.X is not adata_norm.X

    # Check that size factors are stored
    assert "size_factors" in adata_norm.obs.columns

    # Check that data is log1p transformed
    assert np.all(adata_norm.X >= 0)  # log1p ensures non-negative

    # Check that normalization was applied (values should be different from raw)
    assert not np.array_equal(dummy_adata.X, adata_norm.X)


def test_deseq2_normalize_sparse(sparse_adata):
    """Test DESeq2 normalization with sparse data."""
    adata_norm = deseq2_normalize(sparse_adata, min_genes=2)

    # Check that size factors are stored
    assert "size_factors" in adata_norm.obs.columns

    # Check that data is log1p transformed (convert to dense for comparison)
    if sparse.issparse(adata_norm.X):
        X_dense = adata_norm.X.toarray()
    else:
        X_dense = adata_norm.X
    assert np.all(X_dense >= 0)


def test_pydeseq2_normalize_global(dummy_adata):
    """Test PyDESeq2 normalization function."""
    adata_norm, size_factors = pydeseq2_normalize_global(dummy_adata, min_genes=5, apply_log1p=True)

    # Check that original data is not modified
    assert dummy_adata.X is not adata_norm.X

    # Check that size factors are stored
    assert "size_factors" in adata_norm.obs.columns
    assert "size_factors_global" in adata_norm.obs.columns

    # Check that data is log1p transformed
    assert np.all(adata_norm.X >= 0)

    # Check that size factors are returned
    assert len(size_factors) == dummy_adata.n_obs
    assert np.all(size_factors > 0)


def test_pydeseq2_normalize_global_no_log1p(dummy_adata):
    """Test PyDESeq2 normalization without log1p transformation."""
    adata_norm, size_factors = pydeseq2_normalize_global(
        dummy_adata, min_genes=5, apply_log1p=False
    )

    # Check that size factors are stored
    assert "size_factors" in adata_norm.obs.columns

    # Check that data is not log1p transformed (should have raw normalized values)
    # This is harder to test directly, but we can check the data type and range
    assert adata_norm.X.dtype in [np.float64, np.float32]


def test_apply_global_size_factors(dummy_adata):
    """Test applying global size factors to new data."""
    # First compute global size factors
    _, global_size_factors = pydeseq2_normalize_global(dummy_adata, min_genes=5)

    # Apply global size factors to the same data
    adata_norm = apply_global_size_factors(dummy_adata, global_size_factors)

    # Check that global size factors are stored
    assert "size_factors_global" in adata_norm.obs.columns

    # Check that normalization was applied
    assert not np.array_equal(dummy_adata.X, adata_norm.X)


def test_apply_global_size_factors_no_log1p(dummy_adata):
    """Test applying global size factors without log1p transformation."""
    # Create dummy global size factors
    global_size_factors = np.array([1.0, 1.2, 0.8, 1.1, 0.9])

    adata_norm = apply_global_size_factors(dummy_adata, global_size_factors, apply_log1p=False)

    # Check that global size factors are stored
    assert "size_factors_global" in adata_norm.obs.columns

    # Check that normalization was applied
    assert not np.array_equal(dummy_adata.X, adata_norm.X)


def test_prepare_for_de_analysis(dummy_adata):
    """Test preparing data for differential expression analysis."""
    # Create dummy size factors
    size_factors = np.array([1.0, 1.2, 0.8, 1.1, 0.9])

    adata_de = prepare_for_de_analysis(dummy_adata, size_factors)

    # Check that size factors are stored as offsets
    assert "size_factors_offset" in adata_de.obs.columns
    assert "log_size_factors" in adata_de.obs.columns

    # Check that log size factors are computed correctly
    expected_log_sf = np.log(size_factors)
    np.testing.assert_array_almost_equal(adata_de.obs["log_size_factors"].values, expected_log_sf)

    # Check that original data is preserved
    assert np.array_equal(dummy_adata.X, adata_de.X)


def test_plot_size_factor_distribution(dummy_adata):
    """Test plotting size factor distribution."""
    # First normalize the data to get size factors
    adata_norm = deseq2_normalize(dummy_adata, min_genes=5)

    # Test that plotting doesn't raise an error
    # Note: We can't easily test the actual plot output, but we can test that it runs
    try:
        plot_size_factor_distribution(adata_norm)
    except Exception as e:
        pytest.fail(f"plot_size_factor_distribution raised an exception: {e}")


def test_plot_size_factor_distribution_missing_factors(dummy_adata):
    """Test that plotting raises error when size factors are missing."""
    with pytest.raises(KeyError, match="Size factors not found"):
        plot_size_factor_distribution(dummy_adata)


def test_normalization_edge_cases():
    """Test normalization with edge cases."""
    # Test with very small data
    X = np.array([[1, 2], [3, 4]])
    adata = sc.AnnData(X=X)
    adata.var_names = ["Gene1", "Gene2"]
    adata.obs_names = ["Cell1", "Cell2"]

    # Should not raise an error
    adata_norm = deseq2_normalize(adata, min_genes=1)
    assert "size_factors" in adata_norm.obs.columns


def test_normalization_consistency():
    """Test that normalization is consistent across multiple runs."""
    np.random.seed(42)
    X = np.random.poisson(5, size=(10, 20)).astype(float)
    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(20)]
    adata.obs_names = [f"Cell_{i}" for i in range(10)]

    # Run normalization twice
    adata_norm1 = deseq2_normalize(adata, min_genes=5)
    adata_norm2 = deseq2_normalize(adata, min_genes=5)

    # Size factors should be the same
    np.testing.assert_array_almost_equal(
        adata_norm1.obs["size_factors"].values, adata_norm2.obs["size_factors"].values
    )

    # Normalized data should be the same
    np.testing.assert_array_almost_equal(adata_norm1.X, adata_norm2.X)


def test_size_factor_properties(dummy_adata):
    """Test that size factors have expected properties."""
    size_factors = compute_size_factors(dummy_adata, min_genes=5)

    # Size factors should be positive
    assert np.all(size_factors > 0)

    # Size factors should be centered around 1 (median should be close to 1)
    median_sf = np.median(size_factors)
    assert 0.5 < median_sf < 2.0  # Reasonable range for median

    # Size factors should not be all the same
    assert len(np.unique(size_factors)) > 1


def test_geometric_mean_properties(dummy_adata):
    """Test that geometric means have expected properties."""
    geometric_means = geometric_mean_per_gene(dummy_adata)

    # Geometric means should be non-negative
    assert np.all(geometric_means >= 0)

    # Should have some defined geometric means
    n_defined = np.sum(geometric_means > 0)
    assert n_defined > 0
    assert n_defined <= dummy_adata.n_vars
