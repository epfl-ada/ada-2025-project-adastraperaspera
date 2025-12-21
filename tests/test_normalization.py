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

    np.random.seed(42)
    X = np.random.poisson(5, size=(5, 10)).astype(float)

    X[0, 0] = 0
    X[1, 1] = 0

    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(10)]
    adata.obs_names = [f"Cell_{i}" for i in range(5)]

    return adata


@pytest.fixture
def sparse_adata():
    """Create a dummy sparse AnnData object for testing."""
    np.random.seed(42)
    X = np.random.poisson(3, size=(3, 5)).astype(float)
    X[X < 2] = 0

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

    assert len(geometric_means) == dummy_adata.n_vars
    assert geometric_means.shape == (10,)

    assert np.all(geometric_means >= 0)

    assert np.sum(geometric_means > 0) > 0


def test_geometric_mean_per_gene_sparse(sparse_adata):
    """Test geometric mean calculation with sparse data."""
    geometric_means = geometric_mean_per_gene(sparse_adata)

    assert len(geometric_means) == sparse_adata.n_vars
    assert geometric_means.shape == (5,)

    assert np.all(geometric_means >= 0)


def test_geometric_mean_per_gene_all_zeros(all_zeros_adata):
    """Test geometric mean calculation with all zeros."""
    geometric_means = geometric_mean_per_gene(all_zeros_adata)

    assert np.all(geometric_means == 0)


def test_compute_size_factors(dummy_adata):
    """Test size factor computation."""
    size_factors = compute_size_factors(dummy_adata, min_genes=5)

    assert len(size_factors) == dummy_adata.n_obs
    assert size_factors.shape == (5,)

    assert np.all(size_factors > 0)

    assert np.all(size_factors > 0.1)
    assert np.all(size_factors < 10.0)


def test_compute_size_factors_sparse(sparse_adata):
    """Test size factor computation with sparse data."""
    size_factors = compute_size_factors(sparse_adata, min_genes=2)

    assert len(size_factors) == sparse_adata.n_obs
    assert size_factors.shape == (3,)

    assert np.all(size_factors > 0)


def test_compute_size_factors_min_genes(dummy_adata):
    """Test size factor computation with different min_genes parameter."""

    size_factors = compute_size_factors(dummy_adata, min_genes=100)

    assert len(size_factors) == dummy_adata.n_obs
    assert np.all(size_factors > 0)


def test_deseq2_normalize(dummy_adata):
    """Test DESeq2 normalization function."""
    adata_norm = deseq2_normalize(dummy_adata, min_genes=5)

    assert dummy_adata.X is not adata_norm.X

    assert "size_factors" in adata_norm.obs.columns

    assert np.all(adata_norm.X >= 0)

    assert not np.array_equal(dummy_adata.X, adata_norm.X)


def test_deseq2_normalize_sparse(sparse_adata):
    """Test DESeq2 normalization with sparse data."""
    adata_norm = deseq2_normalize(sparse_adata, min_genes=2)

    assert "size_factors" in adata_norm.obs.columns

    if sparse.issparse(adata_norm.X):
        X_dense = adata_norm.X.toarray()
    else:
        X_dense = adata_norm.X
    assert np.all(X_dense >= 0)


def test_pydeseq2_normalize_global(dummy_adata):
    """Test PyDESeq2 normalization function."""
    adata_norm, size_factors = pydeseq2_normalize_global(
        dummy_adata, min_genes=5, apply_log1p=True
    )

    assert dummy_adata.X is not adata_norm.X

    assert "size_factors" in adata_norm.obs.columns
    assert "size_factors_global" in adata_norm.obs.columns

    assert np.all(adata_norm.X >= 0)

    assert len(size_factors) == dummy_adata.n_obs
    assert np.all(size_factors > 0)


def test_pydeseq2_normalize_global_no_log1p(dummy_adata):
    """Test PyDESeq2 normalization without log1p transformation."""
    adata_norm, size_factors = pydeseq2_normalize_global(
        dummy_adata, min_genes=5, apply_log1p=False
    )

    assert "size_factors" in adata_norm.obs.columns

    assert adata_norm.X.dtype in [np.float64, np.float32]


def test_apply_global_size_factors(dummy_adata):
    """Test applying global size factors to new data."""

    _, global_size_factors = pydeseq2_normalize_global(dummy_adata, min_genes=5)

    adata_norm = apply_global_size_factors(dummy_adata, global_size_factors)

    assert "size_factors_global" in adata_norm.obs.columns

    assert not np.array_equal(dummy_adata.X, adata_norm.X)


def test_apply_global_size_factors_no_log1p(dummy_adata):
    """Test applying global size factors without log1p transformation."""

    global_size_factors = np.array([1.0, 1.2, 0.8, 1.1, 0.9])

    adata_norm = apply_global_size_factors(
        dummy_adata, global_size_factors, apply_log1p=False
    )

    assert "size_factors_global" in adata_norm.obs.columns

    assert not np.array_equal(dummy_adata.X, adata_norm.X)


def test_prepare_for_de_analysis(dummy_adata):
    """Test preparing data for differential expression analysis."""

    size_factors = np.array([1.0, 1.2, 0.8, 1.1, 0.9])

    adata_de = prepare_for_de_analysis(dummy_adata, size_factors)

    assert "size_factors_offset" in adata_de.obs.columns
    assert "log_size_factors" in adata_de.obs.columns

    expected_log_sf = np.log(size_factors)
    np.testing.assert_array_almost_equal(
        adata_de.obs["log_size_factors"].values, expected_log_sf
    )

    assert np.array_equal(dummy_adata.X, adata_de.X)


def test_plot_size_factor_distribution(dummy_adata):
    """Test plotting size factor distribution."""

    adata_norm = deseq2_normalize(dummy_adata, min_genes=5)

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

    X = np.array([[1, 2], [3, 4]])
    adata = sc.AnnData(X=X)
    adata.var_names = ["Gene1", "Gene2"]
    adata.obs_names = ["Cell1", "Cell2"]

    adata_norm = deseq2_normalize(adata, min_genes=1)
    assert "size_factors" in adata_norm.obs.columns


def test_normalization_consistency():
    """Test that normalization is consistent across multiple runs."""
    np.random.seed(42)
    X = np.random.poisson(5, size=(10, 20)).astype(float)
    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(20)]
    adata.obs_names = [f"Cell_{i}" for i in range(10)]

    adata_norm1 = deseq2_normalize(adata, min_genes=5)
    adata_norm2 = deseq2_normalize(adata, min_genes=5)

    np.testing.assert_array_almost_equal(
        adata_norm1.obs["size_factors"].values, adata_norm2.obs["size_factors"].values
    )

    np.testing.assert_array_almost_equal(adata_norm1.X, adata_norm2.X)


def test_size_factor_properties(dummy_adata):
    """Test that size factors have expected properties."""
    size_factors = compute_size_factors(dummy_adata, min_genes=5)

    assert np.all(size_factors > 0)

    median_sf = np.median(size_factors)
    assert 0.5 < median_sf < 2.0

    assert len(np.unique(size_factors)) > 1


def test_geometric_mean_properties(dummy_adata):
    """Test that geometric means have expected properties."""
    geometric_means = geometric_mean_per_gene(dummy_adata)

    assert np.all(geometric_means >= 0)

    n_defined = np.sum(geometric_means > 0)
    assert n_defined > 0
    assert n_defined <= dummy_adata.n_vars
