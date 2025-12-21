"""
load_data.py

Loads normalized Xenium CSVs for plaque signature modeling.
Includes modality selection and validation.
"""

import pandas as pd
from pathlib import Path

from .schema import (
    SPATIAL_COLS,
    MORPH_COLS,
    CLUSTER_COL,
    TARGET_COL,
    PLAQUE_COLS,
    get_gene_columns,
)


def load_mouse_csv(path, include_modalities=None):
    """
    Load a mouse dataset and return a dict of feature matrices.

    Parameters
    ----------
    path : str or Path
        Path to the normalized CSV.

    include_modalities : list of str or None
        Modalities to include. Options:
        ["genes", "morph", "spatial", "cluster"]
        If None → all available modalities are returned.

    Returns
    -------
    data : dict
        {
            'df': original dataframe,
            'X': concatenated feature matrix (pd.DataFrame),
            'genes': gene-only dataframe,
            'morph': morphology dataframe,
            'spatial': spatial coords dataframe,
            'cluster': one-hot cluster dataframe,
            'target': distance_to_plaque (if exists),
            'modalities': list of modalities included
        }
    """

    path = Path(path)
    df = pd.read_csv(path)

    gene_cols = get_gene_columns(df)

    if include_modalities is None:
        include_modalities = ["genes", "morph", "spatial", "cluster"]

    data = {"df": df, "modalities": include_modalities}

    if "genes" in include_modalities:
        data["genes"] = df[gene_cols]
    else:
        data["genes"] = None

    if "morph" in include_modalities:
        data["morph"] = df[MORPH_COLS]
    else:
        data["morph"] = None

    if "spatial" in include_modalities:
        data["spatial"] = df[SPATIAL_COLS]
    else:
        data["spatial"] = None

    if "cluster" in include_modalities:

        clusters = pd.get_dummies(df[CLUSTER_COL], prefix="cluster")
        data["cluster"] = clusters
    else:
        data["cluster"] = None

    blocks = []
    for key in ["genes", "morph", "spatial", "cluster"]:
        if data[key] is not None:
            blocks.append(data[key])

    data["X"] = pd.concat(blocks, axis=1)

    if TARGET_COL in df.columns:
        data["target"] = df[TARGET_COL]
    else:
        data["target"] = None

    data["gene_cols"] = gene_cols

    return data
