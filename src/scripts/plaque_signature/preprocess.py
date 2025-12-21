"""
preprocess.py

Preprocessing utilities for plaque signature modeling:
- scaling
- train/test splitting
- modality normalization
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from .schema import TARGET_COL


def split_train_test(df, test_fraction=0.2, random_state=42):
    """
    Only used for Tg 17.9-month mouse (where plaque distance exists).
    Returns train/test indices.
    """
    if TARGET_COL not in df.columns:
        raise ValueError(
            "Target column not found - this mouse has no plaque distances."
        )

    y = df[TARGET_COL]
    idx = np.arange(len(df))

    train_idx, test_idx = train_test_split(
        idx, test_size=test_fraction, random_state=random_state, shuffle=True
    )

    return train_idx, test_idx


class ModalityScaler:
    """
    Fits separate scalers for each modality (genes, morph, spatial, cluster),
    then concatenates scaled data into a single feature matrix.
    """

    def __init__(self):
        self.scalers = {}

    def fit(self, data):
        """
        data: dict from load_mouse_csv()
        """
        for modality in ["genes", "morph", "spatial", "cluster"]:
            block = data.get(modality)
            if block is None:
                continue
            scaler = StandardScaler()
            scaler.fit(block.values)
            self.scalers[modality] = scaler

    def transform(self, data):
        """
        Returns a scaled feature matrix aligned with the fitted modalities.
        """
        blocks = []
        for modality in ["genes", "morph", "spatial", "cluster"]:
            block = data.get(modality)
            if block is None:
                continue
            scaler = self.scalers[modality]
            X_scaled = scaler.transform(block.values)
            blocks.append(
                pd.DataFrame(
                    X_scaled,
                    index=block.index,
                    columns=[f"{modality}_{col}" for col in block.columns],
                )
            )
        return pd.concat(blocks, axis=1)
