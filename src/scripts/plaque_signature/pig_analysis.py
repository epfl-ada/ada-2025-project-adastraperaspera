import pandas as pd


def load_all_mice(mouse_paths):
    """
    Load all mouse CSVs into a dict of DataFrames.
    Ensures cluster_leiden is integer and drops missing clusters.
    """
    mice = {}
    for name, path in mouse_paths.items():
        df = pd.read_csv(path)
        df = df.dropna(subset=["cluster_leiden"]).copy()
        df["cluster_leiden"] = df["cluster_leiden"].astype(int)
        mice[name] = df
    return mice


def zscore_mouse(df, gene_cols):
    """
    Z-score all gene columns within a single mouse.
    """
    df_z = df.copy()
    for g in gene_cols:
        mu = df[g].mean()
        sd = df[g].std()
        df_z[g] = 0 if sd == 0 else (df[g] - mu) / sd
    return df_z


def zscore_all_mice(mice, gene_cols):
    """
    Apply z-scoring to each mouse independently.
    """
    mice_z = {}
    for name, df in mice.items():
        mice_z[name] = zscore_mouse(df, gene_cols)
    return mice_z


import numpy as np

PIG_GENES = [
    "Hexb",
    "Cst3",
    "Cd63",
    "C4b",
    "Ctsd",
    "B2m",
    "H2-K1",
    "Apoe",
    "Gfap",
    "Nrep",
    "Serpina3n",
    "Cd74",
    "Cxcl10",
    "Vim",
    "S100a6",
    "Ifit3",
]


def compute_pig_scores(mice_z, pig_genes=PIG_GENES):
    """
    Compute PIG (Plaque-Induced Gene) scores per cluster per mouse.

    Returns:
        DataFrame with columns:
        mouse, cluster, mean_pig, std_pig, cells
    """
    rows = []

    for mouse, df in mice_z.items():
        for c in sorted(df["cluster_leiden"].unique()):
            df_c = df[df["cluster_leiden"] == c]

            pig_score = df_c[pig_genes].mean(axis=1)

            rows.append(
                {
                    "mouse": mouse,
                    "cluster": c,
                    "mean_pig": pig_score.mean(),
                    "std_pig": pig_score.std(),
                    "cells": len(pig_score),
                }
            )

    return pd.DataFrame(rows)
