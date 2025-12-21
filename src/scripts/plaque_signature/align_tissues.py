"""
align_tissues.py

Utilities for spatial alignment of WT and Tg Xenium tissues.
Includes:
- Loading files
- Coordinate range inspection
- Flipping WT horizontally
- Centroid-based alignment
- Overlay plotting
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_mouse(tg17_path, wt13_path):
    df_tg17 = pd.read_csv(tg17_path)
    df_wt13 = pd.read_csv(wt13_path)
    return df_tg17, df_wt13


def coord_stats(df, name="mouse"):
    print(f"\n{name} coordinate ranges:")
    print("  x min/max:", df["x_centroid"].min(), df["x_centroid"].max())
    print("  y min/max:", df["y_centroid"].min(), df["y_centroid"].max())


def flip_x(df):
    df_flip = df.copy()
    df_flip["x_centroid"] = df_flip["x_centroid"].max() - df_flip["x_centroid"]
    return df_flip


def compute_translation(df_ref, df_to_align):
    cx_ref = df_ref["x_centroid"].median()
    cy_ref = df_ref["y_centroid"].median()

    cx_align = df_to_align["x_centroid"].median()
    cy_align = df_to_align["y_centroid"].median()

    dx = cx_ref - cx_align
    dy = cy_ref - cy_align

    print("Reference centroid:", cx_ref, cy_ref)
    print("To-align centroid:", cx_align, cy_align)
    print("Translation:", dx, dy)

    return dx, dy


def translate(df, dx, dy):
    df_aligned = df.copy()
    df_aligned["x_centroid"] = df_aligned["x_centroid"] + dx
    df_aligned["y_centroid"] = df_aligned["y_centroid"] + dy
    return df_aligned


def plot_overlay(df_ref, df_aligned, label_ref="Tg17", label_aligned="WT13"):
    plt.figure(figsize=(8, 8))
    plt.scatter(
        df_ref["x_centroid"], df_ref["y_centroid"], s=1, alpha=0.3, label=label_ref
    )
    plt.scatter(
        df_aligned["x_centroid"],
        df_aligned["y_centroid"],
        s=1,
        alpha=0.3,
        label=label_aligned,
    )
    plt.gca().invert_yaxis()
    plt.legend()
    plt.title("Aligned Tissues (flip-X + translation)")
    plt.show()
