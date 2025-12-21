"""
plaque_regions.py

Functions for:
- Creating a plaque mask (distance threshold)
- Summarizing plaques based on nearest_plaque_id
- Selecting large / meaningful plaques
- Extracting bounding boxes for each plaque
- Extracting matching WT regions after alignment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def create_plaque_mask(df_tg, dist_threshold=50):
    df_tg = df_tg.copy()
    df_tg["plaque_region"] = df_tg["distance_to_plaque"] < dist_threshold
    return df_tg


def plot_plaque_map(df_tg, dist_threshold=50):
    df_temp = create_plaque_mask(df_tg, dist_threshold)

    plt.figure(figsize=(7, 7))
    plt.scatter(
        df_temp["x_centroid"],
        df_temp["y_centroid"],
        c=df_temp["plaque_region"],
        cmap="coolwarm",
        s=3,
        alpha=0.5,
    )
    plt.gca().invert_yaxis()
    plt.title(f"Tg17 Plaque Region Map (<{dist_threshold} µm)")
    plt.colorbar(label="Inside Plaque Region")
    plt.show()


def summarize_plaques(df_tg, local_dist=100):
    plaque_ids = df_tg["nearest_plaque_id"].unique()
    plaque_ids = [pid for pid in plaque_ids if pid != -1]

    plaque_summary = []

    for pid in plaque_ids:
        mask = (df_tg["nearest_plaque_id"] == pid) & (
            df_tg["distance_to_plaque"] < local_dist
        )
        count = mask.sum()
        plaque_summary.append({"plaque_id": pid, "local_cells": count})

    plaque_summary = pd.DataFrame(plaque_summary)
    plaque_summary = plaque_summary.sort_values("local_cells", ascending=False)

    return plaque_summary


def select_big_plaques(plaque_summary, min_cells=100):
    big_plaques = plaque_summary[plaque_summary["local_cells"] >= min_cells][
        "plaque_id"
    ].tolist()
    return big_plaques


def plaque_bounding_boxes(df_tg, big_plaques, local_dist=100):
    boxes = []

    for pid in big_plaques:
        pr = df_tg[
            (df_tg["nearest_plaque_id"] == pid)
            & (df_tg["distance_to_plaque"] < local_dist)
        ]

        x_min = pr["x_centroid"].min()
        x_max = pr["x_centroid"].max()
        y_min = pr["y_centroid"].min()
        y_max = pr["y_centroid"].max()

        boxes.append(
            {
                "plaque_id": pid,
                "x_min": x_min,
                "x_max": x_max,
                "y_min": y_min,
                "y_max": y_max,
                "cells": len(pr),
            }
        )

    return pd.DataFrame(boxes).sort_values("cells", ascending=False)


def extract_wt_regions(df_wt_aligned, plaque_boxes_df):
    wt_regions = []

    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]
        x_min, x_max = row["x_min"], row["x_max"]
        y_min, y_max = row["y_min"], row["y_max"]

        wt_patch = df_wt_aligned[
            (df_wt_aligned["x_centroid"] >= x_min)
            & (df_wt_aligned["x_centroid"] <= x_max)
            & (df_wt_aligned["y_centroid"] >= y_min)
            & (df_wt_aligned["y_centroid"] <= y_max)
        ]

        wt_regions.append({"plaque_id": pid, "wt_cells": len(wt_patch)})

    return pd.DataFrame(wt_regions)
