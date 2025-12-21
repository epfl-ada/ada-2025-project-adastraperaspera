"""
plaque_decay_analysis.py

This module analyzes spatial decay of the PLS plaque signature
for Tg17 vs WT13 using:
- plaque bounding boxes (from plaque_region_analysis)
- plaque centers (Tg-based)
- PLS signatures (already computed externally)

It includes:
- Tg and WT distance extraction
- decay plotting utilities
- plaque-wise correlation & slope statistics
- FDR correction
- human-readable interpretation labels
- confidence scoring
- global summary stats
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests


def compute_tg_distances(df_tg17, plaque_boxes_df):
    """
    Extract Tg17 true distances + PLS signatures per plaque.
    """
    records = []
    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]

        tg_patch = df_tg17[
            (df_tg17["nearest_plaque_id"] == pid)
            & (df_tg17["distance_to_plaque"] < 1000)
        ]

        records.append(
            {
                "plaque_id": pid,
                "tg_distances": tg_patch["distance_to_plaque"].values,
                "tg_signatures": tg_patch["pls_signature"].values,
            }
        )

    return records


def compute_wt_distances(df_wt_aligned, plaque_boxes_df, plaque_centers):
    """
    Estimate WT13 distances from Tg plaque centers inside bbox.
    """
    records = []
    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]
        cx, cy = plaque_centers[pid]

        wt_patch = df_wt_aligned[
            (df_wt_aligned["x_centroid"] >= row["x_min"])
            & (df_wt_aligned["x_centroid"] <= row["x_max"])
            & (df_wt_aligned["y_centroid"] >= row["y_min"])
            & (df_wt_aligned["y_centroid"] <= row["y_max"])
        ].copy()

        wt_patch["est_dist"] = np.sqrt(
            (wt_patch["x_centroid"] - cx) ** 2 + (wt_patch["y_centroid"] - cy) ** 2
        )

        records.append(
            {
                "plaque_id": pid,
                "wt_distances": wt_patch["est_dist"].values,
                "wt_signatures": wt_patch["pls_signature"].values,
            }
        )

    return records


def plot_decay(pid, tg_records, wt_records):
    """
    Raw scatter comparison of Tg vs WT signature decay.
    """
    tg = next(r for r in tg_records if r["plaque_id"] == pid)
    wt = next(r for r in wt_records if r["plaque_id"] == pid)

    plt.figure(figsize=(7, 4))
    sns.scatterplot(
        x=tg["tg_distances"], y=tg["tg_signatures"], s=6, alpha=0.35, label="Tg17"
    )
    sns.scatterplot(
        x=wt["wt_distances"], y=wt["wt_signatures"], s=6, alpha=0.35, label="WT13"
    )
    plt.xlabel("Distance to Plaque (µm)")
    plt.ylabel("PLS Signature")
    plt.title(f"Signature Decay — Plaque {pid}")
    plt.legend()
    plt.show()


def plot_binned_profile(pid, tg_records, wt_records, bin_size=20, max_dist=300):
    """
    Binned distance profiles with SEM error bars.
    """
    tg = next(r for r in tg_records if r["plaque_id"] == pid)
    wt = next(r for r in wt_records if r["plaque_id"] == pid)

    tg_df = pd.DataFrame({"dist": tg["tg_distances"], "sig": tg["tg_signatures"]})
    wt_df = pd.DataFrame({"dist": wt["wt_distances"], "sig": wt["wt_signatures"]})

    bins = np.arange(0, max_dist + bin_size, bin_size)

    tg_b = tg_df.groupby(pd.cut(tg_df["dist"], bins))["sig"].agg(["mean", "sem"])
    wt_b = wt_df.groupby(pd.cut(wt_df["dist"], bins))["sig"].agg(["mean", "sem"])

    centers = bins[:-1] + bin_size / 2

    plt.figure(figsize=(7, 5))
    plt.errorbar(centers, tg_b["mean"], yerr=tg_b["sem"], label="Tg17", color="red")
    plt.errorbar(centers, wt_b["mean"], yerr=wt_b["sem"], label="WT13", color="blue")

    plt.xlabel("Distance (µm)")
    plt.ylabel("PLS Signature")
    plt.title(f"Binned Signature Profile — Plaque {pid}")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()


def compute_stats(plaque_boxes_df, tg_records, wt_records):
    """
    Computes:
        - Spearman correlations (Tg, WT)
        - OLS slopes (Tg, WT)
        - inner-zone (0–50 µm) mean signature differences
        - FDR-corrected significance flags
    """
    rows = []

    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]

        tg = next(r for r in tg_records if r["plaque_id"] == pid)
        wt = next(r for r in wt_records if r["plaque_id"] == pid)

        tg_dist, tg_sig = tg["tg_distances"], tg["tg_signatures"]
        wt_dist, wt_sig = wt["wt_distances"], wt["wt_signatures"]

        tg_corr, tg_p = spearmanr(tg_dist, tg_sig)
        wt_corr, wt_p = spearmanr(wt_dist, wt_sig)

        tg_slope = LinearRegression().fit(tg_dist.reshape(-1, 1), tg_sig).coef_[0]
        wt_slope = LinearRegression().fit(wt_dist.reshape(-1, 1), wt_sig).coef_[0]

        tg_inner = tg_sig[tg_dist < 50].mean() if np.any(tg_dist < 50) else np.nan
        wt_inner = wt_sig[wt_dist < 50].mean() if np.any(wt_dist < 50) else np.nan
        inner_diff = tg_inner - wt_inner

        rows.append(
            {
                "plaque_id": pid,
                "tg_corr": tg_corr,
                "tg_corr_p": tg_p,
                "wt_corr": wt_corr,
                "wt_corr_p": wt_p,
                "tg_slope": tg_slope,
                "wt_slope": wt_slope,
                "inner_mean_diff": inner_diff,
            }
        )

    stats_df = pd.DataFrame(rows)

    stats_df["tg_corr_fdr"] = multipletests(stats_df["tg_corr_p"], method="fdr_bh")[1]
    stats_df["wt_corr_fdr"] = multipletests(stats_df["wt_corr_p"], method="fdr_bh")[1]
    stats_df["tg_sig"] = stats_df["tg_corr_fdr"] < 0.05
    stats_df["wt_sig"] = stats_df["wt_corr_fdr"] < 0.05

    return stats_df


def add_significance_and_labels(stats_df):
    """
    Adds:
        - slope class (strong_negative, weak_positive, flat, etc)
        - interpretation text
        - confidence score
    """

    def slope_label(s):
        if s < -0.01:
            return "strong_negative"
        elif s < -0.002:
            return "weak_negative"
        elif s > 0.01:
            return "strong_positive"
        elif s > 0.002:
            return "weak_positive"
        return "flat"

    stats_df["tg_slope_label"] = stats_df["tg_slope"].apply(slope_label)
    stats_df["wt_slope_label"] = stats_df["wt_slope"].apply(slope_label)

    def interpret(r):
        if r["tg_sig"] and not r["wt_sig"]:
            return "Tg shows significant spatial decay; WT shows no structure."
        if r["tg_sig"] and r["wt_sig"]:
            return "Both Tg and WT show spatial structure (likely anatomical or layer mismatch)."
        if not r["tg_sig"] and r["wt_sig"]:
            return (
                "WT shows structure; likely misalignment or cortical layer difference."
            )
        return (
            "Neither Tg nor WT shows strong decay (weak plaque region or noisy area)."
        )

    stats_df["interpretation"] = stats_df.apply(interpret, axis=1)

    def score(r):
        val = 0
        if r["tg_sig"]:
            val += 2
        val += abs(r["tg_slope"]) * 50
        val += abs(r["tg_corr"])
        val += abs(r["inner_mean_diff"]) * 0.2
        return val

    stats_df["confidence"] = stats_df.apply(score, axis=1)

    return stats_df


def summarize_global(stats_df):
    """
    Provides a high-level global summary across all plaques.
    """
    return {
        "mean_tg_slope": stats_df["tg_slope"].mean(),
        "mean_wt_slope": stats_df["wt_slope"].mean(),
        "mean_slope_diff": (stats_df["tg_slope"] - stats_df["wt_slope"]).mean(),
        "mean_tg_corr": stats_df["tg_corr"].mean(),
        "mean_wt_corr": stats_df["wt_corr"].mean(),
        "mean_inner_diff": stats_df["inner_mean_diff"].mean(),
    }
