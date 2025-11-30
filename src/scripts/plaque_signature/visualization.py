import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import umap
from scipy.spatial.distance import jensenshannon


def rotate_coords(df, angle_deg, x_col="x_centroid", y_col="y_centroid"):
    """
    Rotate coordinates in a dataframe by angle_deg degrees around their centroid.
    Positive = counterclockwise rotation.

    Returns a copy with rotated x and y.
    """
    df_rot = df.copy()

    angle = np.radians(angle_deg)

    # Center of rotation = mean position
    cx = df_rot[x_col].mean()
    cy = df_rot[y_col].mean()

    # Shift points so centroid becomes (0,0)
    x = df_rot[x_col] - cx
    y = df_rot[y_col] - cy

    x_new = x * np.cos(angle) - y * np.sin(angle)
    y_new = x * np.sin(angle) + y * np.cos(angle)

    # Shift back to original location
    df_rot[x_col] = x_new + cx
    df_rot[y_col] = y_new + cy

    return df_rot

def plot_spatial_mouse(all_preds, name, cmap="viridis"):
    """
    Scatterplot of tissue coordinates colored by predicted plaque score.
    """
    df = all_preds[name]["df"]
    pred = all_preds[name]["pred"]

    plt.figure(figsize=(7, 7))
    plt.scatter(
        df["x_centroid"], df["y_centroid"],
        c=pred, s=4, cmap=cmap, alpha=0.6
    )
    plt.gca().invert_yaxis()
    plt.colorbar(label="Predicted Plaque Proximity (HGB Spatial-Only)")
    plt.title(f"Spatial Prediction Map — {name}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.tight_layout()
    plt.show()


def mirror_coords(df, mode="horizontal", x_col="x_centroid", y_col="y_centroid"):
    """
    Mirror coordinates in a dataframe around their centroid.
    
    mode:
        - "horizontal": flips x
        - "vertical": flips y
        - "both": flips both x and y
        - None or "none": do nothing

    Returns a copy of the dataframe with mirrored coords.
    """
    df_m = df.copy()

    cx = df_m[x_col].mean()
    cy = df_m[y_col].mean()

    if mode == "horizontal":
        df_m[x_col] = 2*cx - df_m[x_col]
    elif mode == "vertical":
        df_m[y_col] = 2*cy - df_m[y_col]
    elif mode == "both":
        df_m[x_col] = 2*cx - df_m[x_col]
        df_m[y_col] = 2*cy - df_m[y_col]
    else:
        # no mirroring
        return df

    return df_m


def plot_spatial_compare(all_preds, names, rotations=None, mirrors=None, cmap="viridis"):
    """
    Side-by-side spatial maps for multiple mice.

    rotations : dict name → angle_deg
    mirrors   : dict name → {"horizontal", "vertical", "both", None}
    """

    if rotations is None:
        rotations = {name: 0 for name in names}

    if mirrors is None:
        mirrors = {name: None for name in names}

    fig, axes = plt.subplots(1, len(names), figsize=(7 * len(names), 7))

    for ax, name in zip(axes, names):
        df = all_preds[name]["df"]
        pred = all_preds[name]["pred"]

        # 1. mirror first 
        mirror_mode = mirrors.get(name, None)
        df_plot = mirror_coords(df, mirror_mode) if mirror_mode else df

        # 2. rotate after mirroring
        angle = rotations.get(name, 0)
        df_plot = rotate_coords(df_plot, angle) if angle != 0 else df_plot

        sc = ax.scatter(
            df_plot["x_centroid"], df_plot["y_centroid"],
            c=pred, s=4, cmap=cmap, alpha=0.6
        )
        ax.invert_yaxis()
        ax.set_title(f"{name} (rot={angle}°, mirror={mirror_mode})")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

    fig.colorbar(sc, ax=axes.ravel().tolist(), label="Predicted Plaque Score")
    plt.show()


def plot_hist_comparison(all_preds, names):
    """
    Compare score distributions across multiple mice.
    """
    plt.figure(figsize=(10, 6))

    for name in names:
        preds = all_preds[name]["pred"]
        sns.kdeplot(preds, label=name, fill=True, alpha=0.3)

    plt.title("Distribution of Predicted Plaque Score Across Mice")
    plt.xlabel("Predicted Score (HGB Spatial-Only)")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.show()


def compute_umap(df, n_neighbors=30, min_dist=0.1):
    """
    Compute UMAP embedding from spatial coordinates.
    Returns an array of shape (n_cells, 2).
    """
    coords = df[["x_centroid", "y_centroid"]].values

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric="euclidean"
    )
    emb = reducer.fit_transform(coords)
    return emb

def plot_umap_embedding(df, emb, color_values=None, cmap="viridis",
                        title="UMAP Embedding", s=4, alpha=0.6):
    """
    Plot a precomputed UMAP embedding with optional coloring.
    """
    plt.figure(figsize=(7, 7))

    sc = plt.scatter(
        emb[:, 0], emb[:, 1],
        c=color_values, cmap=cmap, s=s, alpha=alpha
    )

    if color_values is not None:
        plt.colorbar(sc, label="Color")

    plt.xlabel("UMAP-1")
    plt.ylabel("UMAP-2")
    plt.title(title)
    plt.tight_layout()
    plt.show()


def compute_and_plot_umap(all_preds, name,
                          n_neighbors=30, min_dist=0.1,
                          color="prediction", cmap="viridis"):
    """
    Fast helper: compute UMAP once, then call plot_umap_embedding.
    """

    df = all_preds[name]["df"]
    pred = all_preds[name]["pred"]

    emb = compute_umap(df, n_neighbors=n_neighbors, min_dist=min_dist)

    if color == "prediction":
        cvals = pred
        title = f"{name} — UMAP Colored by Prediction"
    elif color == "x":
        cvals = df["x_centroid"]
        title = f"{name} — UMAP Colored by X"
    elif color == "y":
        cvals = df["y_centroid"]
        title = f"{name} — UMAP Colored by Y"
    elif color == "none":
        cvals = None
        title = f"{name} — UMAP (No color)"
    else:
        raise ValueError(f"Unknown color option {color}")

    plot_umap_embedding(df, emb, cvals, cmap=cmap, title=title)

    return emb

def plot_ablation_heatmap(results_ablation, cmap="viridis"):
    """
    Heatmap of test R² for modality ablation results.
    """
    heatmap_df = results_ablation.pivot_table(
        index="Modalities",
        columns="Model",
        values="test_r2"
    )

    plt.figure(figsize=(16, 8))
    sns.heatmap(
        heatmap_df,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        linewidths=0.5,
        cbar_kws={"label": "Test R²"}
    )
    plt.title("Modality Ablation Heatmap (Test R²)", fontsize=16)
    plt.xlabel("Model")
    plt.ylabel("Modalities")
    plt.tight_layout()
    plt.show()


def plot_best_model_per_modality(results_ablation):
    """
    Barplot of the best-performing model for each modality.
    """
    best_per_modality = (
        results_ablation
        .sort_values(["Modalities", "test_r2"], ascending=[True, False])
        .groupby("Modalities")
        .first()
        .reset_index()
    )

    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=best_per_modality,
        x="Modalities",
        y="test_r2",
        hue="Model",
        palette="tab10"
    )

    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Best Test R²")
    plt.title("Best Model for Each Modality")
    plt.tight_layout()
    plt.show()


def plot_pred_vs_true(y_true, y_pred, title="Predicted vs True Plaque Distance"):
    """
    Scatter plot for predicted vs. true plaque distances.
    """
    plt.figure(figsize=(7, 7))
    plt.scatter(y_true, y_pred, s=10, alpha=0.4)
    plt.plot(
        [min(y_true), max(y_true)],
        [min(y_true), max(y_true)],
        "r--"
    )
    plt.xlabel("True Distance")
    plt.ylabel("Predicted Distance")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_spatial_score_vs_age_by_genotype(all_preds, age_map, disease_map,
                                          title="Spatial-Only Predicted Score vs Age (TG vs WT)"):
    """
    Plot mean ± std predicted plaque score vs age for TG and WT mice.

    Parameters:
    -----------
    all_preds : dict
        Output from predict_spatial_only().

    age_map : dict
        Example: {"tg2": 2, "tg5": 5, "tg17": 17, "wt2": 2, "wt5": 5, "wt13": 13}

    disease_map : dict
        Maps mouse name → "TG" or "WT".

    title : str
        Title for the plot.
    """
    # Build grouped values
    grouped = {}

    for name, entry in all_preds.items():
        age = age_map[name]
        disease = disease_map[name]
        pred = entry["pred"]

        if disease not in grouped:
            grouped[disease] = {"ages": [], "means": [], "stds": []}

        grouped[disease]["ages"].append(age)
        grouped[disease]["means"].append(pred.mean())
        grouped[disease]["stds"].append(pred.std())

    # Plot
    plt.figure(figsize=(9, 6))

    for disease, d in grouped.items():
        # Sort by age
        ages, means, stds = zip(*sorted(zip(d["ages"], d["means"], d["stds"])))
        plt.errorbar(
            ages, means, yerr=stds,
            fmt="-o", capsize=5,
            label=disease
        )

    plt.title(title)
    plt.xlabel("Age (months)")
    plt.ylabel("Mean Predicted Plaque Score")
    plt.legend(title="Genotype")
    plt.tight_layout()
    plt.show()


def compute_jsd_matrix(all_preds, bins=100):
    """
    Compute a Jensen-Shannon divergence matrix between mice,
    based on predicted score distributions.

    Parameters
    ----------
    all_preds : dict
        From predict_spatial_only(): {mouse: {"df":..., "pred":...}}
    bins : int
        Number of bins for histogram JSD approximation.

    Returns
    -------
    DataFrame
        Square matrix of JSD distances.
    """

    mice = list(all_preds.keys())
    jsd_matrix = pd.DataFrame(index=mice, columns=mice)

    def jsd(p, q):
        min_v, max_v = min(p.min(), q.min()), max(p.max(), q.max())
        hist_p, _ = np.histogram(p, bins=bins, range=(min_v, max_v), density=True)
        hist_q, _ = np.histogram(q, bins=bins, range=(min_v, max_v), density=True)
        return jensenshannon(hist_p, hist_q)

    for m1 in mice:
        p1 = all_preds[m1]["pred"]
        for m2 in mice:
            p2 = all_preds[m2]["pred"]
            jsd_matrix.loc[m1, m2] = jsd(p1, p2)

    return jsd_matrix.astype(float)


def plot_jsd_heatmap(jsd_matrix, title="Jensen–Shannon Divergence Between Mice"):
    """
    Plot a heatmap of the Jensen–Shannon divergence matrix.

    Parameters
    ----------
    jsd_matrix : DataFrame
        Square matrix of JSD distances (output of compute_jsd_matrix).
    """

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        jsd_matrix,
        annot=True,
        fmt=".3f",
        cmap="mako",  # darker → lighter
        linewidths=0.5,
        square=True,
        cbar_kws={"label": "Jensen–Shannon Distance"}
    )
    plt.title(title, fontsize=14)
    plt.xlabel("Mouse")
    plt.ylabel("Mouse")
    plt.tight_layout()
    plt.show()

