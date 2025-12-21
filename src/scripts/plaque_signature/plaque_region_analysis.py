import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import seaborn as sns


METADATA_COLS = [
    "cell_id",
    "x_centroid",
    "y_centroid",
    "transcript_counts",
    "control_probe_counts",
    "control_codeword_counts",
    "unassigned_codeword_counts",
    "total_counts",
    "cell_area",
    "nucleus_area",
    "cluster_leiden",
    "distance_to_plaque",
    "nearest_plaque_center_dist",
    "inside_any_plaque",
    "nearest_plaque_id",
    "nearest_plaque_area",
    "plaque_region",
]


def extract_genes(df):
    gene_cols = [
        c for c in df.columns if c not in METADATA_COLS and df[c].dtype != "object"
    ]
    return df[gene_cols].values, gene_cols


def train_pls_on_tg(df_tg17, gene_cols, n_components=10):
    X = df_tg17[gene_cols].values
    y = df_tg17["distance_to_plaque"].values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    pls = PLSRegression(n_components=n_components)
    pls.fit(X_train_scaled, y_train)

    train_r2 = r2_score(y_train, pls.predict(X_train_scaled))
    test_r2 = r2_score(y_test, pls.predict(X_test_scaled))

    return pls, scaler, train_r2, test_r2


def apply_pls_signature(df, gene_cols, scaler, pls):
    X = df[gene_cols].values
    X_scaled = scaler.transform(X)
    latent = pls.transform(X_scaled)
    sig = latent[:, 0]
    sig_norm = (sig - sig.mean()) / sig.std()
    df["pls_signature"] = sig_norm
    return df


def compute_plaque_centers(df_tg17, plaque_boxes_df):
    centers = {}
    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]
        pr = df_tg17[
            (df_tg17["nearest_plaque_id"] == pid)
            & (df_tg17["distance_to_plaque"] < 100)
        ]
        cx = pr["x_centroid"].mean()
        cy = pr["y_centroid"].mean()
        centers[pid] = (cx, cy)
    return centers


def compute_wt_distances(df_wt, plaque_boxes_df, plaque_centers):
    records = []
    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]
        cx, cy = plaque_centers[pid]

        wt_patch = df_wt[
            (df_wt["x_centroid"] >= row["x_min"])
            & (df_wt["x_centroid"] <= row["x_max"])
            & (df_wt["y_centroid"] >= row["y_min"])
            & (df_wt["y_centroid"] <= row["y_max"])
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


def compute_tg_distances(df_tg, plaque_boxes_df):
    records = []
    for _, row in plaque_boxes_df.iterrows():
        pid = row["plaque_id"]
        tg_patch = df_tg[
            (df_tg["nearest_plaque_id"] == pid) & (df_tg["distance_to_plaque"] < 1000)
        ]
        records.append(
            {
                "plaque_id": pid,
                "tg_distances": tg_patch["distance_to_plaque"].values,
                "tg_signatures": tg_patch["pls_signature"].values,
            }
        )
    return records


def plot_decay(pid, tg_records, wt_records):
    tg = next(r for r in tg_records if r["plaque_id"] == pid)
    wt = next(r for r in wt_records if r["plaque_id"] == pid)

    plt.figure(figsize=(7, 4))
    sns.scatterplot(
        x=tg["tg_distances"], y=tg["tg_signatures"], s=5, alpha=0.4, label="Tg17"
    )
    sns.scatterplot(
        x=wt["wt_distances"], y=wt["wt_signatures"], s=5, alpha=0.4, label="WT13"
    )
    plt.xlabel("Distance to Plaque (µm)")
    plt.ylabel("PLS Signature")
    plt.title(f"Decay Curve – Plaque {pid}")
    plt.legend()
    plt.show()


def plot_binned_profile(pid, tg_records, wt_records, bin_size=20, max_dist=300):
    tg = next(r for r in tg_records if r["plaque_id"] == pid)
    wt = next(r for r in wt_records if r["plaque_id"] == pid)

    tg_df = pd.DataFrame({"dist": tg["tg_distances"], "sig": tg["tg_signatures"]})
    wt_df = pd.DataFrame({"dist": wt["wt_distances"], "sig": wt["wt_signatures"]})

    bins = np.arange(0, max_dist + bin_size, bin_size)

    tg_binned = tg_df.groupby(pd.cut(tg_df["dist"], bins))["sig"].agg(["mean", "sem"])
    wt_binned = wt_df.groupby(pd.cut(wt_df["dist"], bins))["sig"].agg(["mean", "sem"])

    centers = bins[:-1] + bin_size / 2

    plt.figure(figsize=(7, 5))
    plt.errorbar(
        centers, tg_binned["mean"], yerr=tg_binned["sem"], label="Tg17", color="red"
    )
    plt.errorbar(
        centers, wt_binned["mean"], yerr=wt_binned["sem"], label="WT13", color="blue"
    )
    plt.title(f"Binned Signature Profile – Plaque {pid}")
    plt.xlabel("Distance to Plaque (µm)")
    plt.ylabel("PLS Signature")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()
