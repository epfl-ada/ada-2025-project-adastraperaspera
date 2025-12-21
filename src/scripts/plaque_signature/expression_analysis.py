import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def run_pca_per_cluster(mice, gene_cols, glial_clusters, disease_map):
    """
    Compute PCA for each glial cluster using gene expression.
    Returns dict: cluster → {PCs, labels, explained_var}
    """
    results = {}

    for c in glial_clusters:
        df_all = []
        labels = []

        for mouse, df in mice.items():
            df_c = df[df["cluster_leiden"] == c]
            if len(df_c) == 0:
                continue
            df_all.append(df_c[gene_cols])
            labels += [mouse] * len(df_c)

        if not df_all:
            continue

        X = pd.concat(df_all)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        pca = PCA(n_components=3)
        PCs = pca.fit_transform(X_scaled)

        results[c] = {
            "PCs": PCs,
            "labels": labels,
            "explained_var": pca.explained_variance_ratio_,
            "scaler": scaler,
            "pca": pca,
        }

    return results


def compute_de_per_cluster(mice_z, gene_cols, glial_clusters, disease_map):
    """
    MWU-based differential expression (z-scored data).
    Returns dict: cluster → DataFrame of genes with p & logFC.
    """
    de_results = {}

    for c in glial_clusters:
        tg_cells = []
        wt_cells = []

        for mouse, df in mice_z.items():
            df_c = df[df["cluster_leiden"] == c][gene_cols]
            if disease_map[mouse] == "TG":
                tg_cells.append(df_c)
            else:
                wt_cells.append(df_c)

        if not tg_cells or not wt_cells:
            continue

        tg_mat = pd.concat(tg_cells)
        wt_mat = pd.concat(wt_cells)

        rows = []
        for g in gene_cols:
            try:
                _, p = mannwhitneyu(tg_mat[g], wt_mat[g], alternative="two-sided")
            except:
                p = 1.0

            logFC = tg_mat[g].mean() - wt_mat[g].mean()

            rows.append({"gene": g, "p": p, "logFC": logFC})

        de_results[c] = pd.DataFrame(rows).sort_values("logFC", ascending=False)

    return de_results


def extract_top_genes(de_results, n=20):
    """
    For each cluster: top N upregulated and top N downregulated genes.
    """
    top_genes_per_cluster = {}

    for c, df in de_results.items():
        top_up = df.nlargest(n, "logFC")
        top_down = df.nsmallest(n, "logFC")
        top_genes_per_cluster[c] = {"upregulated": top_up, "downregulated": top_down}

    return top_genes_per_cluster


def build_top_genes_heatmap_matrix(de_results, top_n=10):
    """
    Build matrix for heatmap: cluster × gene.
    """
    all_top = []

    for c, df_de in de_results.items():
        df = df_de.copy()
        df["cluster"] = c
        df = df.nlargest(top_n, "logFC")
        all_top.append(df)

    matrix_df = pd.concat(all_top)
    return matrix_df.pivot_table(
        index="cluster", columns="gene", values="logFC", fill_value=0
    )


def compute_ad_specificity(mice_z, gene_cols, glial_clusters, disease_map):
    """
    AD-specificity score = (TG_mean - WT_mean) - WT_mean
    High = strongly induced only in TG.
    """
    ad_results = {}

    for c in glial_clusters:
        tg_cells = []
        wt_cells = []

        for mouse, df in mice_z.items():
            df_c = df[df["cluster_leiden"] == c][gene_cols]
            if disease_map[mouse] == "TG":
                tg_cells.append(df_c)
            else:
                wt_cells.append(df_c)

        if not tg_cells or not wt_cells:
            continue

        TG = pd.concat(tg_cells).mean()
        WT = pd.concat(wt_cells).mean()

        df_out = pd.DataFrame(
            {"gene": gene_cols, "TG_mean": TG.values, "WT_mean": WT.values}
        )

        df_out["logFC"] = df_out["TG_mean"] - df_out["WT_mean"]
        df_out["disease_specific_score"] = df_out["logFC"] - df_out["WT_mean"]

        ad_results[c] = df_out.sort_values("disease_specific_score", ascending=False)

    return ad_results


def build_ad_specific_heatmap(ad_specific, top_n=10):
    """
    Build matrix for AD-specificity heatmap: cluster × gene
    """
    rows = []
    for c, df in ad_specific.items():
        for _, row in df.head(top_n).iterrows():
            rows.append(
                {
                    "cluster": c,
                    "gene": row["gene"],
                    "score": row["disease_specific_score"],
                }
            )

    heat_df = pd.DataFrame(rows)
    return heat_df.pivot_table(
        index="cluster", columns="gene", values="score", fill_value=0
    )


def rank_clusters_by_ad_specificity(ad_specific):
    rows = []
    for c, df in ad_specific.items():
        mean_score = df["disease_specific_score"].abs().mean()
        top_gene = df.loc[df["disease_specific_score"].idxmax(), "gene"]
        top_val = df["disease_specific_score"].max()
        rows.append(
            {
                "cluster": c,
                "mean_AD_specificity": mean_score,
                "top_gene": top_gene,
                "top_gene_score": top_val,
            }
        )
    return pd.DataFrame(rows).sort_values("mean_AD_specificity", ascending=False)
