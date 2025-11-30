import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt

def run_pca_per_cluster(mice, gene_cols, disease_map, glial_clusters):
    results = {}
    for c in glial_clusters:
        df_all = []
        labels = []
        for mouse, df in mice.items():
            df_c = df[df["cluster_leiden"] == c]
            if len(df_c)==0:
                continue
            df_all.append(df_c[gene_cols])
            labels += [mouse]*len(df_c)
        if not df_all:
            continue
        X = pd.concat(df_all)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        pca = PCA(n_components=3)
        PCs = pca.fit_transform(X_scaled)
        results[c] = {
            "PCs": PCs, "labels": labels,
            "explained_var": pca.explained_variance_ratio_
        }
    return results


def compute_de_per_cluster(mice_z, gene_cols, disease_map, glial_clusters):
    de_results = {}
    for c in glial_clusters:
        tg = []
        wt = []
        for mouse, df in mice_z.items():
            df_c = df[df["cluster_leiden"] == c]
            (tg if disease_map[mouse]=="TG" else wt).append(df_c[gene_cols])
        if not tg or not wt:
            continue
        tg_m = pd.concat(tg)
        wt_m = pd.concat(wt)
        rows=[]
        for g in gene_cols:
            try:
                _, p = mannwhitneyu(tg_m[g], wt_m[g])
            except:
                p = 1
            logFC = tg_m[g].mean() - wt_m[g].mean()
            rows.append({"gene":g,"p":p,"logFC":logFC})
        de_results[c] = pd.DataFrame(rows).sort_values("logFC",ascending=False)
    return de_results


def compute_ad_specificity(mice_z, gene_cols, disease_map, glial_clusters):
    out = {}
    for c in glial_clusters:
        tg = []
        wt = []
        for mouse,df in mice_z.items():
            df_c = df[df["cluster_leiden"] == c][gene_cols]
            (tg if disease_map[mouse]=="TG" else wt).append(df_c)
        if not tg or not wt:
            continue
        TG = pd.concat(tg).mean()
        WT = pd.concat(wt).mean()
        df = pd.DataFrame({"gene":gene_cols,
                           "TG_mean":TG.values,
                           "WT_mean":WT.values})
        df["logFC"]=df["TG_mean"]-df["WT_mean"]
        df["disease_specific_score"]=df["logFC"]-df["WT_mean"]
        out[c]=df.sort_values("disease_specific_score",ascending=False)
    return out
