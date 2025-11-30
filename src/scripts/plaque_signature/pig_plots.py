import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import ptitprince as pt


# -------------------------------------------------------------------
# 1. BASIC CLUSTER PLOTS
# -------------------------------------------------------------------

def plot_pig_by_cluster(pig_glial_df):
    """Barplot of mean PIG scores per glial cluster per mouse."""
    plt.figure(figsize=(14,6))
    sns.barplot(data=pig_glial_df, x="cluster", y="mean_pig", hue="mouse")
    plt.title("Z-normalized PIG Scores (Glial Clusters Only)")
    plt.ylabel("Mean PIG z-score")
    plt.show()


# -------------------------------------------------------------------
# 2. VOLCANO PLOT: Disease Effect vs Age Progression
# -------------------------------------------------------------------

def plot_volcano(summary):
    df = summary.reset_index()
    plt.figure(figsize=(10,8))
    sns.scatterplot(
        data=df, x="TG_minus_WT", y="abs_slope",
        hue="cluster", size="r",
        sizes=(50,200), palette="tab20"
    )
    plt.axvline(0, color='grey', ls='--')
    plt.axhline(0, color='grey', ls='--')
    plt.xlabel("Disease Effect (TG − WT)")
    plt.ylabel("|Age Slope| (Magnitude)")
    plt.title("Cluster Volcano: Disease Effect vs Age Progression")
    plt.tight_layout()
    plt.show()


# -------------------------------------------------------------------
# 3. DISEASE MAP: WT baseline vs TG activation
# -------------------------------------------------------------------

def plot_disease_map(summary):
    df = summary.reset_index()
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=df,
        x="WT_mean",
        y="TG_minus_WT",
        hue="cluster",
        palette="tab20",
        s=120
    )

    plt.axhline(0, color="grey", linestyle="--")
    plt.axvline(0, color="grey", linestyle="--")
    plt.xlabel("Baseline WT Expression (PIG Score)")
    plt.ylabel("Disease Activation (TG − WT)")
    plt.title("Disease Map: Baseline vs Activation")
    plt.show()


# -------------------------------------------------------------------
# 4. TG-ONLY AGE PROGRESSION (all clusters on one plot)
# -------------------------------------------------------------------

def plot_age_progression(pig_glial_df, glial_clusters):
    plt.figure(figsize=(12,7))

    for c in glial_clusters:
        df = pig_glial_df[(pig_glial_df["cluster"] == c) & (pig_glial_df["disease"]=="TG")]
        df = df.sort_values("age")
        plt.plot(df["age"], df["mean_pig"], marker="o", label=f"Cluster {c}")

    plt.xlabel("Age (months)")
    plt.ylabel("Mean PIG (TG only)")
    plt.title("Age Progression of PIG Scores — TG Only")
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
    plt.show()


# -------------------------------------------------------------------
# 5. VIOLIN PLOT (TG vs WT, sorted)
# -------------------------------------------------------------------

def plot_violin(pig_glial_df, order):
    plt.figure(figsize=(12,6))
    sns.violinplot(
        data=pig_glial_df,
        x="cluster",
        y="mean_pig",
        hue="disease",
        order=order,
        split=True,
        inner="quartile",
        linewidth=0.8,
        palette={"TG": "#d62728", "WT": "#1f77b4"}
    )
    plt.title("TG vs WT PIG Score per Glial Cluster (Sorted)")
    plt.xlabel("Cluster")
    plt.ylabel("PIG Gene Score (z-normalized)")
    plt.show()


# -------------------------------------------------------------------
# 6. RAINCLOUD PLOT
# -------------------------------------------------------------------

def plot_raincloud(pig_glial_df, order):
    plt.figure(figsize=(12,6))
    pt.RainCloud(
        data=pig_glial_df,
        x="cluster",
        y="mean_pig",
        hue="disease",
        order=order,
        palette={"TG":"#d62728","WT":"#1f77b4"},
        alpha=.8,
        width_viol=.6,
        width_box=.3,
        move=0.2
    )
    plt.title("TG vs WT PIG Score per Cluster (Raincloud)")
    plt.show()


# -------------------------------------------------------------------
# 7. FACETGRID VIOLINS
# -------------------------------------------------------------------

def plot_facet_violin(pig_glial_df, order):
    g = sns.FacetGrid(
        pig_glial_df,
        col="cluster",
        col_wrap=5,
        sharey=False,
        col_order=order,
        hue="disease",
        palette={"TG":"#d62728","WT":"#1f77b4"}
    )
    g.map(sns.violinplot, "disease", "mean_pig", split=False, inner="box")
    g.add_legend()
    plt.show()


# -------------------------------------------------------------------
# 8. PCA PLOT PER CLUSTER
# -------------------------------------------------------------------

def plot_cluster_pcas(pca_results, disease_map, glial_clusters):
    num_clusters = len(glial_clusters)
    cols = 4
    rows = int(np.ceil(num_clusters / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(cols*4, rows*4))
    axes = axes.flatten()

    for i, c in enumerate(glial_clusters):
        ax = axes[i]
        res = pca_results[c]

        df_plot = pd.DataFrame({
            "PC1": res["PCs"][:,0],
            "PC2": res["PCs"][:,1],
            "mouse": res["labels"]
        })
        df_plot["disease"] = df_plot["mouse"].map(disease_map)

        sns.scatterplot(
            data=df_plot,
            x="PC1", y="PC2",
            hue="disease",
            alpha=0.6, s=10, ax=ax,
            palette={"TG":"#d62728","WT":"#1f77b4"}
        )

        ax.set_title(f"Cluster {c}")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

        if i != 0:
            ax.get_legend().remove()

    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle("PCA Disease Axis per Glial Cluster", fontsize=18, y=1.02)
    plt.tight_layout()
    plt.show()


# -------------------------------------------------------------------
# 9. TOP-10 DE HEATMAP
# -------------------------------------------------------------------

def plot_top_de_heatmap(matrix):
    """
    Plot heatmap of top differential genes.
    Expects a pivoted matrix: cluster × gene.
    """
    plt.figure(figsize=(20, 8))
    sns.heatmap(matrix, cmap="coolwarm", center=0, linewidths=0.5)
    plt.title("Top Differential Genes per Glial Cluster (logFC)", fontsize=18)
    plt.xlabel("Gene")
    plt.ylabel("Cluster")
    plt.show()

# -------------------------------------------------------------------
# 10. TG vs WT AGE CURVES (one panel per cluster)
# -------------------------------------------------------------------

def plot_age_curves_by_cluster(pig_glial_df, glial_clusters):
    rows = int(np.ceil(len(glial_clusters)/4))
    cols = 4
    fig, axes = plt.subplots(rows, cols, figsize=(22, 5*rows))
    axes = axes.flatten()

    for i, c in enumerate(glial_clusters):
        ax = axes[i]
        df_c = pig_glial_df[pig_glial_df["cluster"] == c]

        # TG
        df_tg = df_c[df_c["disease"]=="TG"].groupby("age")["mean_pig"].mean().reset_index()
        if len(df_tg):
            ax.plot(df_tg["age"], df_tg["mean_pig"],
                    marker="o", linewidth=2.5, color="#d62728", label="TG")

        # WT
        df_wt = df_c[df_c["disease"]=="WT"].groupby("age")["mean_pig"].mean().reset_index()
        if len(df_wt):
            ax.plot(df_wt["age"], df_wt["mean_pig"],
                    marker="o", linewidth=2.5, color="#1f77b4", label="WT")

        ax.set_title(f"Cluster {c}", fontsize=16)
        ax.set_xlabel("Age")
        ax.set_ylabel("Mean PIG (z-score)")

    for j in range(i+1, rows*cols):
        fig.delaxes(axes[j])

    plt.suptitle("Age Progression — TG vs WT", fontsize=26)
    plt.tight_layout()
    plt.show()


# -------------------------------------------------------------------
# 11. AD-SPECIFIC GENE HEATMAP
# -------------------------------------------------------------------

def plot_ad_specific_heatmap(matrix):
    """
    Plot heatmap of AD-specific gene scores.
    Expects a pivoted matrix: cluster × gene.
    """
    plt.figure(figsize=(20, 8))
    sns.heatmap(matrix, cmap="viridis", linewidths=0.5)
    plt.title("Most AD-Specific Genes (High TG, Low WT)", fontsize=20)
    plt.xlabel("Gene")
    plt.ylabel("Cluster")
    plt.show()



# -------------------------------------------------------------------
# 12. AD-SPECIFIC SCATTERPLOT (per cluster)
# -------------------------------------------------------------------

def plot_ad_scatter_per_cluster(disease_specific_genes):
    for c, df in disease_specific_genes.items():
        plt.figure(figsize=(10,12))
        sns.scatterplot(
            data=df,
            x="WT_mean",
            y="logFC",
            size="disease_specific_score",
            hue="disease_specific_score",
            palette="viridis",
            sizes=(20,200)
        )
        plt.title(f"AD-Specific Gene Activation — Cluster {c}")
        plt.xlabel("WT Expression")
        plt.ylabel("TG − WT logFC")
        plt.axhline(0, color="black", ls="--")
        plt.axvline(0, color="black", ls="--")
        plt.tight_layout()
        plt.show()
