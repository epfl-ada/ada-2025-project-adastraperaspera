import pandas as pd


def identify_cluster_celltypes(df_example, gene_cols, marker_sets):
    """
    Assign a dominant cell type to each cluster using average marker expression.
    """
    cluster_celltype = {}
    for c in sorted(df_example["cluster_leiden"].unique()):
        df_c = df_example[df_example["cluster_leiden"] == c]
        scores = {
            t: df_c[genes].mean().mean() if len(genes) > 0 else 0
            for t, genes in marker_sets.items()
        }
        cluster_celltype[c] = max(scores, key=scores.get)
    return cluster_celltype


def get_glial_clusters(cluster_celltype, glial_types=("micro", "astro", "oligo")):
    return [c for c, t in cluster_celltype.items() if t in glial_types]
