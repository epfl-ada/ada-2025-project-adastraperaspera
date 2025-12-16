import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const CellClusteringSection = () => {
  return (
    <section className="py-24 bg-background">
      <div className="container mx-auto px-6">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Cell Type Analysis
          </h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            We will now explore clustering of cells based on the 347-dimensional gene expression vector.
            We apply Leiden clustering with 15 nearest neighbors on PCA-reduced gene expression space.
            As a result, we obtained $K=19$ clusters. The clusters are then reduced to 2 dimensions using UMAP; the resulting plot is shown in the following figure.
          </p>
          <p className="text-lg text-muted-foreground max-w-3xl mx-auto mb-12 text-center">
          Superimposing the color-coded clusters onto the brain tissue, we can see that the
          gene expression-based clustering strongly correlates with the brain morphology.
        </p>
        </div>
        <div className="grid lg:grid-cols-2 gap-12 items-start">
          <div className="bg-card rounded-2xl border border-border p-8 shadow-md">
            <h3 className="text-lg font-semibold text-foreground mb-4">
              Joint clustering UMAP
            </h3>
            <PlotFrame
              src={`${base}plots/joint_clustering_umap.html`}
              title="Joint clustering UMAP"
              size="md"
              caption="UMAP embedding colored by Leiden clusters (interactive)."
            />
          </div>

          <div className="bg-card rounded-2xl border border-border p-8 shadow-md">
            <h3 className="text-lg font-semibold text-foreground mb-4">
              Spatial overlay of clusters
            </h3>
            <PlotFrame
              src={`${base}plots/joint_clustering_overlayed.html`}
              title="Spatial overlay of clusters"
              size="md"
              caption="Spatial coordinates colored by Leiden clusters (interactive)."
            />
          </div>
        </div>
        <p className="text-lg text-muted-foreground max-w-3xl mx-auto mb-12 text-center">
          Cluster 6 (i.e., vascular cells) is almost exclusively localized to the outer rim of the brain, corresponding to epidural space. Hippocampal formation shows a distinct cluster 14 which follows the elongated shape of the dentate gyrus. Looking at the inferred cell type, we can confirm that cluster 14 corresponds to Dentate gyrus immature neurons (glutamatergic). Meanwhile, the region hosting amygdala and hypothalamus is dominated by cluster 10. Cluster 10 corresponds to Hypothalamic medial mammillary glutamatergic neurons, matching its observed localization. The ventricle cavities are lined with distinct cluster of cells (ID 15). This cluster corresponds to Hypothalamic GnRH1-expressing glutamatergic neurons. This finding is expected since hypothalamus forms the floor and part of the lateral walls of the third ventricle.
        </p>
      </div>
    </section>
  );
};

export default CellClusteringSection;
