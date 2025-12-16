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
            Explore clustering of cells based on the 347-dimensional gene expression vector.
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
      </div>
    </section>
  );
};

export default CellClusteringSection;
