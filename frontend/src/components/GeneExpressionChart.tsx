import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const GeneExpressionChart = () => {
  return (
    <section className="py-24 bg-muted/30">
      {/* Section title + content */}
        <div className="mb-6 max-w-4xl">

          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Gene expression
          </h2>

          <p className="text-lg text-muted-foreground">
            We are dealing with a spatial transcriptomics dataset which contains
            single cell gene expression measurements of 347 genes. The gene
            expression matrix tends to be sparse. For the transgenic mouse at
            17.9 months of age, 302 out of 347 genes have zero median transcript count.
            Out of 347 genes, 248 represent canonical cell-type markers, 83 genes
            relate to activated microglia and astrocytes, and 16 are plaque-induced
            genes curated from literature.
          </p>

          {/* Plot BELOW the text */}
          <div className="mt-8">
            <PlotFrame
              src={`${base}plots/microscopy_cells_unified.html`}
              title="Microscopy cells unified"
              size="sm"
            />
          </div>

        </div>

        {/* Scroll hint (TOP) */}
        <p className="mb-4 text-xs text-muted-foreground">
          Tip: scroll horizontally to explore the analysis →
        </p>

        {/* =========================
            HORIZONTAL SCROLL STRIP
           ========================= */}
        <div
          className="
            relative
            flex
            flex-row
            gap-12
            overflow-x-auto
            overflow-y-hidden
            snap-x snap-mandatory
            pb-6
          "
        >
          <div
            className="
              snap-start
              min-w-[900px]
              max-w-[900px]
              bg-card
              rounded-2xl
              border border-border
              p-8
              flex-shrink-0
              space-y-6
            "
          >
          
            <PlotFrame
              src={`${base}plots/expression_distribution.html`}
              title="Microscopy cells unified"
              size="md"
            />
            <p className="text-xs text-muted-foreground text-center">
              Interactive microscopy-based spatial maps. Pan and zoom enabled.
            </p>
            <p className="text-lg text-muted-foreground">
              First, let us explore the data by visualizing the distribution of log1p-transformed transcript counts for the 16 PIGs against the average distribution of all 347 genes. We will focus on the mouse with the most advanced stage of the Alzheimer's disease (transgenic at 17.9 months of age).
            </p>

          </div>

        </div>
      </div>
    </section>
  );
};

export default GeneExpressionChart;
