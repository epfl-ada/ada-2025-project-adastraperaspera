import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const GeneExpressionChart = () => {
  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">

        {/* Section title (fixed, not scrolling) */}
        <div className="mb-8">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Gene expression
          </h2>
          <p className="text-lg text-muted-foreground max-w-3xl">
            Spatial transcriptomics overview and microscopy-based cell maps.
          </p>
        </div>

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

          {/* PANEL — Text + microscopy plot */}
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
            {/* Text */}
            <p className="text-lg text-muted-foreground">
              We are dealing with a spatial transcriptomics dataset which contains
              single cell gene expression measurements of 347 genes. The gene
              expression matrix tends to be sparse. For the transgenic mouse at
              17.9 months of age, 302 out of 347 genes have zero median transcript count.
            </p>

            <p className="text-lg text-muted-foreground">
              Out of 347 genes, 248 represent canonical cell-type markers, 83 genes
              relate to activated microglia and astrocytes, and 16 are plaque-induced
              genes curated from literature.
            </p>

            {/* Microscopy plot BELOW text */}
            <PlotFrame
              src={`${base}plots/microscopy_cells_unified.html`}
              title="Microscopy cells unified"
              size="lg"
            />

            <p className="text-xs text-muted-foreground text-center">
              Interactive microscopy-based spatial maps. Pan and zoom enabled.
            </p>
          </div>

        </div>

        {/* Scroll hint */}
        <p className="mt-4 text-xs text-muted-foreground">
          Tip: scroll horizontally to explore the analysis →
        </p>

      </div>
    </section>
  );
};

export default GeneExpressionChart;
