import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const GeneExpressionChart = () => {
  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">

        {/* =========================
            TITLE + INTRO + SMALL PLOT
           ========================= */}
        <div className="mb-8 grid grid-cols-1 lg:grid-cols-[1fr_360px] gap-8 items-start">
          {/* Text */}
          <div>
            <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
              Gene expression
            </h2>
            <p className="text-lg text-muted-foreground max-w-3xl">
               We are dealing with a spatial transcriptomics dataset which contains single cell gene expression measurements of 347 genes. The gene expression matrix tends to be sparse. For the transgenic mouse at 17.9 months of age, 302 out of 347 genes have zero median transcript count.
              Looking at the gene selection, out of 347 genes, 248 represent markers for 8 main cell types, canonical neuronal cortical layer markers, and non-neuronal markers; 83 genes related to activated microglia and astrocytes; and 16 PIGs curated from primary literature.
            </p>
          </div>

          {/* Figure (no card) */}
          <figure className="hidden lg:block float-right ml-10 mb-6 w-[520px]">
            <img
              src={`${base}figures/microscopy_cells_unified.png`}
              alt="Microscopy cells unified"
              className="w-full rounded-2xl border border-border/60 shadow-md"
            />
            <figcaption className="mt-2 text-sm text-muted-foreground">
              Microscopy cells unified
            </figcaption>
          </figure>
</div>


        {/* Scroll hint */}
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
          {/* Block 1 */}
          <div
            className="
              snap-start
              min-w-[700px]
              max-w-[700px]
              bg-card
              rounded-2xl
              border border-border
              p-8
              flex-shrink-0
              space-y-6
            "
          >
            <p className="text-lg text-muted-foreground">
              First, let us explore the data by visualizing the distribution of log1p-transformed transcript counts for the 16 PIGs against the average distribution of all 347 genes. We will focus on the mouse with the most advanced stage of the Alzheimer's disease (transgenic at 17.9 months of age).
            </p>

            <PlotFrame
              src={`${base}figures/expression_distribution.png`}
              title="Gene expression distribution (log1p-transformed)"
              size="md"
              fit="contain"
            />

            <p className="text-lg text-muted-foreground">
              We can observe that for 13 out of 16 PIGs, the distribution has a mode at zero. Above zero, the support gradually drops off at higher transcript counts; however, different PIGs have a different rate of the density decay.
            </p>
          </div>

          {/* Block 2 */}
          <div
            className="
              snap-start
              min-w-[500px]
              max-w-[500px]
              bg-card
              rounded-2xl
              border border-border
              p-8
              flex-shrink-0
              space-y-6
            "
          >
            <p className="text-lg text-muted-foreground leading-relaxed">
              Next, we diagnose the PIGs with the most unusual expression patterns. To that end, we combine several diagnostic metrics (among which, zero-inflation, dispersion, and shape) into a composite score.
            </p>

            <table className="min-w-full text-sm">
              <thead className="bg-muted/50">
                <tr>
                  <th className="px-4 py-3 text-left font-semibold text-foreground">
                    gene
                  </th>
                  <th className="px-4 py-3 text-right font-semibold text-foreground">
                    zero_frac
                  </th>
                  <th className="px-4 py-3 text-right font-semibold text-foreground">
                    weird_score
                  </th>
                </tr>
              </thead>
              <tbody>
                {[
                  { gene: "Cxcl10", zero_frac: 1.0, weird_score: 9.35 },
                  { gene: "Cd74", zero_frac: 0.98, weird_score: 8.13 },
                  { gene: "Serpina3n", zero_frac: 0.88, weird_score: 0.79 },
                  { gene: "C4b", zero_frac: 0.9, weird_score: 0.19 },
                  { gene: "Gfap", zero_frac: 0.69, weird_score: 0.05 },
                ].map((row) => (
                  <tr
                    key={row.gene}
                    className="border-t border-border hover:bg-muted/30 transition"
                  >
                    <td className="px-4 py-2 font-mono text-foreground">
                      {row.gene}
                    </td>
                    <td className="px-4 py-2 text-right text-muted-foreground">
                      {row.zero_frac.toFixed(2)}
                    </td>
                    <td className="px-4 py-2 text-right text-muted-foreground">
                      {row.weird_score.toFixed(2)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>

            <p className="text-lg text-muted-foreground leading-relaxed">
              We find that among the 5 most unusual PIGs, Cxcl10 and Cd74 clearly stand out. Both have the weirdness score exceeding 8; the next highest is Serpina3n with a score an order of magnitude lower at 0.79. Looking at the diagnostic statistics above, we can see that both Cxcl10 and Cd74 have a very high zero proportion exceeding 0.97%.
            </p>
          </div>
        </div>
      </div>
    </section>
  );
};

export default GeneExpressionChart;
