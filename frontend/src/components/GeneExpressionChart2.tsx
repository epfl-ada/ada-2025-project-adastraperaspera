import PlotFrame from "@/components/PlotFrame";
import { Dna } from "lucide-react";

const base = import.meta.env.BASE_URL;

const GeneExpressionSection2 = () => {
  return (
    <section id="gene-expression" className="py-24 bg-muted/30 scroll-mt-24">
      <div className="container mx-auto px-6 max-w-7xl">
        {/* =========================
            HEADER
           ========================= */}
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Gene Expression
          </h2>
          <p className="text-lg text-muted-foreground max-w-3xl mx-auto">
            Overview and diagnostics of gene expression patterns in the Xenium
            spatial transcriptomics dataset.
          </p>
        </div>

        {/* =========================
            STICKY LEFT + FLOW RIGHT
           ========================= */}
        <div className="grid grid-cols-1 lg:grid-cols-12 gap-10 items-start">
          {/* =========================
              LEFT STICKY PANEL
             ========================= */}
          <aside className="lg:col-span-4">
            <div className="lg:sticky lg:top-24 rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Dna className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Gene Expression
                  </div>
                  <div className="text-xs text-muted-foreground">
                    Dataset overview & diagnostics
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#gene-expression-intro"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Dataset overview
                </a>
                <a
                  href="#gene-expression-distribution"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Expression distributions
                </a>
                <a
                  href="#gene-expression-weirdness"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Weird genes diagnostics
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* =========================
              RIGHT ANALYSIS FLOW
             ========================= */}
          <div className="lg:col-span-8 space-y-12">
            {/* =========================
                INTRO SECTION
               ========================= */}
            <section
              id="gene-expression-intro"
              className="scroll-mt-24 bg-card rounded-2xl border border-border p-8"
            >
              <div className="max-w-none">
                <figure
                  className="
                    hidden lg:block
                    float-right
                    ml-6
                    mb-4
                    w-[520px]
                    max-w-[60%]
                  "
                >
                  <img
                    src={`${base}figures/microscopy_.png`}
                    alt="Microscopy cells unified"
                    className="w-full h-auto rounded-xl border border-border"
                  />
                  <figcaption className="mt-1 text-xs text-muted-foreground">
                    Microscopy cells unified
                  </figcaption>
                </figure>

                <p className="text-lg text-muted-foreground leading-relaxed">
                  We are dealing with a spatial transcriptomics dataset which contains single cell gene
                  expression measurements of 347 genes. The gene expression matrix tends to be sparse.
                  For the transgenic mouse at 17.9 months of age, 302 out of 347 genes have zero median
                  transcript count. Looking at the gene selection, out of 347 genes, 248 represent markers
                  for 8 main cell types, canonical neuronal cortical layer markers, and non-neuronal markers;
                  83 genes related to activated microglia and astrocytes; and 16 PIGs curated from primary
                  literature.
                </p>

                <div className="clear-both" />
              </div>
            </section>

            {/* =========================
                DISTRIBUTION SECTION
               ========================= */}
            <section
              id="gene-expression-distribution"
              className="scroll-mt-24 bg-card rounded-2xl border border-border p-8 space-y-6"
            >
              <h3 className="text-xl font-semibold text-foreground">
                Expression distributions
              </h3>

              <p className="text-lg text-muted-foreground">
                First, let us explore the data by visualizing the distribution of
                log1p-transformed transcript counts for the 16 PIGs against the
                average distribution of all 347 genes. We will focus on the mouse
                with the most advanced stage of the Alzheimer's disease
                (transgenic at 17.9 months of age).
              </p>

              <PlotFrame
                src={`${base}figures/expression_distribution.png`}
                title="Gene expression distribution (log1p-transformed)"
                size="md"
                fit="contain"
              />

              <p className="text-lg text-muted-foreground">
                We can observe that for 13 out of 16 PIGs, the distribution has a
                mode at zero. Above zero, the support gradually drops off at
                higher transcript counts; however, different PIGs have a
                different rate of the density decay.
              </p>
            </section>

            {/* =========================
                WEIRDNESS SECTION
               ========================= */}
            <section
              id="gene-expression-weirdness"
              className="scroll-mt-24 bg-card rounded-2xl border border-border p-8 space-y-6"
            >
              <h3 className="text-xl font-semibold text-foreground">
                Weird genes diagnostics
              </h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Next, we diagnose the PIGs with the most unusual expression
                patterns. To that end, we combine several diagnostic metrics
                (among which, zero-inflation, dispersion, and shape) into a
                composite score.
              </p>

              <div className="overflow-x-auto">
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
              </div>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We find that among the 5 most unusual PIGs, Cxcl10 and Cd74
                clearly stand out. Both have the weirdness score exceeding 8;
                the next highest is Serpina3n with a score an order of magnitude
                lower at 0.79.
              </p>
            </section>
          </div>
        </div>
      </div>
    </section>
  );
};

export default GeneExpressionSection2;
