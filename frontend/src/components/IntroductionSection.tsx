import { Brain } from "lucide-react";

const RQ_LINK_CLASS =
  "block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition";

const IntroductionSection = () => {
  return (
    <section id="introduction" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-7xl">
        {/* Section title */}
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground">
            Introduction
          </h2>
        </div>

        {/* Two-column layout */}
        <div className="grid gap-10 lg:grid-cols-[1fr_360px] items-start">
          {/* LEFT — Overview text */}
          <div className="space-y-6">
            <h3 className="text-2xl font-semibold text-foreground">
              Overview
            </h3>

            <p className="text-lg text-muted-foreground leading-relaxed">
              Amyloid beta plaques (Aβ) are a known hallmark of Alzheimer's disease (AD)
              with known effects including changes in gene expression, glial activation
              and neuronal death. However, the bulk of the existing research only examines
              this influence across rough distance bins. Building a finer model of
              plaque-induced microenvironment has important downstream applications.
              With the knowledge of <em>which</em> cells and <em>which</em> genes respond
              <em>where</em> around the plaques, drug developers can pre-filter
              therapeutic targets accessible from the vasculature.
            </p>

            <p className="text-lg text-muted-foreground leading-relaxed">
              Using a Xenium dataset spanning{" "}
              <span className="font-medium text-foreground">351,714 cells</span>,{" "}
              <span className="font-medium text-foreground">347 genes</span>, and{" "}
              <span className="font-medium text-foreground">1,736 plaques</span>, we
              analyze how plaque proximity relates to cell-type composition, gene
              expression, neighborhood context, predictive modeling, and age- and
              genotype-specific signatures.
            </p>
          </div>

          {/* RIGHT — Research questions card */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3 mb-4">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Brain className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Research questions
                  </div>
                  <div className="text-xs text-muted-foreground">
                    Guiding the analysis
                  </div>
                </div>
              </div>

              <nav className="space-y-2 text-sm">
                <a href="#rq-1" className={RQ_LINK_CLASS}>
                  <span className="font-semibold text-foreground">RQ1:</span>{" "}
                  How does the cell type composition change in plaque proximity?
                </a>

                <a href="#rq-2" className={RQ_LINK_CLASS}>
                  <span className="font-semibold text-foreground">RQ2:</span>{" "}
                  How are the cell type composition, PIG expression, and plaque distance related?
                </a>

                <a href="#rq-3" className={RQ_LINK_CLASS}>
                  <span className="font-semibold text-foreground">RQ3:</span>{" "}
                  How does the Plaque Induced Gene (PIG) expression change in plaque proximity?
                </a>

                <a href="#rq-4" className={RQ_LINK_CLASS}>
                  <span className="font-semibold text-foreground">RQ4:</span>{" "}
                  When modeling plaque distance, which feature modalities are most important?
                </a>

                <a href="#rq-5" className={RQ_LINK_CLASS}>
                  <span className="font-semibold text-foreground">RQ5:</span>{" "}
                  How does the gene expression change with age for each cell type and mouse group?
                </a>
              </nav>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default IntroductionSection;
