import { Brain, Target,  Activity , CheckCircle2 } from "lucide-react";

const RQ_LINK_CLASS =
  "block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition";

const IntroductionSection = () => {
  return (
    <section id="introduction" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-7xl">
        {/* Section title */}
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Introduction
          </h2>
        </div>

        {/* 3 cards always on one line */}
        <div className="flex flex-nowrap gap-8 items-start">
          {/* Card 1 — Motivation */}
          <div className="flex-1 p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4">
              <div className="w-12 h-12 rounded-xl bg-chart-3/10 flex items-center justify-center flex-shrink-0">
                <Activity  className="w-6 h-6 text-chart-3" />
              </div>

              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">
                  Motivation
                </h3>

                <p className="text-muted-foreground leading-relaxed mb-4">
                  Amyloid beta plaques (Aβ) are a hallmark of Alzheimer&apos;s disease (AD)
                  with known effects including changes in gene expression, glial activation,
                  and neuronal death. However, most studies examine plaque influence using
                  coarse distance bins.
                </p>

                <p className="text-muted-foreground leading-relaxed">
                  Building a finer model of the plaque-induced microenvironment enables
                  actionable downstream insights. Knowing <em>which</em> cells and <em>which</em>{" "}
                  genes respond <em>where</em> around plaques helps drug developers pre-filter
                  targets accessible from the vasculature.
                </p>
              </div>
            </div>
          </div>

          {/* Card 2 — Project goals */}
          <div className="flex-1 p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4">
              <div className="w-12 h-12 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                <Target className="w-6 h-6 text-primary" />
              </div>

              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">
                  Project goals
                </h3>

                <p className="text-muted-foreground leading-relaxed mb-4">
                  We use a Xenium murine spatial transcriptomics dataset to develop and interpret
                  quantitative models of plaque proximity.
                </p>

                <ul className="space-y-3 text-muted-foreground">
                  <li className="flex items-start gap-2">
                    <CheckCircle2 className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>Study temporal and spatial trends in gene expression</span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CheckCircle2 className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>Quantify time- and distance-dependent changes in cell composition</span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CheckCircle2 className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>Benchmark predictive models to infer plaque distance from multiple modalities</span>
                  </li>
                </ul>
              </div>
            </div>
          </div>

          {/* Card 3 — Research questions */}
          <div className="flex-1 p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4">
              <div className="w-12 h-12 rounded-xl bg-accent/10 flex items-center justify-center flex-shrink-0">
                <Brain className="w-6 h-6 text-accent" />
              </div>

              <div className="w-full">
                <h3 className="text-xl font-semibold text-foreground mb-3">
                  Research questions
                </h3>

                <p className="text-muted-foreground leading-relaxed mb-4">
                  In the upcoming sections, we address the following RQs:
                </p>

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

                <p className="text-muted-foreground leading-relaxed mt-4">
                  We report statistical significance and include diagnostics to critically evaluate models.
                </p>
              </div>
            </div>
          </div>
        </div>

        {/* Optional: for horizontal scrolling on small screens instead of overflow
            Replace the flex container above with:
            <div className="flex flex-nowrap gap-8 overflow-x-auto pb-2">
            and add min widths to cards:
            className="min-w-[340px] flex-1 ..."
        */}
      </div>
    </section>
  );
};

export default IntroductionSection;
