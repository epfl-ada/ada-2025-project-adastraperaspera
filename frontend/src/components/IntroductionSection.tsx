import { Target, Brain, FlaskConical, TrendingUp } from "lucide-react";

const IntroductionSection = () => {
  return (
    <section id="introduction" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">

        {/* Section title */}
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Introduction
          </h2>
          <p className="text-lg text-muted-foreground max-w-3xl mx-auto">
            Scientific motivation, application context, and analytical strategy
            underlying this study.
          </p>
        </div>

        {/* =========================
            ROW 1 — Side-by-side cards
           ========================= */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-12">

          {/* Card 1 — Research Objective */}
          <div className="p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4 mb-6">
              <div className="w-12 h-12 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                <Target className="w-6 h-6 text-primary" />
              </div>
              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">
                  Research Objective
                </h3>
                <p className="text-muted-foreground leading-relaxed">
                  This project investigates how Aβ plaques impact the surrounding
                  microenvironment. We aim to develop a quantitative model to
                  precisely describe the influence of plaques on nearby tissue.
                  Identifying which genes and cells are affected at different
                  plaque distances refines our understanding of Alzheimer’s
                  disease progression.
                </p>
              </div>
            </div>
          </div>

          {/* Card 2 — Drug Development */}
          <div className="p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4 mb-6">
              <div className="w-12 h-12 rounded-xl bg-accent/10 flex items-center justify-center flex-shrink-0">
                <FlaskConical className="w-6 h-6 text-accent" />
              </div>
              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">
                  Drug Development Application
                </h3>
                <p className="text-muted-foreground leading-relaxed">
                  Developers of new therapeutics can leverage our model to
                  identify biologically realistic targets within plaque regions
                  that remain accessible from the vasculature.
                </p>
              </div>
            </div>
          </div>
        </div>

        {/* =========================
            ROW 2 — Full-width card
           ========================= */}
        <div className="p-8 rounded-2xl bg-card border border-border">
          <div className="flex items-start gap-4 mb-6">
            <div className="w-12 h-12 rounded-xl bg-chart-3/10 flex items-center justify-center flex-shrink-0">
              <Brain className="w-6 h-6 text-chart-3" />
            </div>
            <div>
              <h3 className="text-xl font-semibold text-foreground mb-3">
                Analysis Approach
              </h3>
              <p className="text-muted-foreground leading-relaxed mb-4">
                Our analysis explores how plaque proximity reshapes the cellular,
                molecular, and tissue environments. We quantify cell-to-plaque
                distances and apply rigorous statistical methods to characterize
                spatial trends in gene expression and cell composition.
              </p>

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <TrendingUp className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Analyze the interplay between gene expression and cell composition
                  </span>
                </li>
                <li className="flex items-start gap-2">
                  <TrendingUp className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Determine which phenomena drive observed spatial patterns
                  </span>
                </li>
                <li className="flex items-start gap-2">
                  <TrendingUp className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Benchmark predictive models to infer plaque distance from multigene expression
                  </span>
                </li>
              </ul>
            </div>
          </div>
        </div>

      </div>
    </section>
  );
};

export default IntroductionSection;
