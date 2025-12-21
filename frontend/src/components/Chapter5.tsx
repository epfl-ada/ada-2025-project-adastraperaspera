import { Activity } from "lucide-react";

const Chapter5 = () => {
  return (
    <section id="ch-5" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + content right */}
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky panel (LEFT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Activity className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 5
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Changing neighborhoods
                  </div>
                </div>
              </div>

              <div className="mt-6 space-y-3 text-sm text-muted-foreground leading-relaxed">
                <p>Cell composition shifts with distance.</p>
                <p>The tissue near plaques forms a distinct ecosystem.</p>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll to follow the transition.
              </div>
            </div>
          </aside>

          {/* Content (RIGHT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 5: Approaching the plaque
              </h2>
            </div>

            <div className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Now we watch what happens as we walk toward a plaque.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We model the relationship between <span className="font-medium text-foreground">cluster membership</span> and <span className="font-medium text-foreground">distance to the nearest plaque</span> using logistic regression. After fitting, we find statistically significant coefficients for <span className="font-medium text-foreground">14 of 19 clusters</span>, using a <span className="font-medium text-foreground">Bonferroni-adjusted p-value threshold of 0.01</span>.
              </p>

              <ul className="space-y-3 text-muted-foreground">
                <p className="text-lg text-muted-foreground leading-relaxed">
                To make coefficients interpretable, we translate them into:
                </p>
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">p(0):</span>{" "}
                  baseline probability of observing a cluster at the plaque surface
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">p(100):</span>{" "}
                  probability of observing a cluster at{" "}
                  <span className="font-medium text-foreground">100 µm</span> from the plaque
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">p(100) − p(0):</span>{" "}
                  relative change over{" "}
                  <span className="font-medium text-foreground">100 µm</span> away from the plaque
                </li>
              </ul>


              <p className="text-lg text-muted-foreground leading-relaxed">
                At the plaque surface, the mix of cells changes. Logistic
                regression makes it quantitative.{" "}
                <span className="font-semibold text-foreground">
                  Immune cells, vascular cells, and astrocytes become more likely
                  near plaques
                </span>
                , while{" "}
                <span className="font-semibold text-foreground">
                  many neuronal types become less likely
                </span>
                .
              </p>

              <PlotFrame
                src={`${base}plots/slopes_types.html`}
                title="Plaque distance by cell type"
                size="md"
                caption="Regression slopes for significant clusters."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                This is the sick mouse’s first clear symptom. The tissue near
                plaques is not the same tissue. It is a different ecosystem.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Some changes are gradual, others are sudden. Certain clusters
                shift sharply only after specific distance thresholds.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default Chapter5;
