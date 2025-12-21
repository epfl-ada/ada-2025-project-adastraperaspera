import { FlaskConical, Wrench } from "lucide-react";

const AppendixSection = () => {
  return (
    <section id="appendix" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Title */}
        <div className="text-center mb-12">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground">
            Appendix
          </h2>
        </div>

        <div className="grid gap-8 lg:grid-cols-2 items-start">
          {/* 7.a Multiple testing correction */}
          <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
            <div className="flex items-start gap-3 mb-4">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                <FlaskConical className="w-5 h-5 text-primary" />
              </div>
              <div>
                <div className="text-sm font-semibold text-foreground">
                  7.a Multiple testing correction
                </div>
                <div className="text-xs text-muted-foreground">
                  Controlling false positives across many tests
                </div>
              </div>
            </div>

            <ul className="space-y-3 text-muted-foreground">
              <li className="leading-relaxed">
                We used Bonferroni correction in RQ1 for per-cluster logistic
                regressions and in RQ3 for ANOVA across distance bins.
              </li>
              <li className="leading-relaxed">
                Additionally, we applied Benjamini–Hochberg FDR correction in RQ3
                when building per-PIG distance regressions and nested model
                comparisons.
              </li>
            </ul>
          </div>

          {/* 7.b Feature engineering */}
          <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
            <div className="flex items-start gap-3 mb-4">
              <div className="w-10 h-10 rounded-xl bg-accent/10 flex items-center justify-center flex-shrink-0">
                <Wrench className="w-5 h-5 text-accent" />
              </div>
              <div>
                <div className="text-sm font-semibold text-foreground">
                  7.b Feature engineering
                </div>
                <div className="text-xs text-muted-foreground">
                  Definitions used in regression models
                </div>
              </div>
            </div>

            <ul className="space-y-3 text-muted-foreground">
              <li className="leading-relaxed">
                When discussing the features for linear regression, we use the
                following terms:
              </li>
              <li className="leading-relaxed">
                Nearest plaque geometry refers to plaque area, perimeter, major
                axis length, and orientation.
              </li>
              <li className="leading-relaxed">
                Multi-plaque proximity refers to the plaque count and the average
                distance within R. In our case, we set R = 61 µm, which is the
                median cell-to-plaque distance.
              </li>
              <li className="leading-relaxed">
                Neighborhood context refers to the average expression of other
                PIGs in 100 nearest neighbors.
              </li>
            </ul>
          </div>
        </div>
      </div>
    </section>
  );
};

export default AppendixSection;
