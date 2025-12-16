import PlotFrame from "@/components/PlotFrame";
import { Brain } from 'lucide-react'

const base = import.meta.env.BASE_URL;

const RQ3Section = () => {
  return (
    <section id="rq-3" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
            {/* Sticky panel */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Brain className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">Research Question 3</div>
                  <div className="text-xs text-muted-foreground">
                    How does PIG expression change in plaque proximity?
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#rq3-mean"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Mean PIG expression vs distance
                </a>
                <a
                  href="#rq3-regression"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Regression analysis
                </a>
                <a
                  href="#rq3-halves"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Half-distance summary
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>
          {/*
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
              <div className="text-xs uppercase tracking-wide text-muted-foreground">
                Research Question 3
              </div>
              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                How does PIG expression change in plaque proximity?
              </h2>

              <div className="mt-5 space-y-3 text-sm">
                <a
                  href="#rq3-mean"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">
                    Mean PIG expression vs distance
                  </div>
                  <div className="text-muted-foreground">
                    Binned means + 95% CI across 5 distance bins
                  </div>
                </a>

                <a
                  href="#rq3-regression"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">
                    Regression analysis
                  </div>
                  <div className="text-muted-foreground">
                    Slopes + distances to halve expression
                  </div>
                </a>

                <a
                  href="#rq3-halves"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">
                    Half-distance summary
                  </div>
                  <div className="text-muted-foreground">
                    d₁/₂ per gene (interactive)
                  </div>
                </a>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: use these links, or scroll the analysis on the right.
              </div>
            </div>
          </aside> Sticky left panel */}

          {/* Main analysis */}
          <div className="space-y-12">
            {/* =========================
                SECTION 1 — Mean PIG expression
               ========================= */}
            <div id="rq3-mean" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                RQ3 Analysis
              </h3>

              <h4 className="text-xl font-semibold text-foreground">
                Mean PIG expression at different plaque distances
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                In this section, we investigate how the expression of the 16
                plaque-induced genes changes with distance to the nearest
                plaque. We group cells into 5 equal-count distance bins and
                compute the mean log1p-normalized transcript count within each
                bin along with the 95% confidence interval. ANOVA confirms that
                all 16 PIGs show significant differences in mean expression
                across distance bins at Bonferroni-corrected FDR set to 0.01.
              </p>

              <PlotFrame
                src={`${base}plots/PIG_expression_vs_distance.html`}
                title="Expression of PIGs vs. distance to plaque"
                size="xl"
                caption="Interactive: mean log1p expression per distance bin with 95% CI."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                From the figure above, <span className="font-medium text-foreground">Gfap</span>{" "}
                shows the biggest difference in log1p expression between the closest
                and furthest distance bins (0.72). This is roughly two orders of
                magnitude larger than the SEM (~0.01) thanks to the large number of
                cells per bin (~10,779). On the natural scale, this corresponds to{" "}
                <span className="font-medium text-foreground">exp(0.72) ≈ 2.05×</span>{" "}
                more expression in the closest bin compared to the furthest bin, with
                a monotonic decrease from near to far.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                The visualization also highlights consistent decreasing gradients for
                microglial markers (e.g., Hexb, Ctsd, Cst3, Apoe) and astrocytic markers
                (e.g., Gfap, Serpina3n, Vim). Regions most proximal to plaques show
                elevated microglial and astrocytic expression that fades with distance.
              </p>
            </div>

            {/* =========================
                SECTION 2 — Regression analysis
               ========================= */}
            <div id="rq3-regression" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">
                Regression analysis: PIG expression vs. plaque distance
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                For each of the 16 PIGs, we regress log1p-normalized transcript counts
                against distance to the nearest plaque. We apply Benjamini–Hochberg
                False Discovery Rate (FDR) adjustment across the 16 independent tests.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                All 16 PIGs exhibit a statistically significant negative slope at the
                0.01 FDR level, but slopes vary considerably. The smallest and largest
                absolute slopes correspond to Cxcl10 and Gfap, respectively. Translating
                slopes into an intuitive distance scale, we compute the distance required
                to halve expression (d₁/₂).
              </p>
            </div>

            {/* =========================
                SECTION 3 — Half-distances
               ========================= */}
            <div id="rq3-halves" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">
                Distances to halve expression (d₁/₂)
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We convert each gene’s slope into a half-distance{" "}
                <span className="font-medium text-foreground">d₁/₂ = ln(2) / |β|</span>.{" "}
                For example, Gfap halves over ~128 µm, while Cxcl10 halves over ~4,415 µm,
                indicating nearly constant expression at the scale of the mouse brain.
                Cxcl10 is also extremely zero-inflated, which reduces apparent spatial variation.
              </p>

              <PlotFrame
                src={`${base}plots/distances_to_halve_expression.html`}
                title="Distances to halve expression for the 16 PIGs"
                size="lg"
                caption="Interactive: d₁/₂ per PIG computed from regression slopes."
              />
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ3Section;
