import PlotFrame from "@/components/PlotFrame";
import { Brain, CircleChevronRight } from "lucide-react";

const base = import.meta.env.BASE_URL;

const PIG_CORR_ROWS = [
  { target: "Gfap", neighbor: "C4b", rho: "0.47" },
  { target: "Gfap", neighbor: "S100a6", rho: "0.47" },
  { target: "Gfap", neighbor: "Vim", rho: "0.45" },
];

function PigNeighborCorrTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">
        Local PIG–PIG correlations
      </h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Target PIG</th>
              <th className="py-2 pr-4 text-left font-medium">Neighbor PIG</th>
              <th className="py-2 text-right font-medium">Pearson correlation</th>
            </tr>
          </thead>

          <tbody>
            {PIG_CORR_ROWS.map((r, i) => (
              <tr
                key={`${r.target}-${r.neighbor}-${i}`}
                className="border-b border-border/60 last:border-b-0"
              >
                <td className="py-2 pr-4 text-foreground">{r.target}</td>
                <td className="py-2 pr-4 text-foreground">{r.neighbor}</td>
                <td className="py-2 text-right tabular-nums">{r.rho}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>


      <ul className="mt-3 space-y-3 text-sm text-muted-foreground">
        <li>
          <p className="text-lg text-muted-foreground leading-relaxed">
            Interpretations:
          </p>
        </li>

        <li className="flex items-start gap-2">
          •
          <span>
            <span className="font-medium text-foreground">Gfap / C4b:</span>{" "}
            reactive gliosis around plaques co-occurs with complement cascade activation;
            C4a/C4b is expressed in astrocytes, making co-variation expected.
          </span>
        </li>

        <li className="flex items-start gap-2">
          •{" "}
          <span>
            <span className="font-medium text-foreground">Gfap / S100a6:</span>{" "}
            S100a6 is upregulated in astrocytes in AD and concentrates around Aβ plaques,
            often alongside Gfap-positive reactive astrocytes.
          </span>
        </li>

        <li className="flex items-start gap-2">
          •{" "}
          <span>
            <span className="font-medium text-foreground">Gfap / Vim:</span>{" "}
            both are intermediate filament proteins upregulated in reactive astrogliosis;
            coordinated induction is expected in plaque-adjacent astrocytes.
          </span>
        </li>
      </ul>
    </div>
  );
}

const RQ3Section = () => {
  return (
    <section id="rq-3" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">
          {/* Main analysis */}
          <div className="space-y-12">
            {/* =========================
                SECTION 1 - Mean PIG expression
               ========================= */}
            <div id="rq3-mean" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                RQ3 : How does the PIG expression change in plaque proximity?
              </h3>
              <h4 className="text-xl font-semibold text-foreground">
                Distance-binned means, confidence intervals, and ANOVA
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We quantify plaque-induced gene behavior directly by binning cells into <span className="font-medium text-foreground">five equal-count distance bins</span> and computing mean log1p expression with uncertainty estimates.
                ANOVA confirms that <span className="font-medium text-foreground">all 16 PIGs differ significantly across distance bins</span>.

              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                The qualitative pattern is consistent:
              </p>

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  PIG expression is highest near plaques
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  Expression generally decreases with distance
                </li>
              </ul>

              <PlotFrame
                src={`${base}plots/PIG_expression_vs_distance.html`}
                title="Expression of PIGs vs. distance to plaque"
                size="xl"
                caption="Distance-binned mean expression of all 16 PIGs with 95% C, showing consistent plaque-proximal elevation for glial/immune markers.."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The magnitude and spatial scale of decay vary by gene:
              </p>

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <em>Gfap</em> shows the strongest and most localized gradient
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <em>Cxcl10</em> appears nearly flat due to extreme zero inflation
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To summarize this continuously, we regress expression on distance for each PIG and convert slopes into a more interpretable metric: <span className="font-medium text-foreground">distance-to-halve-expression</span>. This highlights gene-specific scales, with Gfap showing a steep decline and Cxcl10 appearing almost flat throughout the brain, consistent with its extreme zero inflation.
              </p>

              <PlotFrame
                src={`${base}plots/distances_to_halve_expression.html`}
                title="Distances to halve expression for the 16 PIGs"
                size="md"
                caption="Interactive: d₁/₂ per PIG computed from regression slopes."
              />
            </div>

            {/* =========================
                SECTION 2 - Regression analysis
               ========================= */}
            <div id="rq3-regression" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">
                Per-gene regression slopes and “distance-to-halve-expression”
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Nearest-plaque distance captures proximity to <em>one</em> plaque, but local
                pathology may also depend on plaque <em>crowding</em>. We therefore introduce
                two additional{" "}
                <span className="font-medium text-foreground">
                  multi-plaque proximity features
                </span>{" "}
                computed within a radius{" "}
                <span className="font-medium text-foreground">R = 61 µm</span> (the median
                nearest-plaque distance): the{" "}
                <span className="font-medium text-foreground">count</span> of plaques within
                this radius, and the{" "}
                <span className="font-medium text-foreground">mean distance</span> to plaques
                within the same neighborhood.
              </p>

              <PlotFrame
                src={`${base}plots/geometric_characteristics_of_plaque.html`}
                title="Geometric characteristics of the nearest plaque"
                size="md"
                caption="Distributions of nearest-plaque geometry features (area, perimeter, major axis, orientation) summarized via box plots."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                To validate that R is informative, we check the Count distribution and CDF.
                The distribution shows meaningful variation:{" "}
                <span className="font-medium text-foreground">55.2%</span> of cells have{" "}
                <span className="font-medium text-foreground">0</span> plaques within 61 µm,
                while the remainder have up to{" "}
                <span className="font-medium text-foreground">21</span> plaques within that
                radius-indicating a usable local density signal.
              </p>

              <PlotFrame
                src={`${base}plots/multi_plaque_proximity.html`}
                title="Multi-plaque proximity features"
                size="lg"
                caption="Histogram and CDF of local plaque density (Count within R = 61 µm), validating that multi-plaque proximity provides informative variation beyond nearest-plaque distance."
              />


              <p className="text-lg text-muted-foreground leading-relaxed">
                Plaque-associated transcriptional responses are not purely cell-autonomous and may reflect local microenvironmental context. To quantify this, we computed neighborhood expression summaries by averaging the expression of the other <span className="font-medium text-foreground">15 PIGs</span> across each cell’s <span className="font-medium text-foreground">100 nearest neighbors</span>. This neighborhood size balances locality with statistical stability, avoiding excessive noise from very small neighborhoods and oversmoothing from very large ones. Given typical cell diameters (~10 µm) and tissue density, this choice corresponds to an interaction scale on the order of <span className="font-medium text-foreground">~200 µm</span>, which is relevant for cell–cell signaling and coordinated glial responses. These neighborhood features allow us to model collective plaque-associated activation rather than isolated single-cell effects.
              </p>

              <PlotFrame
                src={`${base}plots/pigs_coexpression.html`}
                title="Multi-plaque proximity features"
                size="lg"
                caption="Histogram and CDF of local plaque density (Count within R = 61 µm), validating that multi-plaque proximity provides informative variation beyond nearest-plaque distance."
              />


              <PigNeighborCorrTable />

              <p className="text-lg text-muted-foreground leading-relaxed">
                To quantify the contribution of different spatial and contextual features to PIG prediction, we evaluated a series of <span className="font-medium text-foreground">nested linear models</span> separately for each PIG. The baseline model included <span className="font-medium text-foreground">distance to the nearest plaque</span> only, followed by successive additions of <span className="font-medium text-foreground">plaque geometry and local multi-plaque proximity </span>features. We then augmented this model with neighborhood transcriptional context, adding the top <span className="font-medium text-foreground">1, 2, 4, 8, or 15</span> most correlated neighboring PIG features. Neighbor PIGs were ranked by their correlation strength with the target gene. This nested design allows direct assessment of the incremental explanatory value of each new feature group.
              </p>


              <PlotFrame
                src={`${base}plots/pig_trajectories.html`}
                title="Multi-plaque proximity features"
                size="md"
                caption="Mean adjusted R² across PIGs for each nested model, including min/max ranges, showing which feature groups add meaningful predictive value."
              />


              <p className="text-lg text-muted-foreground leading-relaxed">
                Key findings:
              </p>

              <ul className="mt-3 space-y-3 text-sm text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    The largest average adjusted R² gain comes from adding the{" "}
                    <span className="font-medium text-foreground">
                      single most correlated neighborhood PIG
                    </span>{" "}
                    (mean improvement{" "}
                    <span className="font-medium text-foreground">+0.074</span>). This is
                    consistent with plaque proximity acting as a common confounder that drives
                    coordinated PIG activation.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Additional neighbor PIG proxies provide diminishing but still significant
                    gains (e.g.,{" "}
                    <span className="font-medium text-foreground">+0.036</span> for the next
                    increment), and the nested analysis favors using up to{" "}
                    <span className="font-medium text-foreground">15</span> neighbor PIGs even
                    though most benefit comes from the first proxy.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Gfap</span> achieves the
                    highest adjusted R² in{" "}
                    <span className="font-medium text-foreground">7 of 8</span> model variants,
                    consistent with it being strongly distance-linked.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Cxcl10</span> is lowest in{" "}
                    <span className="font-medium text-foreground">5 of 8</span> model variants,
                    consistent with extreme zero inflation limiting explainable variance.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Plaque geometry (Model 1) and multi-plaque proximity (Model 2) yield
                    statistically significant (but modest) improvements in adjusted R² for{" "}
                    <span className="font-medium text-foreground">15 of 16 PIGs</span> (average
                    improvements ~{" "}
                    <span className="font-medium text-foreground">0.0021</span> and{" "}
                    <span className="font-medium text-foreground">0.0063</span>, respectively),
                    implying these spatial descriptors are biologically salient but secondary to
                    neighborhood transcriptional context.
                  </span>
                </li>
              </ul>


            </div>
          </div>

          {/* Sticky panel */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Brain className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Research Question 3
                  </div>
                  <div className="text-xs text-muted-foreground">
                    How does the Plaque Induced Gene (PIG) expression change in plaque proximity?
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#rq3-mean"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Means, Intervals & ANOVA
                </a>
                <a
                  href="#rq3-regression"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Regression analysis
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default RQ3Section;
