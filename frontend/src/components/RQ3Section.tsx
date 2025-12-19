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

      {/* keep your content, only fix structure */}
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
                SECTION 1 — Mean PIG expression
               ========================= */}
            <div id="rq3-mean" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                RQ3 : How does PIG expression change in plaque proximity?
              </h3>
              <h4 className="text-xl font-semibold text-foreground">
                Distance-binned means, confidence intervals, and ANOVA
              </h4>

              {/* keep your content, only fix structure */}
              <ul className="space-y-3 text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    We group cells into{" "}
                    <span className="font-medium text-foreground">
                      5 equal-count distance bins
                    </span>{" "}
                    and compute, for each PIG:
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  mean{" "}
                  <span className="font-medium text-foreground">log1p-normalized</span>{" "}
                  transcript count per bin
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  95% confidence interval (SEM-based)
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                An ANOVA confirms that{" "}
                <span className="font-medium text-foreground">all 16 PIGs</span> have
                significant differences in mean expression across distance bins at{" "}
                <span className="font-medium text-foreground">
                  Bonferroni-corrected FDR = 0.01
                </span>
                .
              </p>

              <PlotFrame
                src={`${base}plots/PIG_expression_vs_distance.html`}
                title="Expression of PIGs vs. distance to plaque"
                size="xl"
                caption="Distance-binned mean expression of all 16 PIGs with 95% C, showing consistent plaque-proximal elevation for glial/immune markers.."
              />

              {/* keep your content, only fix structure */}
              <ul className="space-y-3 text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>
                      A prominent example is{" "}
                      <span className="font-medium text-foreground">Gfap</span>, which shows
                      the largest proximal-to-distal mean difference on the log1p scale:
                    </span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Δ(log1p mean) ≈{" "}
                    <span className="font-medium text-foreground">0.72</span> between
                    closest and farthest bins
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    SEM ≈ <span className="font-medium text-foreground">0.01</span> (with{" "}
                    <span className="font-medium text-foreground">10,779 cells per bin</span>)
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    On the natural scale: exp(0.72) ≈{" "}
                    <span className="font-medium text-foreground">2.05×</span> higher
                    expression near plaques
                  </span>
                </li>

                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>Across genes, the most consistent gradients include:</span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Microglial markers (e.g.,{" "}
                    <span className="font-medium text-foreground">Hexb</span>,{" "}
                    <span className="font-medium text-foreground">Ctsd</span>,{" "}
                    <span className="font-medium text-foreground">Cst3</span>,{" "}
                    <span className="font-medium text-foreground">Apoe</span>)
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Astrocytic markers (e.g.,{" "}
                    <span className="font-medium text-foreground">Gfap</span>,{" "}
                    <span className="font-medium text-foreground">Serpina3n</span>,{" "}
                    <span className="font-medium text-foreground">Vim</span>)
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Collectively, tissue within ~0–29 µm of plaques shows elevated glial/immune
                signatures that fade with distance.
              </p>
            </div>

            {/* =========================
                SECTION 2 — Regression analysis
               ========================= */}
            <div id="rq3-regression" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">
                Per-gene regression slopes and “distance-to-half-expression”
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To summarize gradients continuously, we regress log1p-normalized transcript
                count against distance for each PIG and apply Benjamini–Hochberg FDR
                correction across the 16 regressions.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                All 16 PIGs exhibit a statistically significant negative slope at the 0.01
                FDR level, but slopes vary considerably. The smallest and largest absolute
                slopes correspond to Cxcl10 and Gfap, respectively. Translating slopes into
                an intuitive distance scale, we compute the distance required to halve
                expression (d₁/₂).
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Given the mouse brain diameter (~6,000 µm), Cxcl10’s gradient is effectively
                flat, consistent with its extreme zero inflation (low absolute expression
                variation across distance).
              </p>

              <PlotFrame
                src={`${base}plots/distances_to_halve_expression.html`}
                title="Distances to halve expression for the 16 PIGs"
                size="md"
                caption="Interactive: d₁/₂ per PIG computed from regression slopes."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                To increase explanatory power (R²), reduce heteroscedasticity, and test
                whether plaque shape contributes to local responses, we compute geometric
                properties of the nearest plaque for each cell:{" "}
                <span className="font-medium text-foreground">
                  area, perimeter, major axis length and orientation
                </span>
                .
              </p>

              <PlotFrame
                src={`${base}plots/geometric_characteristics_of_plaque.html`}
                title="Geometric characteristics of the nearest plaque"
                size="md"
                caption="Distributions of nearest-plaque geometry features (area, perimeter, major axis, orientation) summarized via box plots."
              />

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

              <p className="text-lg text-muted-foreground leading-relaxed">
                To validate that R is informative, we check the Count distribution and CDF.
                The distribution shows meaningful variation:{" "}
                <span className="font-medium text-foreground">55.2%</span> of cells have{" "}
                <span className="font-medium text-foreground">0</span> plaques within 61 µm,
                while the remainder have up to{" "}
                <span className="font-medium text-foreground">21</span> plaques within that
                radius—indicating a usable local density signal.
              </p>

              <PlotFrame
                src={`${base}plots/multi_plaque_proximity.html`}
                title="Multi-plaque proximity features"
                size="lg"
                caption="Histogram and CDF of local plaque density (Count within R = 61 µm), validating that multi-plaque proximity provides informative variation beyond nearest-plaque distance."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                To capture local cell–cell context and spatial signaling, we compute
                neighborhood mean expression features: for each cell and each PIG, we
                summarize the expression of the{" "}
                <span className="font-medium text-foreground">15 other PIGs</span> across its{" "}
                <span className="font-medium text-foreground">k = 100 nearest neighbors</span>.
              </p>

              {/* keep your content, only fix structure */}
              <ul className="space-y-3 text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    Motivation for k = 100:
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  balances locality with stability (not too small, not too large)
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  approximates neighborhood effects at a scale relevant to cell–cell
                  communication (~100 µm radius depending on density, assuming ~10 µm cell
                  diameter and relatively dense cell packing)
                </li>
              </ul>

              <PlotFrame
                src={`${base}plots/pigs_coexpression.html`}
                title="Multi-plaque proximity features"
                size="lg"
                caption="Histogram and CDF of local plaque density (Count within R = 61 µm), validating that multi-plaque proximity provides informative variation beyond nearest-plaque distance."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The matrix is not symmetric because target/neighbor roles are not
                commutative. The strongest relationships are:
              </p>

              <PigNeighborCorrTable />

              {/* keep your content, only fix structure */}
              <ul className="mt-3 space-y-3 text-sm text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    We systematically quantify how each feature group improves PIG prediction
                    using nested linear models for each PIG:
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Model 0:</span> distance only
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Model 1:</span> distance +
                    plaque geometry
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Model 2:</span> distance +
                    plaque geometry + multi-plaque proximity
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      Models 3₁, 3₂, 3₄, 3₈, 3₁₅:
                    </span>{" "}
                    Model 2 + neighborhood PIG context features using the top 1 / 2 / 4 / 8 /
                    15 neighbor PIGs (ranked by correlation)
                  </span>
                </li>
              </ul>
              <PlotFrame
                src={`${base}plots/pig_trajectories.html`}
                title="Multi-plaque proximity features"
                size="md"
                caption="Mean adjusted R² across PIGs for each nested model, including min/max ranges, showing which feature groups add meaningful predictive value."
              />

              {/* keep your content, only fix structure */}
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
                    How does PIG expression change in plaque proximity?
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
