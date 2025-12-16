import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const RQ5Section = () => {
  return (
    <section id="rq-5" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky left panel */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
              <div className="text-xs uppercase tracking-wide text-muted-foreground">
                Research Question 5
              </div>
              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                When modeling plaque distance, which features are most important?
              </h2>

              <div className="mt-5 space-y-3 text-sm">
                <a
                  href="#rq5-modalities"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Feature modalities</div>
                  <div className="text-muted-foreground">Genes, morphology, spatial, cluster</div>
                </a>

                <a
                  href="#rq5-linear"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Linear models</div>
                  <div className="text-muted-foreground">OLS vs PLS performance</div>
                </a>

                <a
                  href="#rq5-ablation"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Modality ablation</div>
                  <div className="text-muted-foreground">Which modality drives R²?</div>
                </a>

                <a
                  href="#rq5-spatial-bias"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Spatial-only bias</div>
                  <div className="text-muted-foreground">Why R² can be inflated</div>
                </a>

                <a
                  href="#rq5-across-mice"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Across-mice diagnostics</div>
                  <div className="text-muted-foreground">Distribution + JSD</div>
                </a>

                <a
                  href="#rq5-zscores"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Gene signature strategy</div>
                  <div className="text-muted-foreground">Within-cluster z-normalization</div>
                </a>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: click to jump between subsections while the RQ panel stays visible.
              </div>
            </div>
          </aside>

          {/* Main analysis */}
          <div className="space-y-12">
            {/* =========================
                Intro + modalities
               ========================= */}
            <div id="rq5-modalities" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">RQ5 Analysis</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To study feature importance (and improve linear performance), we augment the 347-gene
                expression vector with morphology, spatial coordinates, and a cell-type indicator
                derived from Leiden clustering.
              </p>

              <div className="grid gap-4 md:grid-cols-2">
                <div className="rounded-2xl border border-border bg-card p-5">
                  <div className="text-sm font-semibold text-foreground">Modalities</div>
                  <ul className="mt-3 space-y-2 text-sm text-muted-foreground">
                    <li>• <span className="text-foreground font-medium">Genes</span>: 347 expression values</li>
                    <li>• <span className="text-foreground font-medium">Morphology</span>: cell area + nucleus area</li>
                    <li>• <span className="text-foreground font-medium">Spatial</span>: (x, y) centroid coordinates</li>
                    <li>• <span className="text-foreground font-medium">Cluster</span>: Leiden cluster (cell type)</li>
                  </ul>
                </div>

                <div className="rounded-2xl border border-border bg-card p-5">
                  <div className="text-sm font-semibold text-foreground">Goal</div>
                  <p className="mt-3 text-sm text-muted-foreground leading-relaxed">
                    Compare how much predictive signal comes from each modality, and whether high
                    spatial-only performance reflects true plaque proximity or anatomical bias.
                  </p>
                </div>
              </div>
            </div>

            {/* =========================
                Linear models
               ========================= */}
            <div id="rq5-linear" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Linear models</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Adding non-gene features improves linear prediction. Ordinary linear regression reaches
                test R² ≈ 0.24, while PLS reaches ≈ 0.19. The small train–test gap suggests stable,
                generalizable linear structure.
              </p>

              {/* Replace filename with your interactive export if different */}
              <PlotFrame
                src={`${base}plots/linear_pls_performance.html`}
                title="Linear vs PLS performance"
                size="lg"
                caption="Interactive: train/test R² for Linear vs PLS models."
              />
            </div>

            {/* =========================
                Ablation
               ========================= */}
            <div id="rq5-ablation" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Modality ablation</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Ablation analysis shows strong variation across modalities and model classes.
                Linear models benefit mainly from gene expression, while nonlinear models extract
                much richer structure. Spatial-only models can appear very strong, but this can be
                misleading.
              </p>

              <div className="grid gap-8 lg:grid-cols-2">
                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Modality ablation (heatmap)
                  </h5>
                  <PlotFrame
                    src={`${base}plots/modality_ablation.html`}
                    title="Modality ablation"
                    size="md"
                    caption="Interactive: model performance across modality combinations."
                  />
                </div>

                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Best model per modality
                  </h5>
                  <PlotFrame
                    src={`${base}plots/best_model_per_modality.html`}
                    title="Best model per modality"
                    size="md"
                    caption="Interactive: best-performing model for each modality subset."
                  />
                </div>
              </div>
            </div>

            {/* =========================
                Spatial-only bias
               ========================= */}
            <div id="rq5-spatial-bias" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Why spatial-only can be “too good”</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Spatial-only models can exploit region-dependent plaque density (anatomical vulnerability)
                rather than learning geometric proximity. High R² may reflect anatomical confounding rather
                than a pathology-sensitive plaque-distance predictor.
              </p>
            </div>

            {/* =========================
                Across mice diagnostics
               ========================= */}
            <div id="rq5-across-mice" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Across-mice diagnostics</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Comparing spatial-model predictions across Tg and WT animals reveals nearly identical
                prediction distributions across mice. This supports the conclusion that the model is
                learning conserved anatomical structure rather than disease progression.
              </p>

              <div className="grid gap-8 lg:grid-cols-2">
                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Predicted plaque distance (across mice)
                  </h5>
                  <PlotFrame
                    src={`${base}plots/predicted_plaque_distance.html`}
                    title="Predicted plaque distance comparison"
                    size="md"
                    caption="Interactive: compare predicted-score distributions across mice."
                  />
                </div>

                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Distribution across mice
                  </h5>
                  <PlotFrame
                    src={`${base}plots/distribution_across_mice.html`}
                    title="Distribution across mice"
                    size="md"
                    caption="Interactive: prediction distributions per mouse (shape + overlap)."
                  />
                </div>
              </div>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <h5 className="text-lg font-semibold text-foreground mb-3">
                  Jensen–Shannon divergence
                </h5>
                <PlotFrame
                  src={`${base}plots/JS_Divergence.html`}
                  title="Jensen–Shannon divergence"
                  size="lg"
                  caption="Interactive: pairwise divergence between mouse prediction distributions."
                />
              </div>
            </div>

            {/* =========================
                Z-score framework
               ========================= */}
            <div id="rq5-zscores" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">A gene-signature alternative</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Because cross-mouse normalization is unreliable for sparse Xenium panels, we adopt a
                within-cluster and within-gene z-normalized signature approach. Averaging z-scores across
                the 16 PIGs yields a plaque-induced activation score per cluster, enabling robust
                comparisons across disease status and age progression.
              </p>

              <div className="grid gap-8 lg:grid-cols-2">
                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Mean PIG per mouse
                  </h5>
                  <PlotFrame
                    src={`${base}plots/mean_PIG_per_mouse.html`}
                    title="Mean PIG per mouse"
                    size="md"
                    caption="Interactive: cluster-level PIG activation across mice."
                  />
                </div>

                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Disease effect & age progression
                  </h5>
                  <PlotFrame
                    src={`${base}plots/disease_effect_age_progression.html`}
                    title="Disease effect and age progression"
                    size="md"
                    caption="Interactive: disease-linked activation and aging trends by cluster."
                  />
                </div>
              </div>

              <div className="grid gap-8 lg:grid-cols-2">
                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    Age progression
                  </h5>
                  <PlotFrame
                    src={`${base}plots/age_progression.html`}
                    title="Age progression"
                    size="md"
                    caption="Interactive: age trajectory of plaque-induced activation."
                  />
                </div>

                <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                  <h5 className="text-lg font-semibold text-foreground mb-3">
                    PCA clusters
                  </h5>
                  <PlotFrame
                    src={`${base}plots/PCA_clusters.html`}
                    title="PCA clusters"
                    size="md"
                    caption="Interactive: PCA projections showing WT vs Tg separation within clusters."
                  />
                </div>
              </div>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <h5 className="text-lg font-semibold text-foreground mb-3">
                  AD-specific genes
                </h5>
                <PlotFrame
                  src={`${base}plots/ad_specific_genes.html`}
                  title="AD-specific genes"
                  size="lg"
                  caption="Interactive: heatmap of disease-specific gene activation per cluster."
                />
              </div>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <h5 className="text-lg font-semibold text-foreground mb-3">
                  Top genes per glial cluster
                </h5>
                <PlotFrame
                  src={`${base}plots/top_genes_per_glial.html`}
                  title="Top genes per glial cluster"
                  size="lg"
                  caption="Interactive: differential expression highlights in glial clusters."
                />
              </div>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <h5 className="text-lg font-semibold text-foreground mb-3">
                  Age progression WT vs Tg
                </h5>
                <PlotFrame
                  src={`${base}plots/age_progression_wt_tg.html`}
                  title="Age progression WT vs Tg"
                  size="lg"
                  caption="Interactive: divergence of Tg vs WT trajectories across age."
                />
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ5Section;
