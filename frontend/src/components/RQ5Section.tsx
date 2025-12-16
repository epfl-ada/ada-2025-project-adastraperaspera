import PlotFrame from "@/components/PlotFrame";
import { Brain } from 'lucide-react'
const base = import.meta.env.BASE_URL;

const linearPerf = [
  { model: "Linear", trainR2: 0.25, testR2: 0.24 },
  { model: "PLS", trainR2: 0.19, testR2: 0.19 },
];

const Th = ({ children }: { children: React.ReactNode }) => (
  <th className="px-4 py-3 text-left text-xs font-semibold text-foreground">
    {children}
  </th>
);

const Td = ({ children }: { children: React.ReactNode }) => (
  <td className="px-4 py-3 text-sm text-muted-foreground align-top">
    {children}
  </td>
);

const RQ5Section = () => {
  return (
    <section id="rq-5" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
            <aside className="lg:sticky lg:top-24">
                            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
                                {/* Header (same pattern as GeneExpression / Microscopy) */}
                                <div className="flex items-start gap-3">
                                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                                    <Brain className="w-5 h-5 text-primary" />
                                </div>
                                <div>
                                    <div className="text-sm font-semibold text-foreground">
                                    Research Question 5
                                    </div>
                                    <div className="text-xs text-muted-foreground">
                                    Which features matter most for predicting plaque distance?
                                    </div>
                                </div>
                                </div>
            
                                {/* Navigation (identical hover / spacing / typography) */}
                                <nav className="mt-6 space-y-2 text-sm">
                                <a
                                    href="#rq5-linear"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                    Linear models
                                </a>
            
                                <a
                                    href="#rq5-ablation"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                    Modality ablation
                                </a>
            
                                <a
                                    href="#rq5-spatial-diagnostics"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                    Spatial-only diagnostics
                                </a>
            
                                <a
                                    href="#rq5-zscore-framework"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                    Z-normalized signatures
                                </a>
                                <a
                                    href="#rq5-disease-age"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                   Disease + age progression
                                </a>
                                <a
                                    href="#rq5-ad-specific"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                   AD-specific genes
                                </a>
                                </nav>
            
                                {/* Footer hint */}
                                <div className="mt-6 text-xs text-muted-foreground">
                                Tip: scroll or use the navigation above.
                                </div>
                            </div>
                        </aside>
          {/* 
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
              <div className="text-xs uppercase tracking-wide text-muted-foreground">
                Research Question 5
              </div>
              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                Which features matter most for predicting plaque distance?
              </h2>

              <div className="mt-5 space-y-3 text-sm">
                <a
                  href="#rq5-linear"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Linear models</div>
                  <div className="text-muted-foreground">Modalities + performance table</div>
                </a>

                <a
                  href="#rq5-ablation"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Modality ablation</div>
                  <div className="text-muted-foreground">Ablation heatmap + best model</div>
                </a>

                <a
                  href="#rq5-spatial-diagnostics"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Spatial-only diagnostics</div>
                  <div className="text-muted-foreground">Cross-mouse comparisons</div>
                </a>

                <a
                  href="#rq5-zscore-framework"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Z-normalized signatures</div>
                  <div className="text-muted-foreground">Mean PIG per cluster per mouse</div>
                </a>

                <a
                  href="#rq5-disease-age"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Disease + age progression</div>
                  <div className="text-muted-foreground">Effect sizes, trajectories, PCA</div>
                </a>

                <a
                  href="#rq5-ad-specific"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">AD-specific genes</div>
                  <div className="text-muted-foreground">Heatmaps + top DE + age curves</div>
                </a>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Figures below are expected in{" "}
                <span className="font-mono">frontend/public/plots</span>.
              </div>
            </div>
          </aside>Sticky left panel */}

          {/* Main content */}
          <div className="space-y-12">
            <div className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                RQ5: When modeling plaque distance, which features are most important?
              </h3>
            </div>
             {/* =========================
                Intro + Modalities
                ========================= */}
            <div id="rq5-modalities" className="space-y-4">
                <h3 className="text-2xl font-bold text-foreground">RQ5 Analysis</h3>

                <p className="text-lg text-muted-foreground leading-relaxed">
                    To improve linear performance and study feature importance, we augment the original
                    347-gene expression vector with morphology, spatial coordinates, and a cell-type
                    indicator derived from Leiden clustering.
                </p>

                <div className="grid gap-4 md:grid-cols-2">
                    <div className="rounded-2xl border border-border bg-card p-5">
                    <div className="text-sm font-semibold text-foreground">Modalities</div>
                    <ul className="mt-3 space-y-2 text-sm text-muted-foreground">
                        <li>• <span className="font-medium text-foreground">Genes</span>: 347 expression values</li>
                        <li>• <span className="font-medium text-foreground">Morphology</span>: cell area + nucleus area</li>
                        <li>• <span className="font-medium text-foreground">Spatial</span>: (x, y) centroid coordinates</li>
                        <li>• <span className="font-medium text-foreground">Cluster</span>: Leiden cluster (cell type)</li>
                    </ul>
                    </div>

                    <div className="rounded-2xl border border-border bg-card p-5">
                    <div className="text-sm font-semibold text-foreground">Goal</div>
                    <p className="mt-3 text-sm text-muted-foreground leading-relaxed">
                        Compare how much predictive signal comes from each modality and determine whether
                        high spatial-only performance reflects true plaque proximity or anatomical bias.
                    </p>
                    </div>
                </div>
            </div>


            {/* =========================
                Linear models: modalities + perf table
               ========================= */}
            <div id="rq5-linear" className="space-y-6">

              <p className="text-lg text-muted-foreground leading-relaxed">
                We also evaluated partial least squares (PLS) regression in addition to an ordinary
                linear model. Performance:
              </p>

              <div className="rounded-2xl border border-border bg-card overflow-hidden">
                <div className="px-6 py-4 border-b border-border">
                  <div className="text-sm font-semibold text-foreground">
                    Linear-model performance (reported)
                  </div>
                </div>
                <div className="overflow-x-auto">
                  <table className="min-w-full">
                    <thead className="bg-muted/40">
                      <tr>
                        <Th>Model</Th>
                        <Th>Train R²</Th>
                        <Th>Test R²</Th>
                      </tr>
                    </thead>
                    <tbody>
                      {linearPerf.map((r) => (
                        <tr key={r.model} className="border-t border-border">
                          <Td>{r.model}</Td>
                          <Td>{r.trainR2.toFixed(2)}</Td>
                          <Td>{r.testR2.toFixed(2)}</Td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>

              <p className="text-lg text-muted-foreground leading-relaxed">
                The linear model with extra features closes part of the gap to boosted trees and
                generalizes well (small train–test gap). PLS is slightly lower, consistent with
                learning only a dominant global axis correlated with plaque distance.
              </p>
            </div>

            {/* =========================
                Modality ablation
               ========================= */}
            <div id="rq5-ablation" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Modality ablation</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We next evaluate which modalities drive performance across model classes using
                ablation experiments.
              </p>

              <PlotFrame
                src={`${base}plots/modality_ablation.html`}
                title="Modality ablation"
                size="lg"
                caption="Interactive ablation heatmap (genes / morphology / spatial / cluster)."
              />

              <PlotFrame
                src={`${base}plots/best_model_per_modality.html`}
                title="Best model per modality"
                size="lg"
                caption="Interactive summary of best-performing model for each modality configuration."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Ablation reveals that performance varies strongly by modality and model class. Linear
                models reach modest accuracy when genes are included, and fall near chance when only
                morphology or spatial coordinates are used. Nonlinear models extract richer structure.
                Notably, spatial-only models can reach very high R² in the 17.9-month Tg mouse, which
                motivates deeper diagnostics.
              </p>
            </div>

            {/* =========================
                Spatial-only diagnostics
               ========================= */}
            <div id="rq5-spatial-diagnostics" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Diagnostics: why spatial-only looks strong</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Spatial-only performance can be inflated by anatomical bias: plaques are not uniformly
                distributed across the brain, so coordinates can proxy regional vulnerability rather
                than geometric proximity. We test this by comparing spatial predictions across mice.
              </p>

              <PlotFrame
                src={`${base}plots/predicted_plaque_distance.html`}
                title="Spatial predictions comparison"
                size="lg"
                caption="Interactive comparison of predicted distance distributions across mice."
              />

              <PlotFrame
                src={`${base}plots/distribution_across_mice.html`}
                title="Distribution across mice"
                size="lg"
                caption="Interactive overlay of predicted-score distributions (Tg vs WT, different ages)."
              />

              <PlotFrame
                src={`${base}plots/JS_Divergence.html`}
                title="Jensen–Shannon divergence"
                size="lg"
                caption="Interactive divergence matrix comparing predicted-score distributions across mice."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                These results show predicted-score distributions are extremely similar across mice,
                implying spatial-only models learn conserved anatomy rather than pathology-sensitive
                plaque proximity.
              </p>
            </div>

            {/* =========================
                Switch to z-normalized signatures
               ========================= */}
            <div id="rq5-zscore-framework" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">
                Shift to within-cluster, within-gene z-normalized signatures
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Because cross-mouse transcript normalization is unreliable (batch-like distortions and
                sparse panel constraints), we adopt within-cluster and within-gene z-normalized
                signatures. We average z-scores across the 16 PIGs to obtain a plaque-induced gene
                activation score per cluster.
              </p>

              <PlotFrame
                src={`${base}plots/mean_PIG_per_mouse.html`}
                title="Mean PIG per mouse"
                size="lg"
                caption='Interactive: cluster-level "plaque-induced gene activation score" across mice.'
              />
            </div>

            {/* =========================
                Disease status + age progression
               ========================= */}
            <div id="rq5-disease-age" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">
                Disease specificity and age progression
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                With the z-score framework in place, we analyze disease status (TG vs WT) and age
                progression (2 → 5 → 17 months). We compute disease-specific activation scores per
                cluster and rank clusters by mean disease specificity.
              </p>

              <PlotFrame
                src={`${base}plots/disease_effect_age_progression.html`}
                title="Disease effect and age progression"
                size="lg"
                caption="Interactive summary of TG vs WT effects and progression patterns."
              />

              <PlotFrame
                src={`${base}plots/age_progression.html`}
                title="Age progression"
                size="lg"
                caption="Interactive age progression view across clusters."
              />

              <PlotFrame
                src={`${base}plots/PCA_clusters.html`}
                title="PCA clusters"
                size="lg"
                caption="Interactive PCA-based projections: WT vs TG separation within clusters."
              />
            </div>

            {/* =========================
                AD-specific genes + supporting heatmaps
               ========================= */}
            <div id="rq5-ad-specific" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">AD-specific activation in glial clusters</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Cluster 8 (microglia) shows the strongest AD-specific activation, followed by cluster
                18 (astrocytes). Disease-specific genes (e.g., Syngr1, Gfap, Sparcl1) exhibit strong
                induction in Tg while remaining mostly silent in WT, consistent with glial reactivity
                and inflammatory remodeling.
              </p>

              <PlotFrame
                src={`${base}plots/ad_specific_genes.html`}
                title="AD-specific genes"
                size="lg"
                caption="Interactive heatmap of AD-specific genes by cluster."
              />

              <PlotFrame
                src={`${base}plots/top_genes_per_glial.html`}
                title="Top genes per glial cluster"
                size="lg"
                caption="Interactive heatmap of top differential genes in glial clusters."
              />

              <PlotFrame
                src={`${base}plots/age_progression_wt_tg.html`}
                title="Age progression WT vs Tg"
                size="lg"
                caption="Interactive trajectories: Tg increases with age in key glial clusters; WT remains stable."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Taken together, these analyses show robust and progressive transcriptional activation
                in biologically relevant glial populations driven by amyloid pathology, while
                spatial-only modeling primarily captures conserved anatomy.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ5Section;
