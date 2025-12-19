import PlotFrame from "@/components/PlotFrame";
import { Brain, CircleChevronRight } from "lucide-react";
const base = import.meta.env.BASE_URL;

const MODEL_R2_ROWS = [
  { model: "Linear", train: "0.25", test: "0.24" },
  { model: "PLS", train: "0.19", test: "0.19" },
];

const MODEL_FULL_R2_ROWS = [
  { model: "LASSO", train: "0.1969", test: "0.1955" },
  { model: "ElasticNet", train: "0.1969", test: "0.1955" },
  { model: "Random Forest", train: "0.2252", test: "0.1874" },
  { model: "XGBoost", train: "0.4784", test: "0.2577" },
];

const RESIDUAL_CLUSTER_ROWS = [
  {
    cluster: "15",
    type: "Hypothalamic GnRH1-expressing glutamatergic neurons",
    n: "917",
    meanAbs: "32.2295",
  },
  {
    cluster: "18",
    type: "Pons glutamatergic neurons",
    n: "499",
    meanAbs: "17.8746",
  },
];

const FEATURE_RANK_ROWS = [
  {
    rank: "1",
    en: "Gfap (7.738978)",
    lasso: "Gfap (7.864686)",
    rf: "Gfap (0.228450)",
    xgb: "Spag16 (0.033058)",
  },
  {
    rank: "2",
    en: "Spag16 (3.780172)",
    lasso: "Spag16 (3.801531)",
    rf: "Spag16 (0.084622)",
    xgb: "Lyz2 (0.023342)",
  },
  {
    rank: "3",
    en: "Slc17a6 (3.025042)",
    lasso: "Slc17a6 (3.073279)",
    rf: "Lyz2 (0.066671)",
    xgb: "Gfap (0.016999)",
  },
  {
    rank: "4",
    en: "Slc17a7 (2.953878)",
    lasso: "Slc17a7 (3.057380)",
    rf: "Cabp7 (0.063918)",
    xgb: "Igf2 (0.015752)",
  },
  {
    rank: "5",
    en: "Igf2 (2.907111)",
    lasso: "Igf2 (2.949014)",
    rf: "B2m (0.034584)",
    xgb: "Strip2 (0.013022)",
  },
];

function FeatureRankTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">
        Top predictive genes by model
      </h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-right font-medium">Rank</th>
              <th className="py-2 pr-4 text-left font-medium">ElasticNet</th>
              <th className="py-2 pr-4 text-left font-medium">LASSO</th>
              <th className="py-2 pr-4 text-left font-medium">Random Forest</th>
              <th className="py-2 text-left font-medium">XGBoost</th>
            </tr>
          </thead>

          <tbody>
            {FEATURE_RANK_ROWS.map((r) => (
              <tr key={r.rank} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-right tabular-nums text-foreground">
                  {r.rank}
                </td>
                <td className="py-2 pr-4 text-foreground whitespace-nowrap">{r.en}</td>
                <td className="py-2 pr-4 text-foreground whitespace-nowrap">{r.lasso}</td>
                <td className="py-2 pr-4 text-foreground whitespace-nowrap">{r.rf}</td>
                <td className="py-2 text-foreground whitespace-nowrap">{r.xgb}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ResidualByClusterTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">
        Mean absolute residuals by cluster
      </h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-right font-medium">Cluster (Leiden)</th>
              <th className="py-2 pr-4 text-left font-medium">Inferred cell type</th>
              <th className="py-2 pr-4 text-right font-medium">n</th>
              <th className="py-2 text-right font-medium">Mean |residual| (µm)</th>
            </tr>
          </thead>

          <tbody>
            {RESIDUAL_CLUSTER_ROWS.map((r) => (
              <tr key={r.cluster} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-right tabular-nums text-foreground">
                  {r.cluster}
                </td>
                <td className="py-2 pr-4 text-muted-foreground min-w-[360px]">{r.type}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.n}</td>
                <td className="py-2 text-right tabular-nums">{r.meanAbs}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ModelPerformanceTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">Performance :</h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Model</th>
              <th className="py-2 pr-4 text-right font-medium">Train R²</th>
              <th className="py-2 text-right font-medium">Test R²</th>
            </tr>
          </thead>

          <tbody>
            {MODEL_FULL_R2_ROWS.map((r) => (
              <tr key={r.model} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground">{r.model}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.train}</td>
                <td className="py-2 text-right tabular-nums">{r.test}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ModelR2Table() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">
        Performance (R²):
      </h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Model</th>
              <th className="py-2 pr-4 text-right font-medium">Train R²</th>
              <th className="py-2 text-right font-medium">Test R²</th>
            </tr>
          </thead>

          <tbody>
            {MODEL_R2_ROWS.map((r) => (
              <tr key={r.model} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground">{r.model}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.train}</td>
                <td className="py-2 text-right tabular-nums">{r.test}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

const RQ4Section = () => {
  return (
    <section id="rq-4" className="py-24 bg-background">
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
                    Research Question 4
                  </div>
                  <div className="text-xs text-muted-foreground">
                    How accurately can we infer plaque distance from gene expression?
                  </div>
                </div>
              </div>

              {/* Navigation (identical hover / spacing / typography) */}
              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#rq4-performance"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Model performance
                </a>

                <a
                  href="#rq4-modalities"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Additional feature modalities
                </a>
              </nav>

              {/* Footer hint */}
              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Main content */}
          <div className="space-y-12">
            <div className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                RQ4: When modeling plaque distance, which features are most important?
              </h3>
            </div>

            {/* =========================
                Performance table + plot
               ========================= */}
            <div id="rq4-performance" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Model performance</h4>

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    We benchmark models that predict plaque distance from the{" "}
                    <span className="font-medium text-foreground">347-gene expression vector</span>,
                    using a random{" "}
                    <span className="font-medium text-foreground">80/20 train/test split over cells</span>.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Linear models: <span className="font-medium text-foreground">LASSO</span>,{" "}
                    <span className="font-medium text-foreground">ElasticNet</span> with standardized inputs
                    (log1p-transformed counts)
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Tree models: <span className="font-medium text-foreground">Random Forest</span>,{" "}
                    <span className="font-medium text-foreground">XGBoost</span> with unstandardized inputs
                    (log1p-transformed counts)
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We report train/test R² and inspect residual structure and feature importance.
              </p>

              <ModelPerformanceTable />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Residual diagnostics for XGBoost show clear heteroscedasticity: near plaques the model tends
                to over-predict, while far from plaques it increasingly under-predicts
              </p>

              <PlotFrame
                src={`${base}plots/residual_diagnostics.html`}
                title="Residual diagnostics (residual vs true distance)"
                size="lg"
                caption="Residual diagnostics for the XGBoost distance model, highlighting heteroscedasticity and systematic bias across the true-distance range."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The target distribution is also highly concentrated in the 0–100 µm range, meaning models
                are effectively trained on a narrow band of distances.
              </p>

              <PlotFrame
                src={`${base}plots/true_vs_predicted.html`}
                title="True vs predicted plaque distance"
                size="lg"
                caption="True vs predicted plaque distance and related distributional diagnostics, illustrating concentration of mass at short distances and systematic residual structure."
              />

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>We further categorize residuals:</span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">25.51%</span>: highly negative residuals
                    (<span className="font-medium text-foreground">−17.2 µm</span>), mostly close to plaques
                    (within <span className="font-medium text-foreground">49 µm</span>)
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">28.15%</span>: highly positive residuals
                    (&gt;<span className="font-medium text-foreground">10.6 µm</span>), mostly far from plaques
                    (beyond <span className="font-medium text-foreground">87 µm</span>)
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">20.94%</span>: moderate residuals
                    (<span className="font-medium text-foreground">−17.2 to 10.6 µm</span>), mostly mid-range
                    (<span className="font-medium text-foreground">49–87 µm</span>)
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Spatially mapping residuals reveals non-uniform error patterns aligned with plaque centroids
                (red stars), implying missing covariates and/or anatomical confounding.
              </p>

              <PlotFrame
                src={`${base}plots/residuals_vs_distance.html`}
                title="Spatial diagnostics (bivariate residual × distance bins)"
                size="xl"
                caption="Spatial error maps showing structured residual patterns aligned with plaque locations, indicating that gene-only models miss important spatial/anatomical factors."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Cluster-specific errors reinforce this: mean absolute residuals are largest for
                ventricular-associated cluster 15 (GnRH1-expressing glutamatergic neurons), plausibly because
                these cells occupy regions far from plaques and the model under-utilizes the full distance
                range.
              </p>

              <ResidualByClusterTable />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Top predictive genes are relatively consistent across modeling classes. Notably, **Gfap** and
                **Spag16** appear in the top 5 across all models, and three of four models rank **Gfap** as
                the single most important predictor—reinforcing earlier evidence that Gfap peaks near plaques
                and decays with distance. Other highly ranked genes include **Lyz2** (immune/glial
                association) and **Igf2**, which shows an opposite trend (peaking around ~270 µm and
                declining toward plaques).
              </p>

              <FeatureRankTable />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Agreement across linear and tree models suggests these genes carry robust, biologically
                meaningful information about plaque proximity, even if that information alone is
                insufficient for high-accuracy distance reconstruction.
              </p>
            </div>

            {/* =========================
                Additional feature modalities
               ========================= */}
            <div id="rq4-modalities" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Additional feature modalities</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To improve interpretability and performance while probing feature importance, we add
                non-transcriptomic predictors:{" "}
                <span className="font-medium text-foreground">
                  cell centroid coordinates, cell area , nucleus area and cell type (Leiden cluster ID)
                </span>
                .
              </p>

              <ul className="mt-3 space-y-2 text-sm text-muted-foreground">
                <li className="flex items-start gap-2">
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>
                      <span className="font-medium text-foreground">Genes:</span>{" "}
                      347 expression values
                    </span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Morphology:</span>{" "}
                    cell area, nucleus area
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Spatial:</span>{" "}
                    centroid coordinates
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Cluster:</span>{" "}
                    Leiden cluster identity
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We compare ordinary linear regression (OLS-like linear model) to Partial Least Squares (PLS),
                which constrains predictions through a small number of latent components optimized for
                covariance with plaque distance.
              </p>

              <ModelR2Table />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The linear model improves over gene-only linear baselines and narrows the gap to XGBoost,
                suggesting that a meaningful share of distance variance is linearly attributable to combined
                gene, spatial, morphological, and cell-type predictors. PLS is lower but extremely stable,
                consistent with capturing a dominant plaque-related axis rather than all linear variance.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To test whether the plaque-trained PLS signature reflects biologically meaningful
                plaque-centered decay (rather than arbitrary structure), we compare Tg17 and WT13 profiles
                within matched plaque-centered regions after geometric alignment.
              </p>

              <PlotFrame
                src={`${base}plots/tg17_wt13_alignment.html`}
                title="Tg17-vs-WT13 plaque-centered alignment"
                size="md"
                caption="Plaque-centered alignment between Tg17 and WT13 tissue regions (flip + translation), used to compare signature decay trends"
              />

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    Accross plaques:
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Tg17</span>{" "}
                    shows more negative correlations between signature and distance (
                    <span className="font-medium text-foreground">mean ρ ≈ −0.19</span>) and steeper
                    negative slopes (
                    <span className="font-medium text-foreground">mean slope ≈ −0.0055</span>).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">WT13</span>{" "}
                    shows near-flat or weak trends (
                    <span className="font-medium text-foreground">mean ρ ≈ −0.02</span>,{" "}
                    <span className="font-medium text-foreground">mean slope ≈ −0.0018</span>).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    FDR-corrected results indicate{" "}
                    <span className="font-medium text-foreground">3 plaques</span>{" "}
                    with significant decay in Tg but not WT, supporting a plaque-linked gradient in a subset
                    of locations.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Some plaques show decay in both Tg and WT, plausibly reflecting imperfect cross-animal
                    alignment or shared anatomical gradients rather than true pathology.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    The inner <span className="font-medium text-foreground">0–50 µm</span>{" "}
                    region has, on average, higher signature in Tg than WT, consistent with localized
                    activation near plaques.
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Two example visualizations:
              </p>

              <PlotFrame
                src={`${base}plots/signature_decay_plaque_1794.html`}
                title="Signature Decay for plaque 1794"
                size="md"
                caption="Example plaque (ID 1794): continuous signature values vs distance, illustrating plaque-centered decay behavior."
              />
              <PlotFrame
                src={`${base}plots/binned_signature_decay_plaque_1794.html`}
                title="Binned-signature Decay for plaque 1794"
                size="md"
                caption="Example plaque (ID 1794): binned signature decay vs distance, providing a more robust view of gradient shape."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Overall, the PLS signature appears biologically informative for a subset of plaques, but
                interpretation must remain cautious given alignment error and the absence of a true plaque
                ground truth in WT tissue.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We next compare linear models (Ridge/Lasso/PLS) to nonlinear models (especially gradient
                boosting) under modality ablation. Results show:
              </p>

              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Linear models reach modest accuracy (
                    <span className="font-medium text-foreground">R² ≈ 0.20–0.24</span>) when gene
                    expression is included.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>Morphology-only or spatial-only linear models perform near chance.</span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>Nonlinear models extract substantially richer structure.</span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Surprisingly, <span className="font-medium text-foreground">spatial-only</span>{" "}
                    nonlinear models can reach <span className="font-medium text-foreground">R² ≈ 0.64</span>{" "}
                    (HistGradientBoosting), exceeding even full multimodal models.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Gene-only and gene+cluster inputs improve tree models moderately (
                    <span className="font-medium text-foreground">R² ≈ 0.26–0.28</span>), while morphology
                    adds little.
                  </span>
                </li>
              </ul>

              <PlotFrame
                src={`${base}plots/modality_ablation.html`}
                title="Modality ablation"
                size="md"
                caption="Modality ablation performance across model classes: heatmap and summary bars highlighting that nonlinear models—especially spatial-only boosting—can achieve high R² driven by spatial structure."
              />
              <PlotFrame
                src={`${base}plots/best_model_per_modality.html`}
                title="Modality ablation"
                size="md"
                caption="Modality ablation performance across model classes: heatmap and summary bars highlighting that nonlinear models—especially spatial-only boosting—can achieve high R² driven by spatial structure."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The high R² from spatial-only boosting is not automatically evidence of plaque biology. A
                plausible explanation is anatomical bias: plaque deposition is not spatially uniform, and
                certain anatomical regions accumulate more plaques than others. In that case, coordinates
                predict *regional vulnerability* rather than true plaque distance.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To test whether spatial-only predictions reflect disease progression versus conserved
                anatomy, we evaluate cross-mouse stability of spatial predictions.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We compare spatial-only model outputs across all Tg and WT animals:
              </p>

              <PlotFrame
                src={`${base}plots/predicted_plaque_distance.html`}
                title="Predicted plaque distance"
                size="md"
                caption="Spatial-only prediction comparison across mice, illustrating similarity in predicted distributions despite genotype/age differences."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The model produces nearly identical prediction distributions across all mice: output
                distributions differ by at most <span className="font-medium text-foreground">~6%</span> at any
                point, which is inconsistent with a pathology-sensitive model that would be expected to
                shift with genotype and disease stage.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We quantify distributional similarity using Jensen–Shannon divergence (JSD) and corroborate
                with KS tests and ANOVA:
              </p>

              <PlotFrame
                src={`${base}plots/distribution_across_mice.html`}
                title="Distribution across mice"
                size="md"
                caption="Predicted plaque-distance score distributions across all mice under the spatial-only model, showing striking overlap across Tg and WT cohorts."
              />

              <PlotFrame
                src={`${base}plots/JS_Divergencee.html`}
                title="Jensen-Shannon divergence between mice"
                size="md"
                caption="Jensen–Shannon divergence matrix between spatial-only prediction distributions across mice; values are uniformly low, indicating near-indistinguishable outputs across genotypes and ages."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The spatial-only model primarily learns conserved tissue geometry (e.g., cortical curvature
                and laminar structure), not plaque pathology. The apparent high R² is therefore driven by
                anatomical confounding rather than disease signal.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">Findings :</p>

              <ul className="mt-3 space-y-3 text-sm text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    JSD values are mostly <span className="font-medium text-foreground">0.05–0.10</span>, only
                    slightly higher (~<span className="font-medium text-foreground">0.14–0.16</span>) for
                    comparisons involving WT-13.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    KS statistics are small (mostly <span className="font-medium text-foreground">0.02–0.08</span>)
                    despite extremely significant p-values driven by large sample sizes.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    ANOVA across Tg mice yields a very significant p-value (p ≈{" "}
                    <span className="font-medium text-foreground">3.6e−22</span>) but with trivial effect size
                    and no monotone increase with age.
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

export default RQ4Section;
