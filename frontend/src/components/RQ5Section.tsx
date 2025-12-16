import PlotFrame from "@/components/PlotFrame";

const base = "/"; // adjust if needed

const RQ5Section = () => {
  return (
    <section id="rq-5" className="py-24 bg-background">
      <div className="container mx-auto px-6 space-y-20">

        {/* =========================
           Header
        ========================= */}
        <div className="text-center space-y-4">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground">
            RQ5 — Feature Importance in Plaque Distance Modeling
          </h2>
          <p className="text-lg text-muted-foreground max-w-3xl mx-auto">
            Which modalities truly carry predictive signal for plaque proximity, and which
            reflect anatomical bias rather than disease biology?
          </p>
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
           Linear Models Table
        ========================= */}
        <div className="space-y-4">
          <h3 className="text-2xl font-bold text-foreground">Linear Models</h3>

          <p className="text-muted-foreground leading-relaxed max-w-4xl">
            We benchmark ordinary least squares (OLS) and partial least squares (PLS) regression
            using all four modalities. PLS constrains the model to a small number of latent
            components aligned with plaque distance.
          </p>

          <div className="overflow-x-auto">
            <table className="w-full border border-border rounded-xl text-sm">
              <thead className="bg-muted">
                <tr>
                  <th className="p-3 text-left">Model</th>
                  <th className="p-3 text-right">Train R²</th>
                  <th className="p-3 text-right">Test R²</th>
                </tr>
              </thead>
              <tbody>
                <tr className="border-t border-border">
                  <td className="p-3">Linear (OLS)</td>
                  <td className="p-3 text-right">0.25</td>
                  <td className="p-3 text-right">0.24</td>
                </tr>
                <tr className="border-t border-border">
                  <td className="p-3">PLS</td>
                  <td className="p-3 text-right">0.19</td>
                  <td className="p-3 text-right">0.19</td>
                </tr>
              </tbody>
            </table>
          </div>

          <p className="text-muted-foreground leading-relaxed max-w-4xl">
            The OLS model achieves a test R² of ~0.24, showing that a purely linear combination
            of genes, morphology, spatial coordinates, and cluster identity explains roughly
            one quarter of the variance. The small train–test gap indicates good generalization.
            PLS captures only the dominant plaque-related axis and therefore achieves lower but
            highly stable performance.
          </p>
        </div>

        {/* =========================
           Modality Ablation
        ========================= */}
        <div className="space-y-10">
          <h3 className="text-2xl font-bold text-foreground">Modality Ablation</h3>

          <PlotFrame
            src={`${base}plots/modality_ablation.html`}
            title="Modality ablation heatmap"
            caption="Predictive performance across models and feature subsets."
          />

          <PlotFrame
            src={`${base}plots/best_model_per_modality.html`}
            title="Best model per modality"
            caption="Nonlinear models outperform linear ones, especially for spatial inputs."
          />

          <p className="text-muted-foreground leading-relaxed max-w-5xl">
            Ablation analysis reveals strong modality dependence. Linear models perform near
            chance without gene expression, while nonlinear models extract richer structure.
            Surprisingly, spatial-only models achieve the highest R² (~0.64), exceeding even
            the full multimodal model—suggesting strong anatomical bias rather than true distance
            prediction.
          </p>
        </div>

        {/* =========================
           Diagnostics: Spatial Bias
        ========================= */}
        <div className="space-y-10">
          <h3 className="text-2xl font-bold text-foreground">Diagnostics: Spatial Bias</h3>

          <PlotFrame
            src={`${base}plots/predicted_plaque_distance.html`}
            title="Spatial predictions across mice"
            caption="Predicted plaque distance distributions are nearly identical across mice."
          />

          <PlotFrame
            src={`${base}plots/distribution_across_mice.html`}
            title="Distribution across mice"
            caption="Spatial-only model outputs show minimal variation by genotype or age."
          />

          <PlotFrame
            src={`${base}plots/JS_Divergence.html`}
            title="Jensen–Shannon divergence"
            caption="Predicted distributions are highly similar across all mice."
          />

          <p className="text-muted-foreground leading-relaxed max-w-5xl">
            These diagnostics confirm that the spatial-only model exploits conserved tissue
            geometry rather than plaque-related pathology. Although statistically detectable,
            differences between mice are biologically negligible and do not track disease
            progression.
          </p>
        </div>

        {/* =========================
           Shift to Gene-Based Strategy
        ========================= */}
        <div className="space-y-10">
          <h3 className="text-2xl font-bold text-foreground">Gene-Based, Age-Aware Strategy</h3>

          <p className="text-muted-foreground leading-relaxed max-w-5xl">
            To avoid anatomical confounding and unsafe cross-mouse normalization, we adopt a
            within-gene and within-cluster z-normalization strategy. Averaging z-scores across
            the 16 PIGs yields a unified plaque-induced activation score per cluster.
          </p>

          <PlotFrame
            src={`${base}plots/mean_PIG_per_mouse.html`}
            title="Mean PIG activation per mouse"
            caption="Cluster-level plaque-induced activation without cross-mouse normalization."
          />

          <PlotFrame
            src={`${base}plots/disease_effect_age_progression.html`}
            title="Disease effect vs age"
            caption="AD-specific activation separates Tg from WT mice within clusters."
          />

          <PlotFrame
            src={`${base}plots/age_progression.html`}
            title="Age progression"
            caption="Plaque-induced activation increases monotonically with age in Tg mice."
          />

          <PlotFrame
            src={`${base}plots/PCA_clusters.html`}
            title="Cluster-level PCA"
            caption="WT and Tg mice separate within each cell type."
          />

          <PlotFrame
            src={`${base}plots/ad_specific_genes.html`}
            title="AD-specific genes"
            caption="Microglia and astrocytes show the strongest disease-specific activation."
          />

          <PlotFrame
            src={`${base}plots/top_genes_per_glial.html`}
            title="Top glial genes"
            caption="Differential expression confirms disease-driven activation."
          />

          <PlotFrame
            src={`${base}plots/age_progression_wt_tg.html`}
            title="Age progression WT vs Tg"
            caption="Progressive activation occurs only in Tg mice."
          />
        </div>

        {/* =========================
           Conclusion
        ========================= */}
        <div className="rounded-2xl border border-border bg-card p-6">
          <p className="text-muted-foreground leading-relaxed max-w-5xl">
            In summary, spatial-only modeling yields deceptively high accuracy due to anatomical
            bias, not disease signal. In contrast, gene-level and age-resolved analyses uncover
            robust, AD-specific, and progressive glial activation localized to biologically
            relevant cell types. This strategy provides a stable and interpretable framework for
            understanding plaque-associated transcriptional remodeling.
          </p>
        </div>

      </div>
    </section>
  );
};

export default RQ5Section;
