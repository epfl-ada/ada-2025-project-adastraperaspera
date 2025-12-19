import { CheckCircle2, AlertTriangle, Target } from "lucide-react";

const DiscussionSection2 = () => {
  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">
        <div className="max-w-5xl mx-auto space-y-12">
          {/* Header */}
          <div className="text-center space-y-4">
            <h2 className="text-3xl md:text-4xl font-bold text-foreground">
              Discussion
            </h2>
            <p className="text-lg text-muted-foreground">
              Integrated interpretation, limitations, and biological relevance
            </p>
          </div>

          {/* 6.a */}
          <div className="rounded-2xl border border-border bg-card p-8 shadow-sm">
            <div className="flex items-center gap-3 mb-6">
              <CheckCircle2 className="w-6 h-6 text-primary" />
              <h3 className="text-xl font-semibold text-foreground">
                6.a. What we can conclude robustly (and what we cannot)
              </h3>
            </div>

            <div className="space-y-4">
              <div className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10">
                <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">
                    Plaque proximity reshapes cell-type composition
                  </span>
                  , with immune/vascular/astrocytic enrichment near plaques and broad neuronal depletion (RQ1).
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10">
                <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">
                    PIG gradients are strongly linked to composition shifts
                  </span>
                  , especially immune enrichment and neuronal depletion, but remain distance-associated even after accounting for cell type (RQ2).
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10">
                <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">
                    All 16 PIGs show significant plaque-proximal elevation
                  </span>
                  , with Gfap exhibiting the strongest spatial gradient and Cxcl10 appearing nearly flat due to extreme zero inflation (RQ3).
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10">
                <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">
                    Gene expression contains limited but real distance information
                  </span>{" "}
                  (best gene-only test R² ≈ 0.26 with XGBoost), and error structure is spatially patterned and heteroscedastic (RQ4).
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10">
                <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">
                    Spatial-only models can be misleadingly strong
                  </span>
                  : cross-mouse diagnostics show they largely learn conserved anatomy rather than plaque pathology (RQ4).
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10">
                <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">
                    Age- and genotype-aware z-normalized signatures reveal AD-specific, age-progressive glial activation
                  </span>
                  , concentrated in microglia and astrocytes (RQ5).
                </p>
              </div>
            </div>
          </div>

          {/* 6.b */}
          <div className="rounded-2xl border border-border bg-card p-8 shadow-sm">
            <div className="flex items-center gap-3 mb-6">
              <AlertTriangle className="w-6 h-6 text-chart-4" />
              <h3 className="text-xl font-semibold text-foreground">
                6.b Confounding, alignment error, sparsity/zero inflation, and interpretation risks
              </h3>
            </div>

            <div className="space-y-4">
              <div className="flex items-start gap-4 p-4 rounded-xl bg-muted/40 border border-border">
                <AlertTriangle className="w-6 h-6 flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">Alignment error:</span> IF-to-morphology RMSE is small (3.2 µm) but nonzero; fine-scale (single-digit µm) conclusions remain sensitive.
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-muted/40 border border-border">
                <AlertTriangle className="w-6 h-6 flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">Zero inflation and sparsity:</span> genes like Cxcl10 illustrate that statistical significance can coexist with minimal practical effect size due to near-all-zero distributions.
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-muted/40 border border-border">
                <AlertTriangle className="w-6 h-6 flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">Anatomical confounding:</span> plaque density varies by region; any model using coordinates must be treated as potentially learning anatomy rather than pathology.
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-muted/40 border border-border">
                <AlertTriangle className="w-6 h-6 flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">Cross-mouse comparability:</span> batch-like distortions make global cross-mouse normalization risky; the chosen within-cluster z-score approach mitigates but does not eliminate all comparability concerns.
                </p>
              </div>

              <div className="flex items-start gap-4 p-4 rounded-xl bg-muted/40 border border-border">
                <AlertTriangle className="w-6 h-6 flex-shrink-0 mt-0.5" />
                <p className="text-foreground">
                  <span className="font-medium">Causal direction:</span> composition shifts and PIG changes co-occur; while regression and correlation help disentangle them, they do not establish causality.
                </p>
              </div>
            </div>
          </div>

          {/* 6.c */}
          <div className="rounded-2xl border border-border bg-card p-8 shadow-sm">
            <div className="flex items-center gap-3 mb-6">
              <Target className="w-6 h-6 text-chart-2" />
              <h3 className="text-xl font-semibold text-foreground">
                6.c. Implications for target selection and drug development relevance
              </h3>
            </div>

            <div className="rounded-xl bg-primary/5 border border-primary/10 p-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Despite these limitations, the combined evidence supports a biologically consistent plaque-centered narrative: plaques are surrounded by activated glial niches (microglia/astrocytes) with elevated plaque-induced genes, accompanied by cell-type redistribution and distance-dependent decay. This quantitative framing is directly relevant for identifying plausible therapeutic targets that are (i) proximal to plaques, (ii) cell-type-specific, and (iii) progressive with pathology and age.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default DiscussionSection2;
