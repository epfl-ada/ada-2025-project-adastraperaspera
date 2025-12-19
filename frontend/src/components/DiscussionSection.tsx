import { Brain, CircleChevronRight } from "lucide-react";

const DiscussionSection = () => {
  return (
    <section id="discussion" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Sticky LEFT panel + main content RIGHT */}
        <div className="grid gap-10 items-start lg:grid-cols-[320px_1fr]">
          {/* Sticky panel (left) */}
          <aside className="lg:sticky lg:top-24 lg:self-start">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Brain className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">Discussion</div>
                  <div className="text-xs text-muted-foreground">
                    Conclusions, caveats, and implications
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#discussion-6a"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Robust conclusions
                </a>
                <a
                  href="#discussion-6b"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Confounding & risks
                </a>
                <a
                  href="#discussion-6c"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Implications
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Main discussion */}
          <div className="space-y-12">
            <div className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">Discussion</h3>
            </div>

            {/* 6.a */}
            <div
              id="discussion-6a"
              className="rounded-2xl border border-border bg-card p-6 shadow-sm space-y-4"
            >
              <h4 className="text-xl font-semibold text-foreground">
                What we can conclude robustly (and what we cannot)
              </h4>

              <ul className="mt-3 space-y-3 text-sm text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      Plaque proximity reshapes cell-type composition
                    </span>
                    , with immune/vascular/astrocytic enrichment near plaques and broad neuronal depletion (RQ1).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      PIG gradients are strongly linked to composition shifts
                    </span>
                    , especially immune enrichment and neuronal depletion, but remain distance-associated even after accounting for cell type (RQ2).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      All 16 PIGs show significant plaque-proximal elevation
                    </span>
                    , with Gfap exhibiting the strongest spatial gradient and Cxcl10 appearing nearly flat due to extreme zero inflation (RQ3).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      Gene expression contains limited but real distance information
                    </span>{" "}
                    (best gene-only test R² ≈ 0.26 with XGBoost), and error structure is spatially patterned and heteroscedastic (RQ4).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      Spatial-only models can be misleadingly strong
                    </span>
                    : cross-mouse diagnostics show they largely learn conserved anatomy rather than plaque pathology (RQ4).
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">
                      Age- and genotype-aware z-normalized signatures reveal AD-specific, age-progressive glial activation
                    </span>
                    , concentrated in microglia and astrocytes (RQ5).
                  </span>
                </li>
              </ul>
            </div>

            {/* 6.b */}
            <div
              id="discussion-6b"
              className="rounded-2xl border border-border bg-card p-6 shadow-sm space-y-4"
            >
              <h4 className="text-xl font-semibold text-foreground">
                Confounding, alignment error, sparsity/zero inflation, and interpretation risks
              </h4>

              <ul className="mt-3 space-y-3 text-sm text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Alignment error:</span>{" "}
                    IF-to-morphology RMSE is small (3.2 µm) but nonzero; fine-scale (single-digit µm) conclusions remain sensitive.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Zero inflation and sparsity:</span>{" "}
                    genes like Cxcl10 illustrate that statistical significance can coexist with minimal practical effect size due to near-all-zero distributions.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Anatomical confounding:</span>{" "}
                    plaque density varies by region; any model using coordinates must be treated as potentially learning anatomy rather than pathology.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Cross-mouse comparability:</span>{" "}
                    batch-like distortions make global cross-mouse normalization risky; the chosen within-cluster z-score approach mitigates but does not eliminate all comparability concerns.
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Causal direction:</span>{" "}
                    composition shifts and PIG changes co-occur; while regression and correlation help disentangle them, they do not establish causality.
                  </span>
                </li>
              </ul>
            </div>

            {/* 6.c */}
            <div
              id="discussion-6c"
              className="rounded-2xl border border-border bg-card p-6 shadow-sm space-y-4"
            >
              <h4 className="text-xl font-semibold text-foreground">
                Implications for target selection and drug development relevance
              </h4>

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

export default DiscussionSection;
