import PlotFrame from "@/components/PlotFrame";
import { Brain, CircleChevronRight } from "lucide-react";
const base = import.meta.env.BASE_URL;


const RQ5Section = () => {
  return (
    <section id="rq-5" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">
            
          {/* 

          {/* Main content */}
          <div className="space-y-12">
            <div className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                RQ5: How does the gene expression change with age for each cell type and mouse group?
              </h3>
            </div>
             {/* =========================
                RQ5
                ========================= */}
            <div id="rq5-modalities" className="space-y-4">

                <p className="text-lg text-muted-foreground leading-relaxed">
                    Because spatial-only prediction is dominated by conserved anatomy, we pivot to gene-expression-anchored, age-aware analyses. This shift is motivated by three constraints:
                </p>

                
                    <div className="rounded-2xl border border-border bg-card p-5">
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> 
                          Plaque-induced transcriptional responses are **cell-state specific** (microglia, astrocytes, and some oligodendrocyte populations can change strongly near plaques in ways spatial position alone cannot resolve).
                        </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> 
                          Direct normalization across mice is unreliable: orientation, capture area, imaging depth, and detection efficiency create batch-like distortions, and global scaling/quantile matching can suppress real gradients or introduce artifacts.  
                          </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> 
                          Common harmonization methods (Harmony, MNN, scVI) are ill-suited here: the panel is sparse, cell count is very large, and plaque-associated variance is biological signal—not batch noise to be removed.
                          </li>
                    </ul>
                    </div>
                
            </div>


            

            

            {/* =========================
                Disease status 
               ========================= */}
            <div id="rq5-disease" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">
                Disease specificity (Tg vs WT) by cluster and gene ranking
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To avoid cross-mouse normalization pitfalls while isolating cell-type-specific signals, we construct within-cluster and within-gene z-normalized signatures, then average z-scores across the <span className="font-medium text-foreground">16 PIGs</span> to obtain a single <span className="font-medium text-foreground">plaque-induced gene activation score</span> per cluster.
              </p>

              <PlotFrame
                src={`${base}plots/mean_PIG_per_mouse.html`}
                title="PIG scores per cluster per mouse"
                size="md"
                caption="Mean PIG activation score per mouse (and cluster context), derived from within-cluster, within-gene z-normalization to enable robust across-mouse comparisons."
              />

              <ul className="space-y-3 text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>
                      With this framework, we analyze two biological dimensions:
                    </span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Disease status: <span className="font-medium text-foreground">Tg vs WT</span> 
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    Age progression: <span className="font-medium text-foreground">2 → 5 → 17 months</span>
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                For each Leiden cluster, we compute a disease-specific activation score that highlights genes strongly expressed in Tg mice but minimally expressed in WT controls—capturing plaque-linked induction rather than baseline glial identity or general aging. We then rank clusters by the mean disease specificity across genes.
              </p>

              <PlotFrame
                src={`${base}plots/disease_effect_age_progression.html`}
                title="Disease effect and age progression"
                size="md"
                caption="Cluster-level summary of disease specificity (Tg − WT) alongside age progression, used to rank which cell types show the strongest plaque-linked transcriptional activation."
              />

              <PlotFrame
                src={`${base}plots/age_progression.html`}
                title="Age progression"
                size="md"
                caption="Per-cluster age progression."
              />

              
            </div>

            {/* =========================
                Age trajectories
               ========================= */}
            <div id="rq5-age" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Age trajectories: Tg (2→5→17 months) vs WT stability </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                A complementary perspective is to view the PCA structure within each cluster: WT and Tg cells separate within the same cluster, with WT occupying the bulk and Tg pushed toward the periphery. This indicates anomalous expression patterns even after conditioning on cell type.
              </p>

              <PlotFrame
                src={`${base}plots/PCA_clusters.html`}
                title="Age progression"
                size="lg"
                caption="Per-cluster age progression."
              />

              <ul className="space-y-3 text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>
                      Our analysis identifies:
                    </span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Cluster 8 (microglia)</span> as the strongest AD-specific activation cluster  
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    <span className="font-medium text-foreground">Cluster 18 (astrocytes)</span> as the next strongest
                  </span>
                </li>
              </ul>

              <p className="text-lg text-muted-foreground leading-relaxed">
                In these clusters, disease-specific genes such as **Syngr1, Gfap, and Sparcl1** show large positive specificity scores (high in Tg, mostly silent in WT), consistent with glial reactivity, complement/inflammatory remodeling, and plaque-associated activation programs.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To ensure signals are not simply driven by differences in cell-type abundance, we incorporate Tg vs WT differential expression and age progression comparisons. Tg − WT effect sizes confirm that microglia and astrocytes exhibit the largest positive shifts in PIG expression.
              </p>

              <PlotFrame
                src={`${base}plots/ad_specific_genes.html`}
                title="AD-specific genes"
                size="sm"
                caption="Heatmap of AD-specific genes (high in Tg, low in WT), emphasizing that plaque-linked activation is concentrated in specific clusters and genes."
              />

              <PlotFrame
                src={`${base}plots/top_genes_per_glial.html`}
                title="Top genes per glial cluster"
                size="sm"
                caption="Top differential genes per glial cluster (logFC), highlighting microglial and astrocytic programs most altered in Tg relative to WT."
              />

              <ul className="space-y-3 text-muted-foreground">
                <li>
                  <p className="text-lg text-muted-foreground leading-relaxed">
                    <span>
                      Finally, age trajectories show divergence:
                    </span>
                  </p>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    In Tg mice, microglia and astrocytes show *<span className="font-medium text-foreground">monotone increases</span> from <span className="font-medium text-foreground">2 to 5 to 17 months</span>.  
                  </span>
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span>
                    WT mice remain stable or slightly decline.
                  </span>
                </li>
              </ul>

              <PlotFrame
                src={`${base}plots/age_progression_wt_tg.html`}
                title="Age progression WT vs Tg"
                size="lg"
                caption="Age progression trajectories of cluster-level activation for WT vs Tg, showing AD-specific, age-progressive glial activation in Tg animals"
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Taken together, RQ6 supports a coherent synthesis: specific glial clusters—particularly microglia (cluster 8) and astrocytes (cluster 18)—undergo robust and progressive transcriptional activation driven by amyloid pathology. These signatures intensify with age in Tg mice but remain absent in age-matched WT animals. Compared to spatial-only modeling, this gene-level, age-resolved approach yields a more stable and pathology-driven understanding of how glial states evolve around Aβ plaques..
              </p>
            </div>
          </div>


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
                                    How does the gene expression change with age for each cell type and mouse group?
                                    </div>
                                </div>
                                </div>
            
                                {/* Navigation (identical hover / spacing / typography) */}
                                <nav className="mt-6 space-y-2 text-sm">
                                <a
                                    href="#rq5-disease"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                    Disease specificity
                                </a>
            
                                <a
                                    href="#rq5-age"
                                    className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                                >
                                    Age trajectories
                                </a>
          
                                </nav>
            
                                {/* Footer hint */}
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

export default RQ5Section;
