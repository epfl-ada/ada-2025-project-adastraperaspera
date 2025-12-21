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
                Predicting gene expression only on the cell coordinate has clear limitations:
              </p>


              <div className="rounded-2xl border border-border bg-card p-5">
                <ul className="space-y-3 text-muted-foreground">
                  <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />Plaque-dependent gene expression changes vary by cell type. For instance, microglia, astrocytes, and oligodendrocytes show much stronger changes near plaques than other cell types.
                  </li>
                  <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Aligning the morphology images across mice is complicated and introduces a significant mismatch due to the individual variations in anatomy and brain proportions. Thus, a coordinate-only model is unlikely to generalize well to unseen mice                  </li>
                </ul>
              </div>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To overcome these issues, we take into account the mouse type (Wt/Tg) and age in the following section.
              </p>

            </div>






            {/* =========================
                Disease status 
               ========================= */}
            <div id="rq5-disease" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">
                Disease specificity (Tg vs WT) by cluster and gene ranking
              </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To examine the per-mouse variation in gene expression, we compute a z-normalized signature which is the log1p-transformed expression level averaged across 16 PIGs. The signature is computed for each mouse type and cell type pair (Tg/Wt x Leiden cluster). The results reveal that the average PIG expression in clusters 2, 5, 17 (Hypothalamic GABAergic, Thalamic Glutaergic, and Medulla GABAergic neurons, respectively) is much lower in oldest Tg mice compared to others. This indicates that these cell types can carry key importance in the development of AD.
              </p>

              <PlotFrame
                src={`${base}plots/mean_PIG_per_mouse.html`}
                title="PIG scores per cluster per mouse"
                size="md"
                caption="Mean PIG activation score per mouse (and cluster context), derived from within-cluster, within-gene z-normalization to enable robust across-mouse comparisons."
              />

            {/*<ul className="space-y-3 text-muted-foreground">
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
              </ul> */}

              <p className="text-lg text-muted-foreground leading-relaxed">
                We will now analyze more closely how exactly the gene expression signature changes with genotype (Tg/Wt) and with age. This framework reveals plaque-induced changes in expression since it operates with general Wt aging as the baseline.
                To this end, we will compute the difference in signatures between age-matched Tg/Wt mice and compute the correlation between this difference and age. We encode the absolute value of this relationship in oY while representing the sign as the circle radius. In other words, circles on the right half of the canvas represent cell types whose average PIG expression increases with age. Meanwhile, the height represents the strength of association with age. We can note that the cluster 18 (Pons Glutaergic neurons) is the right uppermost circle corresponding to the largest age-progressive relative increase of PIG expression. Meanwhile, on the other end of the spectrum, we have cluster 6 (vascular cells) with the largest age-progressive relative decrease in PIG expression. This result can be due to opposite effects plaques have on these cell types as plaque proximity is known to induce neural death but vascular enrichment.
              </p>

              <PlotFrame
                src={`${base}plots/disease_effect_age_progression.html`}
                title="Disease effect and age progression"
                size="md"
                caption="Cluster-level summary of disease specificity (Tg − WT) alongside age progression, used to rank which cell types show the strongest plaque-linked transcriptional activation."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The same pattern can be seen in the next line plot where, in Tg mice, cluster 6 shows a steep decline with age whereas cluster 18 shows a steep increase.
              </p>

              <PlotFrame
                src={`${base}plots/age_progression.html`}
                title="Age progression"
                size="md"
                caption="Per-cluster age progression."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Next, we will look at the PCA-based vizualisation of each cell type. We can see that within each cluster, Wt and Tg cells are somewhat separated. Typically, the Wt cells compose the bulk of the point cloud whereas the Tg cells appear on the periphery. This shows systematic differences in the expression patterns between Tg and Wt mice that appear in most cell types.
              </p>

              <PlotFrame
                src={`${base}plots/PCA_clusters.html`}
                title="Age progression"
                size="lg"
                caption="Per-cluster age progression."
              />

              


            </div>

            {/* =========================
                Age trajectories
               ========================= */}
            <div id="rq5-age" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Age trajectories: Tg (2→5→17 months) vs WT stability </h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We will now review the age trajectories in both mouse types. The following figure reveals that in transgenic mice, cluster 8 (immune cells) shows a monotone increase in average PIG expression from 2 to 5 to 17 months. The same cell type in the Wt mice remains stable.
              </p>

              <PlotFrame
                src={`${base}plots/age_progression_wt_tg.html`}
                title="Age progression WT vs Tg"
                size="lg"
                caption="Age progression trajectories of cluster-level activation for WT vs Tg, showing AD-specific, age-progressive glial activation in Tg animals"
              />

              {/*<ul className="space-y-3 text-muted-foreground">
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
                In these clusters, disease-specific genes such as <span className="font-medium text-foreground">Syngr1, Gfap, and Sparcl1</span> show large positive specificity scores (high in Tg, mostly silent in WT), consistent with glial reactivity, complement/inflammatory remodeling, and plaque-associated activation programs.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To ensure signals are not simply driven by differences in cell-type abundance, we incorporate Tg vs WT differential expression and age progression comparisons. Tg − WT effect sizes confirm that microglia and astrocytes exhibit the largest positive shifts in PIG expression.
              </p>*/}
              
              <p className="text-lg text-muted-foreground leading-relaxed">
                Examining the heatmap of AD-specific genes, we can see that <span className="font-medium text-foreground">Cluster 15 (Hypothalamic Gnrh1, Glutaergic)</span> is the strongest AD-specific activation cluster with <span className="font-medium text-foreground">Cluster 6 (vascular)</span> as the next strongest. In these clusters, disease-specific genes such as <span className="font-medium text-foreground">Syngr1, and Sparcl1</span> show large positive specificity scores (high in Tg, mostly silent in WT). This matches the existing knowledge on inflammatory remodeling pathways in AD as well as plaque-induced gene expression programs.
              </p>

              <PlotFrame
                src={`${base}plots/ad_specific_genes.html`}
                title="AD-specific genes"
                size="sm"
                caption="Heatmap of AD-specific genes (high in Tg, low in WT), emphasizing that plaque-linked activation is concentrated in specific clusters and genes."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Shifting our focus to glial cells, we can see that the genes C3, Nme8, and Lyz2 show the biggest relative increase in expression in Tg relative to Wt.
              </p>

              <PlotFrame
                src={`${base}plots/top_genes_per_glial.html`}
                title="Top genes per glial cluster"
                size="sm"
                caption="Top differential genes per glial cluster (logFC), highlighting microglial and astrocytic programs most altered in Tg relative to WT."
              />

              

              <p className="text-lg text-muted-foreground leading-relaxed">
                Overall, our findings demonstrate that clusters 6 and 15 undergo progressive changes in gene expression patterns as a result of amyloid pathology. These signatures become increasingly more pronounced with age in Tg mice but not in the age-matched Wt mice.
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
