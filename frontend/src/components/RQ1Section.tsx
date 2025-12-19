import PlotFrame from "@/components/PlotFrame";
import { Brain, CircleChevronRight} from 'lucide-react'

const base = import.meta.env.BASE_URL;

const NEG_ROWS = [
  { cluster: "Immune", slope: "−0.01427", p0: "0.159", p100: "0.043", change: "−72.7%" },
  { cluster: "Vascular", slope: "−0.01116", p0: "0.100", p100: "0.035", change: "−64.9%" },
  { cluster: "Pons, Gluta", slope: "−0.00772", p0: "0.015", p100: "0.007", change: "−53.4%" },
  {
    cluster: "Intra/Extratelencephalic, Gluta",
    slope: "−0.00361",
    p0: "0.075",
    p100: "0.054",
    change: "−28.7%",
  },
  { cluster: "Cortex medial, GABA", slope: "−0.00282", p0: "0.041", p100: "0.031", change: "−23.8%" },
  { cluster: "Astrocyte", slope: "−0.00114", p0: "0.136", p100: "0.123", change: "−9.5%" },
];

const POS_ROWS = [
  { cluster: "Hypothalamic Gnrh1, Gluta", slope: "+0.02346", p0: "0.002", p100: "0.019", change: "+927.0%" },
  { cluster: "Medulla, GABA", slope: "+0.01815", p0: "0.002", p100: "0.012", change: "+508.2%" },
  { cluster: "Cerebral LGE, GABA", slope: "+0.00781", p0: "0.016", p100: "0.035", change: "+114.3%" },
  { cluster: "Olfactory bulb, Gluta", slope: "+0.00678", p0: "0.011", p100: "0.021", change: "+94.9%" },
  { cluster: "Dentate, Gluta", slope: "+0.00591", p0: "0.014", p100: "0.025", change: "+78.5%" },
  { cluster: "Corticothalamic, Gluta", slope: "+0.00384", p0: "0.029", p100: "0.042", change: "+44.8%" },
  { cluster: "Hypothalamic medial, Gluta", slope: "+0.00256", p0: "0.033", p100: "0.042", change: "+27.9%" },
  { cluster: "Oligodendrocyte", slope: "+0.00208", p0: "0.161", p100: "0.191", change: "+18.7%" },
];

function StatTable({
  title,
  rows,
}: {
  title: string;
  rows: { cluster: string; slope: string; p0: string; p100: string; change: string }[];
}) {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <div className="flex items-center justify-between gap-3 mb-4">
        <h4 className="text-base font-semibold text-foreground">{title}</h4>
      </div>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Cluster</th>
              <th className="py-2 pr-4 text-right font-medium">Slope βₖ</th>
              <th className="py-2 pr-4 text-right font-medium">p(0)</th>
              <th className="py-2 pr-4 text-right font-medium">p(100)</th>
              <th className="py-2 text-right font-medium">Relative change</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => (
              <tr key={r.cluster} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground">{r.cluster}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.slope}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.p0}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.p100}</td>
                <td className="py-2 text-right tabular-nums">{r.change}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

const RQ1Section = () => {
  return (
    <section id="rq-1" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + scrollable content right */}
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">
          
          {/* Analysis content */}
          <div className="space-y-12">
            {/* Intro / method */}
            <div id="rq1-cluster-dist" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">RQ1 : How does cell type composition change in plaque proximity?</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We model the relationship between <span className="font-medium text-foreground">cluster membership</span> and <span className="font-medium text-foreground">distance to the nearest plaque</span> using logistic regression. After fitting, we find statistically significant coefficients for <span className="font-medium text-foreground">14 of 19 clusters</span>, using a <span className="font-medium text-foreground">Bonferroni-adjusted p-value threshold of 0.01</span>.
              </p>
            </div>
            
             {/*<div className="rounded-2xl border border-border bg-card p-5">*/}
              <ul className="mt-3 space-y-2 text-sm text-muted-foreground">
                <p className="text-lg text-muted-foreground leading-relaxed">
                To make coefficients interpretable, we translate them into:
                </p>
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">p(0):</span>{" "}
                  baseline probability of observing a cluster at the plaque surface
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">p(100):</span>{" "}
                  probability of observing a cluster at{" "}
                  <span className="font-medium text-foreground">100 µm</span> from the plaque
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">p(100) − p(0):</span>{" "}
                  relative change over{" "}
                  <span className="font-medium text-foreground">100 µm</span> away from the plaque
                </li>
              </ul>
             {/*</div>*/}
            <p className="text-lg text-muted-foreground leading-relaxed">
              Key effects (selected):
            </p>

            {/* Tables */}
            <div className="space-y-6">

              <div className="grid gap-6 lg:grid-cols-2">
                <StatTable title="Clusters depleted with distance (negative slopes; enriched near plaques)" rows={NEG_ROWS} />
                <StatTable title="Clusters enriched with distance (positive slopes; depleted near plaques)" rows={POS_ROWS} />
              </div>

              <p className="text-lg text-muted-foreground leading-relaxed">
                These coefficients are summarized visually below.
              </p>

              <PlotFrame
                src={`${base}plots/slopes_types.html`}
                title="Plaque distance by cell type"
                size="md"
                caption="Regression slopes for significant clusters."
              />
            </div>

            <p className="text-lg text-muted-foreground leading-relaxed">
                Overall, the results show a clear composition shift near plaques:
            </p>
            <div className="rounded-2xl border border-border bg-card p-5">
                    <ul className="mt-3 space-y-2 text-sm text-muted-foreground">
                      <li className="flex items-start gap-2">
                        <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                        <span className="font-medium text-foreground">
                          Neuronal cell types are generally depleted
                        </span>{" "}
                        around plaques.
                      </li>

                      <li className="flex items-start gap-2">
                        <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                        <span className="font-medium text-foreground">
                          Immune, vascular, and astrocytic
                        </span>{" "}
                        populations are{" "}
                        <span className="font-medium text-foreground">enriched</span>{" "}
                        at the smallest distances.
                      </li>
                    </ul>
            </div>

            <p className="text-lg text-muted-foreground leading-relaxed">
              This aligns with known AD mechanisms: plaques are associated with neuronal
              degeneration, reactive astrocytosis and microglial activation, and vascular
              remodeling.
            </p>

            <p className="text-lg text-muted-foreground leading-relaxed">
              A particularly strong effect appears for{" "}
              <span className="font-medium text-foreground">
                Hypothalamic GnRH1-expressing glutamatergic neurons (cluster 15)
              </span>
              , whose frequency increases{" "}
              <span className="font-medium text-foreground">~10×</span> over{" "}
              <span className="font-medium text-foreground">100 µm</span>. Spatial overlays
              suggest these cells largely reside in ventricle cavities and are therefore
              typically far from plaques (cerebrospinal fluid regions are effectively
              plaque-free), explaining the steep distance-associated rise.
            </p>

            <p className="text-lg text-muted-foreground leading-relaxed">
              We also observe a more surprising pattern: several neuronal clusters (e.g.,
              intra-/extratelencephalic glutamatergic, pons glutamatergic, cortex-medial
              GABAergic) appear relatively enriched near plaques compared to other neurons.
              One plausible hypothesis is differential resilience: these neurons may be less
              vulnerable to plaque-associated damage, yielding a weaker depletion pattern
              than other neuronal populations.
            </p>

            {/* Frequency vs distance */}
            <div id="rq1-frequency" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Cluster frequency vs distance</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To complement the regression summary, we examine frequency vs. distance using binned distance profiles.
              </p>

              <PlotFrame
                src={`${base}plots/cluster_frequency_distance_to_plaque.html`}
                title="Cluster frequency vs. distance to plaque"
                size="lg"
                caption="Distance to plaque vs % of cells in each bin (interactive)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                  Several clusters show changes primarily at extreme distances rather than gradual shifts
              </p>
                <div className="rounded-2xl border border-border bg-card p-5">
                
                      <ul className="mt-3 space-y-2 text-sm text-muted-foreground">
                        <li className="flex items-start gap-2">
                          <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                          <span className="font-medium text-foreground">
                            Cluster 0 (oligodendrocyte lineage)
                          </span>{" "}
                          shows a clear upward trend away from plaques.
                        </li>

                        <li className="flex items-start gap-2">
                          <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                          <span className="font-medium text-foreground">Cluster 15</span>{" "}
                          increases sharply only after{" "}
                          <span className="font-medium text-foreground">~113 µm</span>.
                        </li>

                        <li className="flex items-start gap-2">
                          <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                          <span className="font-medium text-foreground">Cluster 14</span>{" "}
                          declines sharply after{" "}
                          <span className="font-medium text-foreground">~138 µm</span>.
                        </li>

                        <li className="flex items-start gap-2">
                          <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                          <span className="font-medium text-foreground">Cluster 8 (immune)</span>{" "}
                          drops sharply after the first bin, then levels off.
                        </li>
                      </ul>

                </div>
                <p className="text-lg text-muted-foreground leading-relaxed">
                  These non-linearities motivate correlation and regression analyses that do not assume strict linear response across distance.
                </p>

            </div>

            {/* Marker enrichment heatmap 
            <div id="rq1-markers" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Marker gene enrichment</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Looking at the marker gene enrichment diagram, we see strong over- and under-expression patterns that reflect
                the relationship between brain morphology and gene expression programs.
              </p>

              <PlotFrame
                src={`${base}plots/expression_per_cluster.html`}
                title="Marker gene enrichment (z-scored across clusters)"
                size="lg"
                caption="Z-scored marker expression per cluster (interactive)."
              />




            </div>*/}
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
                        Research Question 1
                        </div>
                        <div className="text-xs text-muted-foreground">
                        How does cell type composition change in plaque proximity?
                        </div>
                    </div>
                    </div>

                    {/* Navigation (identical hover / spacing / typography) */}
                    <nav className="mt-6 space-y-2 text-sm">
                    <a
                        href="#rq1-cluster-dist"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Depletion & enrichment by cluster
                    </a>

                    <a
                        href="#rq1-frequency"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Frequency vs distance
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

export default RQ1Section;
