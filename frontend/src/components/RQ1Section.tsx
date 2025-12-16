import PlotFrame from "@/components/PlotFrame";
import { Brain } from 'lucide-react'

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
    <div id="rq-1" className="rounded-2xl border border-border bg-card p-6 shadow-sm">
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
                        href="#rq1-method"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Method
                    </a>

                    <a
                        href="#rq1-slopes"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Key coefficients
                    </a>

                    <a
                        href="#rq1-frequency"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Frequency vs distance
                    </a>

                    <a
                        href="#rq1-markers"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Marker enrichment
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
              <div className="text-xs uppercase tracking-wide text-muted-foreground">Research Question 1</div>
              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                How does cell type composition change in plaque proximity?
              </h2>

              <div className="mt-5 space-y-3 text-sm">
                <a
                  href="#rq1-method"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Method</div>
                  <div className="text-muted-foreground">Logistic regression vs distance</div>
                </a>

                <a
                  href="#rq1-slopes"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Key coefficients</div>
                  <div className="text-muted-foreground">p(0), p(100), Δp</div>
                </a>

                <a
                  href="#rq1-frequency"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Frequency vs distance</div>
                  <div className="text-muted-foreground">Binned cluster proportions</div>
                </a>

                <a
                  href="#rq1-markers"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Marker enrichment</div>
                  <div className="text-muted-foreground">Z-scored heatmap</div>
                </a>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: use the links above, or just scroll the analysis on the right.
              </div>
            </div>
          </aside>Sticky RQ panel */}

          {/* Analysis content */}
          <div className="space-y-12">
            {/* Intro / method */}
            <div id="rq1-method" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">RQ1 Analysis</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We use logistic regression to model the relationship between cluster membership and distance to the nearest
                plaque. As a result of running regression, we found statistically significant coefficients for 14 out of 19
                clusters. We used a Bonferroni-adjusted p-value threshold of 0.01.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We interpret the regression coefficients as (1) baseline probability of seeing a cluster at the plaque
                surface <span className="font-medium text-foreground">p(0)</span>, (2) probability at 100 µm{" "}
                <span className="font-medium text-foreground">p(100)</span>, and (3) relative change{" "}
                <span className="font-medium text-foreground">p(100) − p(0)</span>.
              </p>
            </div>

            {/* Tables */}
            <div id="rq1-slopes" className="space-y-6">
              <div className="flex items-end justify-between gap-6">
                <h3 className="text-xl font-semibold text-foreground">Significant slopes</h3>
                <div className="text-sm text-muted-foreground">Bonferroni-adjusted p-value threshold = 0.01</div>
              </div>

              <div className="grid gap-6 lg:grid-cols-2">
                <StatTable title="Enriched near plaques (negative slope)" rows={NEG_ROWS} />
                <StatTable title="Enriched far from plaques (positive slope)" rows={POS_ROWS} />
              </div>

              <p className="text-lg text-muted-foreground leading-relaxed">
                The results reveal changes in cell composition with varying distance to the nearest plaque. Most neural
                clusters are depleted around plaques, while astrocytes, vascular, and immune subtypes are enriched near
                plaques.
              </p>

              <PlotFrame
                src={`${base}plots/slopes_types.html`}
                title="Plaque distance by cell type"
                size="md"
                caption="Regression slopes for significant clusters (interactive)."
              />
            </div>

            {/* Frequency vs distance */}
            <div id="rq1-frequency" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Cluster frequency vs distance</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We next examine how cluster frequency varies with binned distance to the nearest plaque. Some clusters
                change mainly at extreme distances, while immune cells show a sharp drop after the first bin and then level
                off.
              </p>

              <PlotFrame
                src={`${base}plots/cluster_frequency_distance_to_plaque.html`}
                title="Cluster frequency vs. distance to plaque"
                size="lg"
                caption="Distance to plaque vs % of cells in each bin (interactive)."
              />
            </div>

            {/* Marker enrichment heatmap */}
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
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ1Section;
