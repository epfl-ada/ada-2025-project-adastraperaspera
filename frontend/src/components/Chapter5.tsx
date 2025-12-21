import { Activity, CircleChevronRight } from "lucide-react";
import PlotFrame from "@/components/PlotFrame";
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

const Chapter5 = () => {
  return (
    <section id="ch-5" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + content right */}
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky panel (LEFT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Activity className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 5
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Changing neighborhoods
                  </div>
                </div>
              </div>

              <div className="mt-6 space-y-3 text-sm text-muted-foreground leading-relaxed">
                <p>Cell composition shifts with distance.</p>
                <p>The tissue near plaques forms a distinct ecosystem.</p>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll to follow the transition.
              </div>
            </div>
          </aside>

          {/* Content (RIGHT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 5: Approaching the plaque
              </h2>
            </div>

            <div className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Now we watch what happens as we walk toward a plaque.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We model the relationship between <span className="font-medium text-foreground">cluster membership</span> and <span className="font-medium text-foreground">distance to the nearest plaque</span> using logistic regression. After fitting, we find statistically significant coefficients for <span className="font-medium text-foreground">14 of 19 clusters</span>, using a <span className="font-medium text-foreground">Bonferroni-adjusted p-value threshold of 0.01</span>.
              </p>

              <ul className="space-y-3 text-muted-foreground">
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

              <div className="grid gap-6 lg:grid-cols-2">
                <StatTable title="Clusters depleted with distance (negative slopes; enriched near plaques)" rows={NEG_ROWS} />
                <StatTable title="Clusters enriched with distance (positive slopes; depleted near plaques)" rows={POS_ROWS} />
              </div>


              <p className="text-lg text-muted-foreground leading-relaxed">
                At the plaque surface, the mix of cells changes. Logistic
                regression makes it quantitative.{" "}
                <span className="font-semibold text-foreground">
                  Immune cells, vascular cells, and astrocytes become more likely
                  near plaques
                </span>
                , while{" "}
                <span className="font-semibold text-foreground">
                  many neuronal types become less likely
                </span>
                .
              </p>

              <PlotFrame
                src={`${base}plots/slopes_types.html`}
                title="Plaque distance by cell type"
                size="md"
                caption="Regression slopes for significant clusters."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                This is the sick mouse’s first clear symptom. The tissue near
                plaques is not the same tissue. It is a different ecosystem.
              </p>

              <PlotFrame
                src={`${base}plots/cluster_frequency_distance_to_plaque.html`}
                title="Cluster frequency vs. distance to plaque"
                size="lg"
                caption="Distance to plaque vs % of cells in each bin (interactive)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Some changes are gradual, others are sudden. Certain clusters
                shift sharply only after specific distance thresholds.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default Chapter5;
