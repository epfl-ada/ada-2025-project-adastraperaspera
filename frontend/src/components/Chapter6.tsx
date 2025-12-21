import { Shuffle, CircleChevronRight } from "lucide-react";
import PlotFrame from "@/components/PlotFrame";
const base = import.meta.env.BASE_URL;

type CorrRow = {
  idLabel: string;
  description: string;
  sign: "Positive" | "Negative";
  rho: string;
  pigs: string;
};

const CORR_ROWS: CorrRow[] = [
  {
    idLabel: '07 - “Corticothalamic, Gluta”',
    description: "Near-projecting corticothalamic layer 6b glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '08 - “Immune”',
    description: "Immune cells (microglia, macrophages, etc.)",
    sign: "Positive",
    rho: "+1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '10 - “Hypothalamic medial, Gluta”',
    description: "Hypothalamic medial mammillary glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '12 - “Cerebral LGE, GABA”',
    description: "Cerebral nuclei LGE-derived GABAergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '14 - “Dentate, Gluta”',
    description: "Dentate gyrus immature neurons (glutamatergic)",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '15 - “Hypothalamic Gnrh1, Gluta”',
    description: "Hypothalamic GnRH1-expressing glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '16 - “Olfactory bulb, Gluta”',
    description: "Olfactory bulb Cajal–Retzius glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '17 - “Medulla, GABA”',
    description: "Medulla GABAergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '03 - “Intra/Extratelencephalic, Gluta”',
    description: "Intratelencephalic–Extratelencephalic glutamatergic neurons",
    sign: "Positive",
    rho: "+1.00",
    pigs: "Nrep only",
  },
];

const DIR_ROWS = [
  {
    direction: "Negative (ρgc = −1.00)",
    n: "7",
    reps: "07, 10, 12, 14, 15, 16, 17 (Core 14)",
  },
  {
    direction: "Positive (ρgc = +1.00)",
    n: "2",
    reps: "08 (Core 14); 03 (Nrep)",
  },
];

const REG_ROWS = [
  {
    term: "Intercept (d = 0, baseline)",
    coef: "1.8328",
    p: "< 1e−300",
    factor: "exp(1.8328) − 1 ≈ 5.25",
    interp: "Baseline expected expression at plaque surface",
  },
  {
    term: "Glia, Astrocyte Ependymal",
    coef: "+0.9494",
    p: "< 1e−300",
    factor: "2.584×",
    interp: "158% higher than baseline",
  },
  {
    term: "Glia, Oligodendrocyte Lineage",
    coef: "−0.0700",
    p: "6.5e−06",
    factor: "0.932×",
    interp: "6.8% lower than baseline",
  },
  {
    term: "Immune, Microglia Macrophage",
    coef: "+0.3562",
    p: "< 1e−300",
    factor: "1.428×",
    interp: "43% higher than baseline",
  },
  {
    term: "Neuron, GABAergic",
    coef: "−0.0271",
    p: "0.068",
    factor: "0.973×",
    interp: "2.7% lower (not significant at 0.05)",
  },
  {
    term: "Neuron, Glutamatergic",
    coef: "−0.1924",
    p: "< 1e−300",
    factor: "0.825×",
    interp: "17.5% lower than baseline",
  },
  {
    term: "Distance to plaque (per µm)",
    coef: "−0.0012",
    p: "< 1e−300",
    factor: "0.9988× per µm",
    interp: "Each µm reduces expected expression by ~0.12%",
  },
];

function RegressionCoefTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">
        Regression coefficients
      </h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Term</th>
              <th className="py-2 pr-4 text-right font-medium">Coef (log)</th>
              <th className="py-2 pr-4 text-right font-medium">p-value</th>
              <th className="py-2 pr-4 text-right font-medium">
                Natural-scale factor
              </th>
              <th className="py-2 text-left font-medium">
                Interpretation vs baseline
              </th>
            </tr>
          </thead>

          <tbody>
            {REG_ROWS.map((r) => (
              <tr key={r.term} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground whitespace-nowrap">
                  {r.term}
                </td>
                <td className="py-2 pr-4 text-right tabular-nums">
                  {r.coef}
                </td>
                <td className="py-2 pr-4 text-right tabular-nums">
                  {r.p}
                </td>
                <td className="py-2 pr-4 text-right tabular-nums">
                  {r.factor}
                </td>
                <td className="py-2 text-muted-foreground min-w-[320px]">
                  {r.interp}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}


function CorrTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">
        Cell types vs PIG expression (significant Spearman correlations)
      </h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Cell type (ID & label)</th>
              <th className="py-2 pr-4 text-left font-medium">Description</th>
              <th className="py-2 pr-4 text-right font-medium">Sign</th>
              <th className="py-2 pr-4 text-right font-medium">ρgc</th>
              <th className="py-2 text-left font-medium">PIG genes involved</th>
            </tr>
          </thead>
          <tbody>
            {CORR_ROWS.map((r) => (
              <tr key={r.idLabel} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground whitespace-nowrap">{r.idLabel}</td>
                <td className="py-2 pr-4 text-muted-foreground min-w-[360px]">{r.description}</td>
                <td className="py-2 pr-4 text-right">{r.sign}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.rho}</td>
                <td className="py-2 text-foreground">{r.pigs}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function DirectionTable() {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <h4 className="text-base font-semibold text-foreground mb-4">Counts by correlation direction</h4>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Direction</th>
              <th className="py-2 pr-4 text-right font-medium"># Cell types</th>
              <th className="py-2 text-left font-medium">Representative IDs</th>
            </tr>
          </thead>
          <tbody>
            {DIR_ROWS.map((r) => (
              <tr key={r.direction} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground">{r.direction}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.n}</td>
                <td className="py-2 text-muted-foreground">{r.reps}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

const Chapter6 = () => {
  return (
    <section id="ch-6" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: content left + sticky panel right */}
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">
          {/* Content (LEFT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 6: State, composition, or both
              </h2>
            </div>

            <div className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Here is the central mystery in the sick mouse. Are plaque induced
                genes high near plaques because cells change their state, or
                because <span className="font-semibold text-foreground">different cells</span>{" "}
                are present there?
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Our heatmap of gene enrichment shows the danger. Some clusters
                naturally express some genes more than others, so composition
                alone could create apparent gradients.
              </p>

              <PlotFrame
                src={`${base}plots/expression_per_cluster.html`}
                title="Cluster-by-gene enrichment heatmap"
                size="lg"
                caption="Cluster-by-gene enrichment (z-score) heatmap illustrating marker structure and motivating composition–expression coupling analyses.."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                We test this directly by looking across distance bins. When
                immune cells rise, a core set of PIGs rises with them. When
                neuronal clusters rise, those PIGs drop. That is not a guess, it
                is a measured correlation structure.
              </p>

              <PlotFrame
                src={`${base}plots/PIG_type_spearman.html`}
                title="PIG type correlation vs. cellular type proportion"
                size="lg"
                caption="Interactive Spearman correlation heatmap (significant cells highlighted)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                We identify 9 <span className="font-medium text-foreground">cell types</span> with significant correlations to PIG expression. Notably:
              </p>
              <ul className="space-y-3 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  <span className="font-medium text-foreground">8/9</span> are strongly correlated with a <span className="font-medium text-foreground">core set of 14 PIGs</span>.
                </li>

                <li className="flex items-start gap-2">
                  <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                  The remaining cell type is strongly correlated with <span className="font-medium text-foreground">Nrep</span> alone, suggesting a complementary pattern rather than redundancy with the core PIG set.
                </li>
              </ul>



              <p className="text-lg text-muted-foreground leading-relaxed">
                But then comes the crucial turn. Even after accounting for broad
                cell type in a joint regression,{" "}
                <span className="font-semibold text-foreground">
                  distance still predicts expression
                </span>
                , illustrated with{" "}
                <span className="italic text-foreground">Apoe</span>. So the mouse
                is telling us two things at once:
              </p>

              <PlotFrame
                src={`${base}plots/interactive_PIG_by_broad_type.html`}
                title="PIG expression by type"
                size="xl"
                caption="Mean Apoe expression across distance bins stratified by broad cell type, with 95% confidence intervals (SEM-based)."
              />

              <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
                <ol className="space-y-3 text-muted-foreground">
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>
                      The cast changes near plaques
                    </span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>
                      Even within a cast member, proximity still matters
                    </span>
                  </li>
                </ol>
              </div>

              <div className="space-y-6">
                <CorrTable />
                <DirectionTable />
              </div>


            </div>
          </div>

          {/* Sticky panel (RIGHT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Shuffle className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 6
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Composition vs state
                  </div>
                </div>
              </div>

              <div className="mt-6 space-y-3 text-sm text-muted-foreground leading-relaxed">
                <p>
                  Composition can fake gradients.
                </p>
                <p>
                  Joint models can separate who is present from how they change.
                </p>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: read the two-part conclusion at the end.
              </div>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default Chapter6;
