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
    idLabel: '07 — “Corticothalamic, Gluta”',
    description: "Near-projecting corticothalamic layer 6b glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '08 — “Immune”',
    description: "Immune cells (microglia, macrophages, etc.)",
    sign: "Positive",
    rho: "+1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '10 — “Hypothalamic medial, Gluta”',
    description: "Hypothalamic medial mammillary glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '12 — “Cerebral LGE, GABA”',
    description: "Cerebral nuclei LGE-derived GABAergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '14 — “Dentate, Gluta”',
    description: "Dentate gyrus immature neurons (glutamatergic)",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '15 — “Hypothalamic Gnrh1, Gluta”',
    description: "Hypothalamic GnRH1-expressing glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '16 — “Olfactory bulb, Gluta”',
    description: "Olfactory bulb Cajal–Retzius glutamatergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '17 — “Medulla, GABA”',
    description: "Medulla GABAergic neurons",
    sign: "Negative",
    rho: "−1.00",
    pigs: "Core 14 PIGs",
  },
  {
    idLabel: '03 — “Intra/Extratelencephalic, Gluta”',
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

function CorrTable() {
  return (
    <div id="rq-2" className="rounded-2xl border border-border bg-card p-6 shadow-sm">
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

const RQ2Section = () => {
  return (
    <section id="rq-2" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky left panel */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
              <div className="text-xs uppercase tracking-wide text-muted-foreground">Research Question 2</div>
              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                How are changes in cell type composition related to PIG expression across plaque distance?
              </h2>

              <div className="mt-5 space-y-3 text-sm">
                <a
                  href="#rq2-motivation"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Motivation</div>
                  <div className="text-muted-foreground">Do composition shifts explain PIG changes?</div>
                </a>

                <a
                  href="#rq2-spearman"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Spearman correlations</div>
                  <div className="text-muted-foreground">Cell-type proportions × PIG expression</div>
                </a>

                <a
                  href="#rq2-tables"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Interpretation tables</div>
                  <div className="text-muted-foreground">Which types correlate with which PIGs?</div>
                </a>

                <a
                  href="#rq2-linear"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Joint regression</div>
                  <div className="text-muted-foreground">Broad type + distance → expression</div>
                </a>

                <a
                  href="#rq2-bytype"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Expression by type</div>
                  <div className="text-muted-foreground">Mean ± 95% CI across distance bins</div>
                </a>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: use these links, or scroll the analysis on the right.
              </div>
            </div>
          </aside>

          {/* Main analysis */}
          <div className="space-y-12">
            <div id="rq2-motivation" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">RQ2 Analysis</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Next, we analyze cell-type proportions and PIG expression across distance. Our goal is to see if cell type
                composition can explain changes in PIG expression. This could reveal a mechanism for how amyloid beta plaques
                influence gene expression patterns in the surrounding area—potentially by selectively killing some cells while
                sparing or recruiting others.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Having observed a strong relationship between cell type proportions and distance to the nearest plaque, we test
                whether PIG expression changes could be driven by changes in cell composition.
              </p>
            </div>

            <div id="rq2-spearman" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Spearman rank correlation (PIG × cell type)</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We compute a Spearman rank correlation matrix between cell type proportions and average PIG expression, aligned
                on the same distance bins. We choose Spearman over Pearson because we suspect non-linear relationships (e.g.,
                step-like patterns from “mega-expressing” grouped cells). We also apply multiple testing correction to obtain a
                p-value for each gene–cell-type pair.
              </p>

              <PlotFrame
                src={`${base}plots/PIG_type_spearman.html`}
                title="PIG type correlation vs. cellular type proportion"
                size="lg"
                caption="Interactive Spearman correlation heatmap (significant cells highlighted)."
              />
            </div>

            <div id="rq2-tables" className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                We found 9 cell types showing significant Spearman rank correlations with PIG expression. Eight of these are
                strongly correlated with the core 14 PIGs. The remaining cell type is strongly correlated with a single PIG,
                Nrep.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Overall, the 14 core PIGs are overexpressed in immune cells and underexpressed in most neuronal cell types.
                Meanwhile, Nrep complements the core PIG set: it correlates with a distinct glutamatergic type that is not
                strongly tied to other PIGs.
              </p>

              <div className="space-y-6">
                <CorrTable />
                <DirectionTable />
              </div>
            </div>

            <div id="rq2-linear" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Joint regression (broad type + distance)</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We then fit a joint linear regression model to predict PIG expression from broad cell type and distance to the
                nearest plaque. We discuss one example gene, <span className="font-medium text-foreground">Apoe</span>, setting
                the baseline cell type to Vascular Endothelial Pericyte cells.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                After fitting the model, we explain 20.3% of the variance of normalized transcript count (R² = 0.203), with a
                highly significant overall F-test. The coefficients confirm earlier patterns: neuronal populations show reduced
                Apoe expression relative to the vascular reference, while astrocyte and microglia populations are enriched near
                plaques.
              </p>
            </div>

            <div id="rq2-bytype" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Mean PIG expression by distance and type (95% CI)</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We inspect mean PIG expression across distance bins and broad cell types, with 95% confidence intervals based on
                SEM. The interactive plot confirms that Glia/Astrocyte/Ependymal cells exhibit higher Apoe transcript counts than
                the vascular reference, while glutamatergic neurons show much lower expression.
              </p>

              <PlotFrame
                src={`${base}plots/interactive_PIG_by_broad_type.html`}
                title="PIG expression by type"
                size="xl"
                caption="Mean expression by distance bin and broad cell type (± 1.96×SEM)."
              />
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ2Section;
