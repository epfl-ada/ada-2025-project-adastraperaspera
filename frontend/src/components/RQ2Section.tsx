import PlotFrame from "@/components/PlotFrame";
import { Brain,CircleChevronRight} from 'lucide-react'

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

const RQ2Section = () => {
  return (
    <section id="rq-2" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky panel */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Brain className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">Research Question 2</div>
                  <div className="text-xs text-muted-foreground">
                    How are the cell type composition, PIG expression, and plaque distance related?
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#rq2-spearman"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Correlation
                </a>
                <a
                  href="#rq2-linear"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Joint regression (Apoe)
                </a>
                <a
                  href="#rq2-bytype"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Expression by type
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Main analysis */}
          <div className="space-y-12">
            <div id="rq2-motivation" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">RQ2 : How are the cell type composition, PIG expression, and plaque distance related?</h3>

               <p className="text-lg text-muted-foreground leading-relaxed">
                Cluster-level marker enrichment shows that some clusters exhibit strong over- or under-expression of specific genes. For example:
              </p>
              {/*<div className="rounded-2xl border border-border bg-card p-5">*/}
                <ul className="space-y-3 text-muted-foreground">
                 <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Cluster 8 (immune) over-expresses <span className="font-medium text-foreground">Hexb</span> (z-score ~4). 
                  </li>

                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Cluster 14 (dentate gyrus immature glutamatergic) under-expresses <span className="font-medium text-foreground">Cst3</span> (z-score ~−2).
                  </li>
                </ul>
              {/*</div>*/}

              <p className="text-lg text-muted-foreground leading-relaxed">
                This reinforces that anatomical/cell-type structure and gene expression patterns are tightly coupled and motivates a central question: <span className="font-medium text-foreground">are plaque-associated PIG gradients direct effects, or are they mediated by cell-type composition shifts?</span>
              </p>

              <PlotFrame
                src={`${base}plots/expression_per_cluster.html`}
                title="Cluster-by-gene enrichment heatmap"
                size="lg"
                caption="Cluster-by-gene enrichment (z-score) heatmap illustrating marker structure and motivating composition–expression coupling analyses.."
              />

            </div>

            <div id="rq2-spearman" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Correlation: cell-type proportions vs mean PIG expression</h4>

              {/*<div className="rounded-2xl border border-border bg-card p-5">*/}
                <ul className="space-y-3 text-muted-foreground">
                  <p className="text-lg text-muted-foreground leading-relaxed">
                We test whether distance-dependent PIG expression could be explained by changing cell-type composition. Concretely:
              </p>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                     Bin cells by plaque distance. 
                  </li>

                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    For each bin, compute <span className="font-medium text-foreground">cell-type proportions</span> and <span className="font-medium text-foreground">mean PIG expression</span>. 
                  </li>

                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Compute a <span className="font-medium text-foreground">Spearman rank correlation matrix</span> between cell-type proportions and PIG expression across bins.  
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                     Use Spearman (rather than Pearson) to accommodate plausible non-linear/step-like behaviors.  
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Perform significance testing with multiple-testing correction per gene–cell-type pair.
                  </li>
                </ul>
              {/*</div>*/}
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
              {/*</div>*/}

              <div className="space-y-6">
                <CorrTable />
                <DirectionTable />
              </div>
              <p className="text-lg text-muted-foreground leading-relaxed">
                Interpreted biologically, the core PIGs are most aligned with immune enrichment and neuronal depletion, consistent with a glial activation signature that strengthens in plaque-proximal bins.
              </p>
          </div>

            <div id="rq2-linear" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Joint regression (Apoe) </h4>

                <ul className="space-y-3 text-muted-foreground">
                  <p className="text-lg text-muted-foreground leading-relaxed">
                To quantify how much cell composition explains PIG expression gradients, we fit a joint linear regression predicting PIG expression from:
                </p>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Broad cell type (with <span className="font-medium text-foreground">Vascular Endothelial Pericyte</span> as the baseline) 
                  </li>

                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    Distance to the nearest plaque
                  </li>
                </ul>

              {/*}
              <p className="text-lg text-muted-foreground leading-relaxed">
                We illustrate results for <span className="font-medium text-foreground">Apoe</span>.
              </p>*/}

              <p className="text-lg text-muted-foreground leading-relaxed">
                The model explains <span className="font-medium text-foreground">20.3% of the variance in normalized transcript counts (R² = 0.203)</span>, with a highly significant overall fit (F-test p &lt; 2.13 × 10⁻¹⁷⁴), rejecting the joint null hypothesis.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
               Coefficients:
              </p>
              
              <RegressionCoefTable />

              <p className="text-lg text-muted-foreground leading-relaxed">
               These results reinforce two earlier observations. Apoe expression is elevated in astrocytic and immune populations and reduced in neuronal populations, particularly glutamatergic neurons. Moreover, even after accounting for cell type, a significant negative distance effect persists, consistent with a true plaque-centered expression gradient
              </p>

            </div>

            <div id="rq2-bytype" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Stratified mean expression by cell type and distance</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We also visualize mean Apoe expression by distance bin and cell type, with 95% confidence intervals based on SEM.
              </p>

              <PlotFrame
                src={`${base}plots/interactive_PIG_by_broad_type.html`}
                title="PIG expression by type"
                size="xl"
                caption="Mean Apoe expression across distance bins stratified by broad cell type, with 95% confidence intervals (SEM-based)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The stratified plot corroborates the regression interpretation: astrocyte/ependymal cells show higher Apoe than vascular baseline, and glutamatergic neurons show markedly lower Apoe, with differences that remain meaningful relative to uncertainty.
              </p>

            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ2Section;
