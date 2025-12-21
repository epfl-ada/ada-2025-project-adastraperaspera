import PlotFrame from "@/components/PlotFrame";
import { Brain, CircleChevronRight } from 'lucide-react'

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
                Because clusters exhibit strong marker-gene structure, apparent plaque-associated expression gradients could arise solely from changes in cell-type composition rather than within-cell-type regulation. To illustrate this, we examined cluster-level gene enrichment patterns. <span className="font-medium text-foreground">Cluster 8</span> (immune cells) strongly over-expresses <span className="font-medium text-foreground">Hexb</span> (z-score ≈ 4), whereas <span className="font-medium text-foreground">cluster 14</span> (dentate gyrus immature glutamatergic neurons) under-expresses <span className="font-medium text-foreground">Cst3</span> (z-score ≈ −2). These examples highlight the need to disentangle compositional effects from true distance-dependent gene regulation.
              </p>


              <PlotFrame
                src={`${base}plots/expression_per_cluster.html`}
                title="Cluster-by-gene enrichment heatmap"
                size="lg"
                caption="Cluster-by-gene enrichment (z-score) heatmap illustrating marker structure and motivating composition–expression coupling analyses.."
              />
              <p className="text-lg text-muted-foreground leading-relaxed">
                The results confirm that expression is strongly cluster-dependent. This establishes the need for analyses that separate composition from distance.
              </p>

            </div>

            <div id="rq2-spearman" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Correlation: cell-type proportions vs mean PIG expression</h4>


              <p className="text-lg text-muted-foreground leading-relaxed">
                To test whether distance-dependent PIG expression can be explained by shifts in cell-type composition, we performed a correlation analysis across plaque-distance bins. Cells were first grouped by distance to the nearest plaque, and for each bin we computed both <span className="font-medium text-foreground">cell-type proportions</span> and <span className="font-medium text-foreground">mean PIG expression</span>. We then calculated a <span className="font-medium text-foreground">Spearman rank correlation matrix</span> relating cell-type proportions to PIG expression across bins. Spearman correlation was used to accommodate potential nonlinear or threshold-like relationships between distance and composition. Statistical significance was assessed for each gene–cell-type pair with appropriate multiple-testing correction.
              </p>

              <PlotFrame
                src={`${base}plots/PIG_type_spearman.html`}
                title="PIG type correlation vs. cellular type proportion"
                size="lg"
                caption="Interactive Spearman correlation heatmap (significant cells highlighted)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                We identify <span className="font-medium text-foreground">9 cell types</span> whose proportions are significantly correlated with PIG expression across distance bins. <span className="font-medium text-foreground">8 of these cell types</span> show strong correlations with a <span className="font-medium text-foreground">core set of 14 PIGs</span>, indicating a shared plaque-associated transcriptional program. The remaining cell type shows a strong association with <span className="font-medium text-foreground">Nrep</span> alone. This pattern suggests complementary, gene-specific relationships rather than redundancy across all PIGs.
              </p>


              <div className="space-y-6">
                <CorrTable />
                <DirectionTable />
              </div>


              <p className="text-lg text-muted-foreground leading-relaxed">
                This suggests composition is a major driver of plaque-associated PIG trends.
              </p>

            </div>

            <div id="rq2-linear" className="space-y-4">
              <h4 className="text-xl font-semibold text-foreground">Joint regression (Apoe) </h4>
              <p className="text-lg text-muted-foreground leading-relaxed">
                To quantify the contribution of cell-type composition to PIG expression gradients, we fit a joint linear regression model predicting PIG expression from <span className="font-medium text-foreground">broad cell type</span> (with vascular endothelial pericytes as the baseline) and distance to the nearest plaque. The model explains <span className="font-medium text-foreground">20.3%</span> of the variance in normalized transcript counts (<span className="font-medium text-foreground">R² = 0.203</span>). The overall model fit is highly significant (<span className="font-medium text-foreground">F-test p &lt; 2.13 × 10⁻¹⁷⁴</span>), rejecting the joint null hypothesis. These results indicate that both cell identity and plaque proximity contribute substantially to PIG expression variation.
              </p>


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
