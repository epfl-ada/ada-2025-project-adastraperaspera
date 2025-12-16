import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const performanceRows = [
  { model: "LASSO", trainR2: 0.1969, testR2: 0.1955 },
  { model: "ElasticNet", trainR2: 0.1969, testR2: 0.1955 },
  { model: "Random Forest", trainR2: 0.2252, testR2: 0.1874 },
  { model: "XGBoost", trainR2: 0.4784, testR2: 0.2577 },
];

const residualClusterRows = [
  {
    cluster: 15,
    cellType: "Hypothalamic GnRH1-expressing glutamatergic neurons",
    n: 917,
    mar: 32.2295,
  },
  {
    cluster: 18,
    cellType: "Pons glutamatergic neurons",
    n: 499,
    mar: 17.8746,
  },
];

const topGenes = [
  {
    rank: 1,
    elasticNetGene: "Gfap",
    elasticNetImp: 7.738978,
    lassoGene: "Gfap",
    lassoImp: 7.864686,
    rfGene: "Gfap",
    rfImp: 0.22845,
    xgbGene: "Spag16",
    xgbImp: 0.033058,
  },
  {
    rank: 2,
    elasticNetGene: "Spag16",
    elasticNetImp: 3.780172,
    lassoGene: "Spag16",
    lassoImp: 3.801531,
    rfGene: "Spag16",
    rfImp: 0.084622,
    xgbGene: "Lyz2",
    xgbImp: 0.023342,
  },
  {
    rank: 3,
    elasticNetGene: "Slc17a6",
    elasticNetImp: 3.025042,
    lassoGene: "Slc17a6",
    lassoImp: 3.073279,
    rfGene: "Lyz2",
    rfImp: 0.066671,
    xgbGene: "Gfap",
    xgbImp: 0.016999,
  },
  {
    rank: 4,
    elasticNetGene: "Slc17a7",
    elasticNetImp: 2.953878,
    lassoGene: "Slc17a7",
    lassoImp: 3.05738,
    rfGene: "Cabp7",
    rfImp: 0.063918,
    xgbGene: "Igf2",
    xgbImp: 0.015752,
  },
  {
    rank: 5,
    elasticNetGene: "Igf2",
    elasticNetImp: 2.907111,
    lassoGene: "Igf2",
    lassoImp: 2.949014,
    rfGene: "B2m",
    rfImp: 0.034584,
    xgbGene: "Strip2",
    xgbImp: 0.013022,
  },
];

const Th = ({ children }: { children: React.ReactNode }) => (
  <th className="px-4 py-3 text-left text-xs font-semibold text-foreground">
    {children}
  </th>
);

const Td = ({ children }: { children: React.ReactNode }) => (
  <td className="px-4 py-3 text-sm text-muted-foreground align-top">
    {children}
  </td>
);

const RQ4Section = () => {
  return (
    <section id="rq-4" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky left panel */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
              <div className="text-xs uppercase tracking-wide text-muted-foreground">
                Research Question 4
              </div>
              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                How accurately can we infer plaque distance from gene expression?
              </h2>

              <div className="mt-5 space-y-3 text-sm">
                <a
                  href="#rq4-performance"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Model performance</div>
                  <div className="text-muted-foreground">Train/Test R² + barplot</div>
                </a>

                <a
                  href="#rq4-residuals"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Residual diagnostics</div>
                  <div className="text-muted-foreground">Heteroscedasticity check</div>
                </a>

                <a
                  href="#rq4-truepred"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">True vs predicted</div>
                  <div className="text-muted-foreground">KDE + Q-Q (interactive)</div>
                </a>

                <a
                  href="#rq4-spatial"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Spatial diagnostics</div>
                  <div className="text-muted-foreground">Bivariate residual×distance</div>
                </a>

                <a
                  href="#rq4-importance"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">Top gene importance</div>
                  <div className="text-muted-foreground">Heatmap + top-5 table</div>
                </a>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                All figures below use the HTML interactive exports in{" "}
                <span className="font-mono">frontend/public/plots</span>.
              </div>
            </div>
          </aside>

          {/* Main content */}
          <div className="space-y-12">
            <div className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">RQ4: Analysis</h3>
              <p className="text-lg text-muted-foreground leading-relaxed">
                We benchmark predictive models to estimate distance to plaque from a 347-gene
                expression vector. We compare LASSO, ElasticNet, Random Forest, and XGBoost with an
                80/20 train/test split. Linear models use standardized log1p counts; tree models use
                unstandardized log1p counts.
              </p>
            </div>

            {/* =========================
                Performance table + plot
               ========================= */}
            <div id="rq4-performance" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Model performance</h4>

              <div className="rounded-2xl border border-border bg-card overflow-hidden">
                <div className="px-6 py-4 border-b border-border">
                  <div className="text-sm font-semibold text-foreground">
                    Train/Test R² (reported)
                  </div>
                </div>
                <div className="overflow-x-auto">
                  <table className="min-w-full">
                    <thead className="bg-muted/40">
                      <tr>
                        <Th>Model</Th>
                        <Th>Train R²</Th>
                        <Th>Test R²</Th>
                      </tr>
                    </thead>
                    <tbody>
                      {performanceRows.map((r) => (
                        <tr key={r.model} className="border-t border-border">
                          <Td>{r.model}</Td>
                          <Td>{r.trainR2.toFixed(4)}</Td>
                          <Td>{r.testR2.toFixed(4)}</Td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>

              <PlotFrame
                src={`${base}plots/model_performance_predict_dist.html`}
                title="Model performance (Train vs Test R²)"
                size="lg"
                caption="Interactive bar plot comparing train/test R² across models."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                XGBoost shows the largest train–test gap (overfitting signature) but the best test
                performance. LASSO and ElasticNet have nearly identical train and test R², consistent
                with lower model capacity.
              </p>
            </div>

            {/* =========================
                Residual diagnostics
               ========================= */}
            <div id="rq4-residuals" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Residual diagnostics</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                The residual distribution reveals heteroscedasticity: near plaques the model tends to
                over-predict, while further away it under-predicts.
              </p>

              <PlotFrame
                src={`${base}plots/residuals_diagnostics.html`}
                title="Residual diagnostics (residual vs true distance)"
                size="lg"
                caption="Interactive residual scatter (with smoother) for the selected model."
              />
            </div>

            {/* =========================
                True vs Predicted
               ========================= */}
            <div id="rq4-truepred" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">True vs predicted distance</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Most mass is concentrated in the 0–100 µm range. The KDE/Q-Q view makes it clear the
                model does not fully utilize the target range.
              </p>

              <PlotFrame
                src={`${base}plots/true_vs_predicted.html`}
                title="True vs predicted plaque distance"
                size="lg"
                caption="Interactive KDE (and Q-Q) comparison between true and predicted distances."
              />
            </div>

            {/* =========================
                Spatial diagnostics + table
               ========================= */}
            <div id="rq4-spatial" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Spatial diagnostics</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Residuals cluster spatially and align with plaque centroids, indicating missing
                covariates and structured error.
              </p>

              <PlotFrame
                src={`${base}plots/residuals_vs_distance.html`}
                title="Spatial diagnostics (bivariate residual × distance bins)"
                size="xl"
                caption="Interactive: left legend heatmap + right spatial map; plaques overlaid as stars."
              />

              <div className="rounded-2xl border border-border bg-card overflow-hidden">
                <div className="px-6 py-4 border-b border-border">
                  <div className="text-sm font-semibold text-foreground">
                    Mean absolute residual by cluster (excerpt)
                  </div>
                </div>
                <div className="overflow-x-auto">
                  <table className="min-w-full">
                    <thead className="bg-muted/40">
                      <tr>
                        <Th>Cluster (Leiden)</Th>
                        <Th>Inferred cell type</Th>
                        <Th>n</Th>
                        <Th>Mean absolute residual (µm)</Th>
                      </tr>
                    </thead>
                    <tbody>
                      {residualClusterRows.map((r) => (
                        <tr key={r.cluster} className="border-t border-border">
                          <Td>{r.cluster}</Td>
                          <Td>{r.cellType}</Td>
                          <Td>{r.n.toLocaleString()}</Td>
                          <Td>{r.mar.toFixed(4)}</Td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>

              <div className="text-lg text-muted-foreground leading-relaxed">
                Residuals fall into three broad categories (highly negative near plaques, highly
                positive far away, and moderate in between), consistent with heteroscedasticity and
                missing covariates.
              </div>
            </div>

            {/* =========================
                Importance heatmap + top genes table
               ========================= */}
            <div id="rq4-importance" className="space-y-6">
              <h4 className="text-xl font-semibold text-foreground">Top gene importance</h4>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Top predictive genes are fairly consistent across models: Gfap and Spag16 appear
                repeatedly, and multiple models rank Gfap as a top predictor of plaque proximity.
              </p>

              <PlotFrame
                src={`${base}plots/top_gene_importance.html`}
                title="Top predictive genes across models"
                size="lg"
                caption="Interactive heatmap of normalized gene importances per model."
              />

              <div className="rounded-2xl border border-border bg-card overflow-hidden">
                <div className="px-6 py-4 border-b border-border">
                  <div className="text-sm font-semibold text-foreground">
                    Top 5 genes per model (reported)
                  </div>
                </div>

                <div className="overflow-x-auto">
                  <table className="min-w-full">
                    <thead className="bg-muted/40">
                      <tr>
                        <Th>Rank</Th>
                        <Th>ElasticNet</Th>
                        <Th>Imp.</Th>
                        <Th>LASSO</Th>
                        <Th>Imp.</Th>
                        <Th>Random Forest</Th>
                        <Th>Imp.</Th>
                        <Th>XGBoost</Th>
                        <Th>Imp.</Th>
                      </tr>
                    </thead>
                    <tbody>
                      {topGenes.map((r) => (
                        <tr key={r.rank} className="border-t border-border">
                          <Td>{r.rank}</Td>

                          <Td>{r.elasticNetGene}</Td>
                          <Td>{r.elasticNetImp.toFixed(6)}</Td>

                          <Td>{r.lassoGene}</Td>
                          <Td>{r.lassoImp.toFixed(6)}</Td>

                          <Td>{r.rfGene}</Td>
                          <Td>{r.rfImp.toFixed(6)}</Td>

                          <Td>{r.xgbGene}</Td>
                          <Td>{r.xgbImp.toFixed(6)}</Td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>

            {/* End */}
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ4Section;
