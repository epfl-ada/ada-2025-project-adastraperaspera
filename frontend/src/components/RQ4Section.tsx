import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const RQ4Section = () => {
  return (
    <section id="rq-4" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* =====================
              Sticky left panel
             ===================== */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
              <div className="text-xs uppercase tracking-wide text-muted-foreground">
                Research Question 4
              </div>

              <h2 className="mt-2 text-xl font-bold text-foreground leading-snug">
                Inferring plaque distance from gene expression
              </h2>

              <p className="mt-4 text-sm text-muted-foreground leading-relaxed">
                Can plaque proximity be inferred from a 347-gene expression vector?
                We benchmark linear and nonlinear models and analyze where predictions fail.
              </p>

              <div className="mt-6 space-y-3 text-sm">
                <a
                  href="#rq4-residuals"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">
                    Residual diagnostics
                  </div>
                  <div className="text-muted-foreground">
                    Heteroscedasticity vs distance
                  </div>
                </a>

                <a
                  href="#rq4-truepred"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">
                    True vs predicted
                  </div>
                  <div className="text-muted-foreground">
                    Target range usage
                  </div>
                </a>

                <a
                  href="#rq4-spatial"
                  className="block rounded-xl border border-border bg-background/40 px-4 py-3 hover:bg-background/70 transition"
                >
                  <div className="font-medium text-foreground">
                    Spatial diagnostics
                  </div>
                  <div className="text-muted-foreground">
                    Residual structure in tissue
                  </div>
                </a>
              </div>
            </div>
          </aside>

          {/* =====================
              Main content
             ===================== */}
          <div className="space-y-14">
            {/* ---------------------
                Residual diagnostics
               --------------------- */}
            <div id="rq4-residuals" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                Residual diagnostics
              </h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                The residuals of the XGBoost model reveal clear heteroscedasticity:
                predictions tend to overestimate plaque distance close to plaques
                and underestimate it farther away.
              </p>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <PlotFrame
                  src={`${base}plots/residuals_diagnostics.html`}
                  title="Residual diagnostics"
                  size="lg"
                  caption="Residuals (true − predicted) vs true plaque distance."
                />
              </div>
            </div>

            {/* ---------------------
                True vs predicted
               --------------------- */}
            <div id="rq4-truepred" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                True vs predicted plaque distance
              </h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Most predictions are concentrated in the 0–100 µm range, indicating
                that the model does not fully exploit the target distance spectrum.
              </p>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <PlotFrame
                  src={`${base}plots/true_vs_predicted.html`}
                  title="True vs predicted plaque distance"
                  size="lg"
                  caption="Distribution of true and predicted plaque distances."
                />
              </div>
            </div>

            {/* ---------------------
                Spatial diagnostics
               --------------------- */}
            <div id="rq4-spatial" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                Spatial diagnostics
              </h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Residuals fall into three broad categories and show strong spatial
                alignment with plaque centroids, revealing unmodeled spatial structure.
              </p>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-md">
                <PlotFrame
                  src={`${base}plots/residuals_vs_distance.html`}
                  title="Spatial residual diagnostics"
                  size="xl"
                  caption="Bivariate residual × distance bins (left) and spatial distribution (right)."
                />
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQ4Section;
