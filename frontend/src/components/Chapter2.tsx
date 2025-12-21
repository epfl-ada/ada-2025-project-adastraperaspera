
import PlotFrame from "@/components/PlotFrame";
import { Ruler, CircleChevronRight } from "lucide-react";
const base = import.meta.env.BASE_URL;

const Chapter2 = () => {
  return (
    <section id="ch-2" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: content left + sticky panel right */}
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">

          {/* Content (LEFT) */}
          <div className="space-y-10">
            <div className="space-y-4">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 2: Distance matters
              </h2>

              <p className="text-lg text-foreground leading-relaxed">
                A sick brain is not sick everywhere equally. The key is{" "}
                <span className="italic">proximity</span>, how far a cell is from
                the plaque surface.
              </p>
            </div>

            <div className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                So for every cell, we compute distance to the{" "}
                <span className="font-semibold text-foreground">
                  nearest plaque boundary
                </span>
                , not the plaque center. That makes near a plaque a real
                geometric statement, not a coarse bin.
              </p>

              <PlotFrame
                      src={`${base}plots/cell_to_plaque_distance_map.html`}
                      title="Brain regions by distance to nearest plaque"
                      size="md"
                    />

              <p className="text-lg text-muted-foreground leading-relaxed">
                We also make sure we are not lying to ourselves with misalignment.
                The IF-to-morphology mapping error is{" "}
                <span className="font-semibold text-foreground">3.2 µm</span>,
                which is small compared with the{" "}
                <span className="font-semibold text-foreground">
                  typical 61 µm cell-to-plaque distance
                </span>
                .
              </p>

              <PlotFrame
                      src={`${base}plots/cell_to_plaque_distance_distribution.html`}
                      title="Cell-to-plaque distance distribution (linear scale)"
                      size="md"
                    />

              <p className="text-lg text-muted-foreground leading-relaxed">
                When we paint the tissue by proximity, we see it immediately.
                Most of the brain lives in the shadow of plaques.{" "}
                <span className="font-semibold text-foreground">
                  More than 99% of cells are within 200 µm
                </span>
                .
              </p>
            </div>
          </div>

          {/* Sticky panel (RIGHT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Ruler className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 2
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Distance to plaques
                  </div>
                </div>
              </div>

              <div className="mt-6 text-sm text-muted-foreground leading-relaxed">
                Distance is treated as a continuous geometric quantity, measured
                from each cell to the nearest plaque surface.
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll to read the analysis.
              </div>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default Chapter2;
