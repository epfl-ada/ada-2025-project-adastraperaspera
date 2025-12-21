import { Shapes,CircleChevronRight} from "lucide-react";
import PlotFrame from "@/components/PlotFrame";
const base = import.meta.env.BASE_URL;

const Chapter8 = () => {
  return (
    <section id="ch-8" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: content left + sticky panel right */}
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">
          {/* Content (LEFT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 8: Not just the nearest scar, crowding and shape
              </h2>
            </div>

            <div className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                A plaque is not just a point. It has a shape, and plaques can
                cluster.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                So we test whether plaque geometry and local plaque density add
                explanatory power.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To increase explanatory power (R²), reduce heteroscedasticity, and test
                whether plaque shape contributes to local responses, we compute geometric
                properties of the nearest plaque for each cell:{" "}
                <span className="font-medium text-foreground">
                  area, perimeter, major axis length and orientation
                </span>
                .
              </p>

              <PlotFrame
                src={`${base}plots/geometric_characteristics_of_plaque.html`}
                title="Geometric characteristics of the nearest plaque"
                size="md"
                caption="Distributions of nearest-plaque geometry features (area, perimeter, major axis, orientation) summarized via box plots."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Nearest-plaque distance captures proximity to <em>one</em> plaque, but local
                pathology may also depend on plaque <em>crowding</em>. We therefore introduce
                two additional{" "}
                <span className="font-medium text-foreground">
                  multi-plaque proximity features
                </span>{" "}
                computed within a radius{" "}
                <span className="font-medium text-foreground">R = 61 µm</span> (the median
                nearest-plaque distance): the{" "}
                <span className="font-medium text-foreground">count</span> of plaques within
                this radius, and the{" "}
                <span className="font-medium text-foreground">mean distance</span> to plaques
                within the same neighborhood.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To validate that R is informative, we check the Count distribution and CDF.
                The distribution shows meaningful variation:{" "}
                <span className="font-medium text-foreground">55.2%</span> of cells have{" "}
                <span className="font-medium text-foreground">0</span> plaques within 61 µm,
                while the remainder have up to{" "}
                <span className="font-medium text-foreground">21</span> plaques within that
                radius—indicating a usable local density signal.
              </p>

              <PlotFrame
                src={`${base}plots/multi_plaque_proximity.html`}
                title="Multi-plaque proximity features"
                size="lg"
                caption="Histogram and CDF of local plaque density (Count within R = 61 µm), validating that multi-plaque proximity provides informative variation beyond nearest-plaque distance."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Overall, geometry and crowding behave like a second layer. Biologically real
                and statistically detectable, but smaller than the dominant drivers.
              </p>
              
            </div>
          </div>

          {/* Sticky panel (RIGHT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Shapes className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 8
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Geometry and crowding
                  </div>
                </div>
              </div>

              <div className="mt-6 space-y-3 text-sm text-muted-foreground leading-relaxed">
                <p>Plaques have shape, not just location.</p>
                <p>Local density adds signal, but it is not the main driver.</p>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: treat this as a second-order effect.
              </div>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default Chapter8;
