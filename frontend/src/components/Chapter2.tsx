
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
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 2: Proximity
              </h2>
              <p className="text-lg text-muted-foreground leading-relaxed">
                A sick brain is not sick everywhere equally. The key is{" "}
                <span className="italic">proximity</span>, how far a cell is from
                the plaque surface.
              </p>
            </div>

            <div id="ch-2-distance" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                Measure distance to the boundary
              </h3>
              <p className="text-lg text-muted-foreground leading-relaxed">
                For every cell, we compute distance to the{" "}
                <span className="font-semibold">nearest plaque boundary</span>,
                not the plaque center. That makes{" "}
                <span className="italic">near a plaque</span> a geometric
                statement, not a coarse bin.
              </p>
            </div>

            <div id="ch-2-alignment" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                Check alignment error
              </h3>
              <p className="text-lg text-muted-foreground leading-relaxed">
                We also make sure we are not lying to ourselves with
                misalignment. The IF to morphology mapping error is{" "}
                <span className="font-semibold">3.2 µm</span>, which is small
                compared with the{" "}
                <span className="font-semibold">
                  typical 61 µm cell to plaque distance
                </span>
                .
              </p>
            </div>

            <div id="ch-2-shadow" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                Most cells live near plaques
              </h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                When we paint the tissue by proximity, we see it immediately.
                Most of the brain lives in the shadow of plaques.{" "}
                <span className="font-semibold">
                  More than 99% of cells are within 200 µm
                </span>
                .
              </p>

              <div className="grid gap-4 md:grid-cols-2">
                <div className="rounded-2xl border border-border bg-card p-5">
                  <div className="text-sm text-muted-foreground">
                    Alignment RMSE
                  </div>
                  <div className="mt-1 text-2xl font-semibold text-foreground tabular-nums">
                    3.2 µm
                  </div>
                </div>

                <div className="rounded-2xl border border-border bg-card p-5">
                  <div className="text-sm text-muted-foreground">
                    Typical cell to plaque distance
                  </div>
                  <div className="mt-1 text-2xl font-semibold text-foreground tabular-nums">
                    61 µm
                  </div>
                </div>
              </div>
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
                    Proximity is the signal
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#ch-2-distance"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Distance to boundary
                </a>
                <a
                  href="#ch-2-alignment"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Alignment check
                </a>
                <a
                  href="#ch-2-shadow"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  The plaque shadow
                </a>
              </nav>

              <div className="mt-6 rounded-xl bg-muted/40 p-4">
                <div className="text-xs text-muted-foreground mb-2">
                  Key takeaways
                </div>
                <ul className="space-y-2 text-sm text-muted-foreground">
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-0.5 flex-shrink-0" />
                    <span>
                      Use <span className="font-semibold text-foreground">boundary</span>{" "}
                      distance, not center distance
                    </span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-0.5 flex-shrink-0" />
                    <span>
                      Mapping error is{" "}
                      <span className="font-semibold text-foreground">3.2 µm</span>
                    </span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-0.5 flex-shrink-0" />
                    <span>
                      <span className="font-semibold text-foreground">99%</span>{" "}
                      of cells are within{" "}
                      <span className="font-semibold text-foreground">200 µm</span>
                    </span>
                  </li>
                </ul>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default Chapter2;
