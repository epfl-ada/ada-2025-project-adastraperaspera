import { Shapes } from "lucide-react";

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
                So you test whether plaque geometry and local plaque density add
                explanatory power. They do, but modestly.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                It is a secondary layer. Biologically real, statistically
                detectable, but smaller than the dominant signals.
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
