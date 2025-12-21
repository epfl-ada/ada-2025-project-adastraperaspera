import { Flame, CircleChevronRight } from "lucide-react";

const Chapter7 = () => {
  return (
    <section id="ch-7" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + content right */}
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky panel (LEFT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Flame className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 7
                  </div>
                  <div className="text-sm text-muted-foreground">
                    The inflammatory halo (RQ3)
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#ch-7-vocabulary"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  The sickness vocabulary
                </a>
                <a
                  href="#ch-7-gradients"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  A common gradient
                </a>
                <a
                  href="#ch-7-heterogeneity"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Not all genes fade alike
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Content (RIGHT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 7: The mouse’s inflammatory halo (RQ3)
              </h2>
            </div>

            <div id="ch-7-vocabulary" className="space-y-4">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Now we focus on the{" "}
                <span className="font-semibold text-foreground">
                  16 plaque induced genes
                </span>{" "}
                as the mouse’s sickness vocabulary.
              </p>
            </div>

            <div id="ch-7-gradients" className="space-y-4">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Across distance bins, every PIG changes significantly with
                proximity. A consistent pattern appears, higher expression near
                plaques that fades outward.{" "}
                <span className="italic text-foreground">Gfap</span> is the
                clearest example, an astrocytic reactivity marker that shows the
                strongest gradient.
              </p>
            </div>

            <div id="ch-7-heterogeneity" className="space-y-4">
              <p className="text-lg text-muted-foreground leading-relaxed">
                But not all PIGs behave equally. Your{" "}
                <span className="font-semibold text-foreground">
                  distance to half expression
                </span>{" "}
                measure makes the heterogeneity tangible. Some genes drop quickly
                with distance; others, like{" "}
                <span className="italic text-foreground">Cxcl10</span>, appear
                almost flat, not because biology is absent, but because the gene
                is nearly always zero in this panel.
              </p>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
                <div className="text-sm font-semibold text-foreground mb-3">
                  What the halo tells us
                </div>
                <ul className="space-y-3 text-muted-foreground">
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>
                      Every PIG responds to proximity across bins
                    </span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>
                      <span className="italic text-foreground">Gfap</span>{" "}
                      shows the strongest fading gradient
                    </span>
                  </li>
                  <li className="flex items-start gap-2">
                    <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>
                      Flat profiles can reflect sparsity, not absence of biology
                    </span>
                  </li>
                </ul>
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default Chapter7;
