import { Droplets } from "lucide-react";
import PlotFrame from "@/components/PlotFrame";
const base = import.meta.env.BASE_URL;

const Chapter1 = () => {
  return (
    <section id="ch-1" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + scrollable content right */}
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky panel (LEFT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Droplets className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 1
                  </div>
                  <div className="text-sm text-muted-foreground">
                    The red scars
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#ch-1-scars"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  What the stain reveals
                </a>
                <a
                  href="#ch-1-objects"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Treat plaques as objects
                </a>
                <a
                  href="#ch-1-count"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  The final map
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Content (RIGHT) */}
          <div className="space-y-12">
            {/* Header */}
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 1: The red scars
              </h2>
              <p className="text-lg text-muted-foreground leading-relaxed">
                In the oldest transgenic mouse (17.9 months), the
                immunofluorescence channel shows what morphology alone does not.
                Amyloid beta plaques appear as bright red islands of pathology.
              </p>
            </div>

            {/* Section 1 */}
            <div id="ch-1-scars" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                What the stain reveals
              </h3>
              <p className="text-lg text-muted-foreground leading-relaxed">
                At a glance, the tissue still looks orderly. Layers and curves
                hold their shape. But the red signal marks damage that is already
                established, scattered across the slice in dense clusters and
                isolated points.
              </p>
            </div>

            {/* Section 2 */}
            <div id="ch-1-objects" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">
                Treat plaques as objects in space
              </h3>
              <p className="text-lg text-muted-foreground leading-relaxed">
                We do not treat plaques as vague regions. We treat them as
                geometric objects. First we label a small set by hand, then we
                train a model to find the rest. After detection, we clean the
                result to keep it biologically plausible. Overlapping detections
                are merged, artifacts outside the brain are removed, and tiny
                fragments are filtered out.
              </p>
            </div>

            {/* Section 3 */}
            <div id="ch-1-count" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">The final map</h3>

              <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
                <p className="text-lg text-foreground leading-relaxed">
                  The mouse now has a map of scars:{" "}
                  <span className="font-semibold">1,736 plaques</span>.
                </p>

                <PlotFrame
                    src={`${base}plots/plaque_geometries.html`}
                    title="Detected Aβ plaques after transformation"
                    size="md"
                  />

              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default Chapter1;
