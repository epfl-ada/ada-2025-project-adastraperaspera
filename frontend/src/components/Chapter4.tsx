import { Users } from "lucide-react";
import PlotFrame from "@/components/PlotFrame";
const base = import.meta.env.BASE_URL;

const Chapter4 = () => {
  return (
    <section id="ch-4" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: content left + sticky panel right */}
        <div className="grid gap-10 lg:grid-cols-[1fr_320px] items-start">
          {/* Content (LEFT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 4: The brain’s cast of characters
              </h2>
            </div>

            <div className="space-y-6">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Before we talk about sickness, we need to know who the cells are.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We cluster cells using expression alone. And something important
                happens. The clusters fall into place along the tissue like anatomy
                remembered. Expression recovers structure.
              </p>

              <PlotFrame
              src={`${base}plots/joint_clustering_umap.html`}
              title="Joint clustering UMAP"
              size="lg"
              caption="UMAP embedding colored by Leiden clusters (interactive)."
                />

              <p className="text-lg text-muted-foreground leading-relaxed">
                This matters because plaques do not land randomly. The brain has
                neighborhoods, and so do the plaques. Any story about plaque
                proximity has to acknowledge that it is also a story about{" "}
                <span className="font-semibold text-foreground">
                  where in the brain you are
                </span>
                .
              </p>

              <PlotFrame
              src={`${base}plots/joint_clustering_overlayed.html`}
              title="Spatial overlay of clusters"
              size="lg"
              caption="Spatial coordinates colored by Leiden clusters (interactive)."
                />

            </div>
          </div>

          {/* Sticky panel (RIGHT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Users className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 4
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Cell identities
                  </div>
                </div>
              </div>

              <div className="mt-6 space-y-3 text-sm text-muted-foreground leading-relaxed">
                <p>
                  Clustering from expression alone can reconstruct anatomy.
                </p>
                <p>
                  Plaque proximity is inseparable from brain location.
                </p>
              </div>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll to read the chapter.
              </div>
            </div>
          </aside>
        </div>
      </div>
    </section>
  );
};

export default Chapter4;
