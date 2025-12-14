import { Microscope, Image, AlignHorizontalJustifyCenter } from "lucide-react";
import PlotFrame from "@/components/PlotFrame";

const base = import.meta.env.BASE_URL;

const MicroscopySection = () => {
  return (
    <section id="microscopy" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">

        {/* Header */}
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Microscopy Data
          </h2>
          <p className="text-lg text-muted-foreground max-w-3xl mx-auto">
            Interactive exploration of morphology imaging, plaque detection,
            and spatial relationships in the Xenium dataset.
          </p>
        </div>

        {/* =========================
            SINGLE MICROSCOPY CARD
           ========================= */}
        <div className="relative p-10 rounded-2xl bg-card border border-border space-y-20">

          {/* =========================
              Section 1 — Study design + alignment
             ========================= */}
          <section>
            <div className="flex items-center gap-3 mb-6">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <Microscope className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-xl font-semibold text-foreground">
                Study Design & Spatial Alignment
              </h3>
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-10 items-start">
              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
                <p>
                  We analyze sagittal brain slices from six mice in the Xenium dataset,
                  stained with DAPI. Three mice are wild-type controls (2.5, 5.7, 13.4 months),
                  while three are transgenic (2.5, 5.7, 17.9 months).
                </p>
                <p>
                  Inter-mouse morphological variation is substantial. An attempted rigid
                  alignment between the two most comparable transgenic samples (5.7 → 17.9 months)
                  reveals large residual errors, particularly in the dentate gyrus.
                </p>
                <p>
                  The alignment RMSE reaches <strong>3,390 µm</strong>, over 50× larger than
                  the median cell-to-plaque distance.
                </p>
              </div>

              <PlotFrame
                src={`${base}plots/tg5_to_tg17_alignment.html`}
                title="Attempted alignment: Tg 5.7 → Tg 17.9"
                size="md"
                caption={
                  <>
                    <span className="font-medium text-foreground">Blue:</span> Tg 17.9 reference.
                    {" "}
                    <span className="font-medium text-foreground">Red:</span> Tg 5.7 after rigid
                    alignment. Residual mismatches highlight alignment limits.
                  </>
                }
              />
            </div>
          </section>

          {/* =========================
              Section 2 — WT vs TG overview
             ========================= */}
          <section>
            <div className="flex items-center gap-3 mb-6">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <Image className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-xl font-semibold text-foreground">
                WT vs Transgenic Spatial Overview
              </h3>
            </div>

            <PlotFrame
              src={`${base}plots/wt_tg_age_grid.html`}
              title="WT vs TG spatial maps by age"
              size="lg"
              caption="Spatial distribution of cells across ages (2.5, 5.7, 13+ months).
              All six mice are visible without scrolling; pan and zoom enabled."
            />
          </section>

          {/* =========================
              Section 3 — Plaque detection
             ========================= */}
          <section>
            <div className="flex items-center gap-3 mb-6">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-xl font-semibold text-foreground">
                Plaque Detection & Coordinate Transformation
              </h3>
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-10 items-start">
              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
                <p>
                  Plaques were detected by training a random forest classifier on
                  manually labeled plaque and plaque-free regions.
                </p>
                <p>
                  Plaque coordinates were transformed from IF image space into
                  morphology space using a RANSAC-based transform trained on
                  26 anatomical landmark pairs.
                </p>
                <p>
                  This achieved an RMSE of <strong>3.2 µm</strong>, well below the
                  median cell-to-plaque distance (61 µm). After filtering, 1,736
                  Aβ plaques remained.
                </p>
              </div>

              <PlotFrame
                src={`${base}plots/plaques_detected.html`}
                title="Detected Aβ plaques"
                size="md"
                caption="Green: convex plaques. Red dashed: non-convex plaques.
                Orange dotted outlines show sampled convex hulls."
              />
            </div>
          </section>

          {/* =========================
              Section 4 — Cell-to-plaque distances
             ========================= */}
          <section>
            <div className="flex items-center gap-3 mb-6">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-xl font-semibold text-foreground">
                Cell–Plaque Distance Analysis
              </h3>
            </div>

            <div className="space-y-10">
              <PlotFrame
                src={`${base}plots/cell_to_plaque_distance_map.html`}
                title="Brain regions by distance to nearest plaque"
                size="md"
              />

              <PlotFrame
                src={`${base}plots/cell_to_plaque_distance_distribution.html`}
                title="Cell-to-plaque distance distribution"
                size="md"
                scrollable
                caption="Linear-scale histogram. Dashed lines indicate proximity thresholds
                (30 µm, 100 µm); dotted lines show distribution quantiles."
              />
            </div>
          </section>

        </div>

        {/* =========================
            Summary (outside the card)
           ========================= */}
        <div className="max-w-4xl mx-auto mt-16">
          <h3 className="text-xl font-semibold text-foreground mb-4">
            Microscopy Summary
          </h3>
          <p className="text-sm text-muted-foreground leading-relaxed">
            Over 99% of cells lie within 200 µm of the nearest plaque, with a median
            distance of 61 µm. The distribution is strongly right-skewed, indicating
            rare but biologically distinct regions far from plaques. Together, these
            analyses motivate distance-based stratification in downstream molecular
            analyses.
          </p>
        </div>

      </div>
    </section>
  );
};

export default MicroscopySection;
