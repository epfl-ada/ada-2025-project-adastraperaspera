import { Image, Microscope, AlignHorizontalJustifyCenter } from "lucide-react";
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
            Morphology imaging, plaque detection, and spatial context for the Xenium
            spatial transcriptomics dataset.
          </p>
        </div>

        {/* =========================
            SINGLE SCROLLABLE CARD
           ========================= */}
        <div className="relative p-10 rounded-2xl bg-card border border-border space-y-20">

          {/* =========================
              SECTION 1 — TEXT | PLOT
             ========================= */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
            <div>
              <div className="flex items-center gap-3 mb-4">
                <Microscope className="w-5 h-5 text-primary" />
                <h3 className="text-xl font-semibold text-foreground">
                  Study Design & Spatial Context
                </h3>
              </div>

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
                <p>
                  We analyze morphology images accompanying the Xenium spatial
                  transcriptomics dataset from 10X Genomics. The data consists of
                  sagittal brain slices from six mice stained with DAPI, enabling
                  precise localization of individual cells.
                </p>
                <p>
                  Three mice are wild-type controls at 2.5, 5.7, and 13.4 months.
                  Three additional mice are transgenic and develop amyloid pathology
                  at matched or advanced ages.
                </p>
                <p>
                  Large inter-mouse morphological variability complicates spatial
                  alignment across samples.
                </p>
              </div>
            </div>

            <PlotFrame
              src={`${base}plots/tg5_to_tg17_alignment.html`}
              title="Tg 5.7 → Tg 17.9 alignment attempt"
              size="md"
              caption={
                <>
                  <span className="font-medium text-foreground">Blue:</span> Tg 17.9 months
                  reference tissue.{" "}
                  <span className="font-medium text-foreground">Red:</span> Tg 5.7 months
                  tissue after rigid alignment. Large residual mismatches highlight
                  the limits of cross-mouse alignment.
                </>
              }
            />
          </div>

          {/* =========================
              SECTION 2 — PLOT | TEXT
             ========================= */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
            <PlotFrame
              src={`${base}plots/wt_tg_age_grid.html`}
              title="WT vs TG spatial maps by age"
              size="lg"
              caption="Cell spatial distributions across wild-type and transgenic mice at 2.5, 5.7, and 13+ months."
            />

            <div>
              <div className="flex items-center gap-3 mb-4">
                <Image className="w-5 h-5 text-primary" />
                <h3 className="text-xl font-semibold text-foreground">
                  Mutation & Immunofluorescence Plaques
                </h3>
              </div>

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
                <p>
                  Transgenic mice express the amyloid precursor protein (App) carrying
                  Alzheimer’s disease–associated mutations. This induces early and
                  aggressive amyloid beta (Aβ) plaque deposition.
                </p>
                <p>
                  Plaques are revealed using immunofluorescence staining and are only
                  visible in transgenic mice at advanced age.
                </p>
              </div>
            </div>
          </div>

          {/* =========================
              SECTION 3 — TEXT | PLOT
             ========================= */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
            <div>
              <div className="flex items-center gap-3 mb-4">
                <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
                <h3 className="text-xl font-semibold text-foreground">
                  Plaque Detection & Coordinate Transformation
                </h3>
              </div>

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
                <p>
                  We hand-labeled plaque-free regions and plaques to train a random
                  forest classifier, enabling automated plaque detection.
                </p>
                <p>
                  Plaque coordinates were transformed from immunofluorescence image
                  space to morphology image space using a RANSAC-based transform trained
                  on anatomical landmarks.
                </p>
                <p>
                  The transformation achieved an RMSE of 3.2 µm, well below the median
                  cell-to-plaque distance.
                </p>
              </div>
            </div>

            <PlotFrame
              src={`${base}plots/plaques_detected.html`}
              title="Detected Aβ plaques after transformation"
              size="md"
              caption="Final plaque geometries after RANSAC alignment and post-processing."
            />
          </div>

          {/* =========================
              SECTION 4 — PLOT | TEXT
             ========================= */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-start">
            <PlotFrame
              src={`${base}plots/cell_to_plaque_distance_map.html`}
              title="Brain regions by distance to nearest plaque"
              size="md"
            />

            <div>
              <div className="flex items-center gap-3 mb-4">
                <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
                <h3 className="text-xl font-semibold text-foreground">
                  Cell-to-Plaque Distances
                </h3>
              </div>

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
                <p>
                  For each cell, we compute the Euclidean distance from the cell
                  centroid to the nearest plaque boundary.
                </p>
                <p>
                  Distance maps reveal regions that are either proximal or distal
                  to amyloid pathology, providing spatial context for downstream
                  transcriptomic analyses.
                </p>
              </div>

              <div className="mt-6">
                <PlotFrame
                  src={`${base}plots/cell_to_plaque_distance_distribution.html`}
                  title="Cell-to-plaque distance distribution"
                  size="sm"
                  caption="Right-skewed distribution with median distance of 61 µm."
                />
              </div>
            </div>
          </div>

        </div>

        {/* =========================
            Conclusion (outside card)
           ========================= */}
        <div className="max-w-4xl mx-auto mt-20">
          <h3 className="text-xl font-semibold text-foreground mb-4">
            Microscopy Summary
          </h3>
          <p className="text-sm text-muted-foreground leading-relaxed">
            Over 99% of cells lie within 200 µm of the nearest plaque, with a median
            distance of 61 µm. While some brain regions remain relatively distant,
            the overall distribution is strongly right-skewed, indicating a long
            tail of infrequent, highly distal cells.
          </p>
        </div>

      </div>
    </section>
  );
};

export default MicroscopySection;
