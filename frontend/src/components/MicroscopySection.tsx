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
            Context and visual overview of morphology imaging and immunofluorescence
            in the Xenium spatial transcriptomics dataset.
          </p>
        </div>

        {/* =========================
            SCROLLABLE SUPER-CARD
           ========================= */}
        <div
          className="
            relative
            rounded-2xl bg-card border border-border
            max-h-[85vh] overflow-y-auto overscroll-contain
          "
        >
          {/* Sticky scroll hint */}
          <div className="sticky top-0 z-10 bg-card/90 backdrop-blur border-b border-border px-8 py-3">
            <p className="text-xs text-muted-foreground text-center">
              Scroll inside this panel to explore microscopy analyses ↓
            </p>
          </div>

          {/* Card content */}
          <div className="p-8 space-y-16">

            {/* =========================
                SECTION 1 — Study design + alignment
               ========================= */}
            <section>
              <div className="flex items-center gap-3 mb-6">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                  <Microscope className="w-5 h-5 text-primary" />
                </div>
                <h3 className="text-xl font-semibold text-foreground">
                  Study Design & Spatial Context
                </h3>
              </div>

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed mb-8 max-w-4xl">
                <p>
                  In this project, we analyze the Xenium dataset from 10X Genomics which contains trascriptomic data accompanied by morphology images. The data comes from sagittal brain slices of 6 mice stained with 4′,6-diamidino-2-phenylindole (DAPI) fluorescent DNA-binding nucleus dye. Three mice constitute healthy controls (wild type, no induced mutations) at 2.5, 5.7, and 13.4 months of age. The remaining mice are mutated (i.e., transgenic) at 2.5, 5.7, and 17.9 months of age.
                </p>
                <p>
                  Since each brain slice comes from a different mouse, the inter-mouse variation in brain morphology is very significant. Our best attempt to align a pair of most similar mice in terms of age and disease status(transgenic at 17.9 and 5.7 months) reveals significant divergences in the brain geometry, especially the dentate gyrus. Overall, the Root Mean Square Error (RMSE) for the 8 key point pairs reached 3,390 µm, which is over 50 times larger than the median cell to plaque distance.
                </p>
              </div>

              <PlotFrame
                src={`${base}plots/tg5_to_tg17_alignment.html`}
                title="Attempted alignment of Tg 5.7 → Tg 17.9 months"
                size="md"
                caption={
                  <>
                    <span className="font-medium text-foreground">Blue:</span> Tg 17.9 months reference tissue.{" "}
                    <span className="font-medium text-foreground">Red:</span> Tg 5.7 months tissue after rigid alignment.
                  </>
                }
              />
            </section>

            <hr className="border-border" />

            {/* =========================
                SECTION 2 — WT vs TG spatial overview
               ========================= */}
            <section>
              <div className="flex items-center gap-3 mb-6">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                  <Image className="w-5 h-5 text-primary" />
                </div>
                <h3 className="text-xl font-semibold text-foreground">
                  Mutation & Immunofluorescence Plaques
                </h3>
              </div>

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed mb-8 max-w-4xl">
                <p>
                  The induced mutation forces the murine cells to express the amyloid precursor protein (App) carrying known Alzheimer's disease familial mutations. As a result of mutations, the transgenic mice express up to 5 times more of the endogenous App. This leads to early and aggressive cerebral amyloid beta (Aβ) plaque deposition as soon as 3 months of age. The Aβ plaques are revealed with immunofluorescence (IF) staining, but only in transgenic mice at 17.9 months of age; in the figure below, the plaques appear in red.
                </p>
              </div>

              <PlotFrame
                src={`${base}plots/wt_tg_age_grid.html`}
                title="WT vs TG spatial maps by age"
                size="lg"
                caption="Spatial distribution of cells in wild-type and transgenic mice across ages."
              />
            </section>

            <hr className="border-border" />

            {/* =========================
                SECTION 3 — Plaque detection & transformation
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

              <div className="space-y-4 text-sm text-muted-foreground leading-relaxed mb-8 max-w-4xl">
                <p>
                   To extract the coordinates of the stained plaques, we hand-labeled 11 plaque-free regions and 9 plaques of various sizes; we then trained a random forest classifier to mark the remaining plaques. Next, we transformed the plaque coordinates from the space of the IF image to that of the morphology image. For that, we used a RANSAC-transform trained on 26 visually aligned pairs of points around important anatomical landmarks. After the transformation, we achieved a Root Mean Square Error (RMSE) of 3.2 µm, which compares favorably to the median cell to plaque distance at 61 µm. Following this, we merged intersecting plaques, plaques outside of the brain boundary, and plaques with areas below the 5th percentile. This resulted in 1736 Aβ plaques visualized in the following figure:
                </p>
                
              </div>

              <PlotFrame
                src={`${base}plots/plaques_detected.html`}
                title="Detected Aβ plaques after transformation"
                size="md"
              />
            </section>

            <hr className="border-border" />

            {/* =========================
                SECTION 4 — Cell-to-plaque distances
               ========================= */}
            <section>
              <div className="flex items-center gap-3 mb-6">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                  <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
                </div>
                <h3 className="text-xl font-semibold text-foreground">
                  Cell-to-Plaque Distances
                </h3>

                <p>
                   Next, for each cell, we computed the distance to the nearest plaque. Namely, we calculate the Euclidean distance between the cell centroid and the nearest plaque boundary. In the following figure, we highlight which brain regions are far away from Aβ-plaques and which are located nearby. We show this by overlaying the plaque polygons onto the morphology image.
                </p>
              </div>

              <div className="space-y-10">
                <PlotFrame
                  src={`${base}plots/cell_to_plaque_distance_map.html`}
                  title="Brain regions by distance to nearest plaque"
                  size="md"
                />
                <PlotFrame
                  src={`${base}plots/cell_to_plaque_distance_distribution.html`}
                  title="Cell-to-plaque distance distribution (linear scale)"
                  size="md"
                  scrollable
                />
              </div>
            </section>

          </div>
        </div>

        {/* =========================
            Conclusion (outside scroll)
           ========================= */}
        <div className="max-w-4xl mx-auto mt-16">
          <h3 className="text-xl font-semibold text-foreground mb-4">
            Microscopy Summary
          </h3>
          <p>
            Over 99% of cells lie within 200 µm of a plaque, with a median distance of 61 µm.
            The distribution is strongly right-skewed, highlighting a small population of
            plaque-distant cells.
          </p>
        </div>

      </div>
    </section>
  );
};

export default MicroscopySection;
