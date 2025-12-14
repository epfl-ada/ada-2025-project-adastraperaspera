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
        <div className="relative p-8 rounded-2xl bg-card border border-border max-h-[85vh] overflow-y-auto space-y-16">

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
                In this project, we analyze the Xenium dataset from 10X Genomics which contains
                transcriptomic data accompanied by morphology images. The data comes from sagittal
                brain slices of 6 mice stained with DAPI.
              </p>
              <p>
                Due to strong inter-mouse morphological variability, aligning different animals is
                challenging. The best alignment between transgenic mice at 5.7 and 17.9 months yielded
                an RMSE of 3,390 µm, far exceeding biologically relevant spatial scales.
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
                Transgenic mice express Alzheimer’s disease–associated App mutations, leading to
                early and aggressive amyloid beta (Aβ) plaque deposition.
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
                Plaques were detected using a random forest classifier trained on manually labeled
                examples, followed by RANSAC-based coordinate transformation into morphology space.
              </p>
              <p>
                The transformation achieved an RMSE of 3.2 µm, comparable to the median
                cell-to-plaque distance (61 µm).
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

        {/* =========================
            Conclusion (outside scroll)
           ========================= */}
        <div className="max-w-4xl mx-auto mt-16">
          <h3 className="text-xl font-semibold text-foreground mb-4">
            Microscopy Summary
          </h3>
          <p className="text-sm text-muted-foreground leading-relaxed">
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
