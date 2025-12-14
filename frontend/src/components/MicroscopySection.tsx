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
            CARD 1 — Study design + alignment limits (TG5 → TG17)
           ========================= */}
        <div className="relative p-8 mb-12 rounded-2xl bg-card border border-border">
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
              brain slices of 6 mice stained with 4′,6-diamidino-2-phenylindole (DAPI) fluorescent
              DNA-binding nucleus dye. Three mice constitute healthy controls (wild type, no induced
              mutations) at 2.5, 5.7, and 13.4 months of age. The remaining mice are mutated
              (i.e., transgenic) at 2.5, 5.7, and 17.9 months of age.
            </p>
            <p>
              Since each brain slice comes from a different mouse, the inter-mouse variation in brain
              morphology is very significant. Our best attempt to align a pair of most similar mice in
              terms of age and disease status (transgenic at 17.9 and 5.7 months) reveals significant
              divergences in brain geometry, especially around the dentate gyrus. Overall, the Root Mean
              Square Error (RMSE) for the 8 key point pairs reached 3,390 µm, which is over 50 times
              larger than the median cell-to-plaque distance.
            </p>
          </div>

          <PlotFrame
            src={`${base}plots/tg5_to_tg17_alignment.html`}
            title="Attempted alignment of Tg 5.7 months old mouse onto Tg 17.9 months old mouse"
            size="md"
            caption={
              <>
                <span className="font-medium text-foreground">Blue:</span> Tg 17.9 months reference tissue.{" "}
                <span className="font-medium text-foreground">Red:</span> Tg 5.7 months tissue after rigid alignment
                (flip + translation) into the Tg 17.9 coordinate frame. Large residual mismatches highlight the limits
                of cross-mouse spatial alignment.
              </>
            }
          />
        </div>

        {/* =========================
            CARD 2 — Mutation & IF plaques (WT vs TG grid)
           ========================= */}
        <div className="relative p-8 mb-12 rounded-2xl bg-card border border-border">
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
              The induced mutation forces the murine cells to express the amyloid precursor protein (App)
              carrying known Alzheimer&apos;s disease familial mutations. As a result of mutations, the transgenic
              mice express up to 5 times more of the endogenous App. This leads to early and aggressive cerebral
              amyloid beta (Aβ) plaque deposition as soon as 3 months of age. The Aβ plaques are revealed with
              immunofluorescence (IF) staining, but only in transgenic mice at 17.9 months of age; in the figure
              below, the plaques appear in red.
            </p>
          </div>

          <PlotFrame
            src={`${base}plots/wt_tg_age_grid.html`}
            title="WT vs TG spatial maps by age"
            size="lg"
            caption="Spatial distribution of cells in wild-type and transgenic mice across 2.5, 5.7, and 13+ months of age. Interactive: pan and zoom enabled."
          />
        </div>

        {/* =========================
            CARD 3 — Plaque coordinate extraction + RANSAC transform (1 plot)
           ========================= */}
        <div className="relative p-8 mb-12 rounded-2xl bg-card border border-border">
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
              To extract the coordinates of the stained plaques, we hand-labeled 11 plaque-free regions and 9 plaques
              of various sizes; we then trained a random forest classifier to mark the remaining plaques. Next, we
              transformed the plaque coordinates from the space of the IF image to that of the morphology image. For
              that, we used a RANSAC transform trained on 26 visually aligned pairs of points around important
              anatomical landmarks.
            </p>
            <p>
              After the transformation, we achieved a Root Mean Square Error (RMSE) of 3.2 µm, which compares favorably
              to the median cell-to-plaque distance at 61 µm. Following this, we merged intersecting plaques, removed
              plaques outside of the brain boundary, and filtered plaques with areas below the 5th percentile. This
              resulted in 1,736 Aβ plaques visualized in the following figure.
            </p>
          </div>

          {/* 1 plot placeholder */}
          <PlotFrame
            src={`${base}plots/plaques_detected.html`}
            title="Detected Aβ plaques after transformation and filtering"
            size="md"
            caption="Aβ plaque geometries in morphology coordinate space after RANSAC-based transformation and post-processing. Green outlines denote convex plaques, red dashed outlines denote non-convex plaques, and orange dotted outlines indicate sampled convex hulls."
          />
        </div>

        {/* =========================
            CARD 4 — Distance to nearest plaque (2 plots side-by-side)
           ========================= */}
        <div className="relative p-8 mb-12 rounded-2xl bg-card border border-border">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
              <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
            </div>
            <h3 className="text-xl font-semibold text-foreground">
              Cell-to-Plaque Distances
            </h3>
          </div>

          <div className="space-y-4 text-sm text-muted-foreground leading-relaxed mb-8 max-w-4xl">
            <p>
              Next, for each cell, we computed the distance to the nearest plaque. Namely, we calculate the Euclidean
              distance between the cell centroid and the nearest plaque boundary. In the following figure, we highlight
              which brain regions are far away from Aβ plaques and which are located nearby by overlaying plaque
              polygons onto the morphology image.
            </p>
          </div>

          {/* Two plots side by side */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 items-start">
            <PlotFrame
              title="Brain regions by distance to nearest plaque"
              size="sm"
              placeholder="Distance heatmap / map plot (to be inserted)"
              src={`${base}plots/cell_to_plaque_distance_map.html`}
            />
            <PlotFrame
              src={`${base}plots/cell_to_plaque_distance_distribution.html`}
              title="Cell-to-plaque distance distribution (linear scale)"
              size="md"
              caption="Histogram of Euclidean distances from cell centroids to the nearest plaque boundary (linear scale). Dashed lines indicate proximity thresholds (30 µm and 100 µm); dotted lines show distribution quantiles."
            />
          </div>
        </div>

        {/* =========================
            Conclusion (not in a card)
           ========================= */}
        <div className="max-w-4xl mx-auto mt-16">
          <h3 className="text-xl font-semibold text-foreground mb-4">
            Microscopy Summary
          </h3>
          <p className="text-sm text-muted-foreground leading-relaxed">
            We can see that the maximum distance from any plaque is 457 µm; however, over 99% of all cells are located
            at most 200 µm from the nearest plaque, with the median being 61 µm. The standard deviation is very
            significant at 44.4 µm. This is supported by the distance map showing that some regions are very close to
            plaques while others are relatively far. Further, the distribution of cell-to-plaque distances is
            right-skewed, with a long tail of infrequent cells that are very far from the nearest plaque.
          </p>
        </div>
      </div>
    </section>
  );
};

export default MicroscopySection;
