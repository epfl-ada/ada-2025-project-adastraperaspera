import { Image, Microscope, AlignHorizontalJustifyCenter } from "lucide-react";
const base = import.meta.env.BASE_URL;

const MiscroscopySection = () => {
  return (
    <section id="microscopy" className="py-24 bg-background">
      <div className="container mx-auto px-6">
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

        {/* NOTE:
            Put your static assets in:
              public/figures/*.png
              public/plots/*.html
            Then reference them as:
              /figures/...
              /plots/...
        */}

        {/* ROW 1: Card 1 - Plot */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 items-start mb-12">
          {/* Card 1: Dataset description */}
          <div className="relative p-6 rounded-2xl bg-card border border-border">
            <div className="flex items-center gap-3 mb-4">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <Microscope className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-lg font-semibold text-foreground">Study Design</h3>
            </div>

            <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
              <p>
                In this project, we analyze the Xenium dataset from 10X Genomics which
                contains transcriptomic data accompanied by morphology images. The data
                comes from sagittal brain slices of 6 mice with 4′,6-diamidino-2-phenylindole
                (DAPI) fluorescent DNA-binding nucleus dye.
              </p>

              <p>
                Three mice constitute healthy controls (wild type, no induced mutations)
                at 2.5, 5.7, and 13.4 months of age. The remaining mice are mutated
                (i.e., transgenic) at 2.5, 5.7, and 17.9 months of age.
              </p>
            </div>
          </div>

          {/* Plot 1 placeholder */}
          <div className="relative p-6 rounded-2xl bg-muted/30 border border-dashed border-border">
            <iframe
              src={`${base}plots/wt_tg_age_grid.html`}
              className="w-full h-[720px] rounded-xl border border-border"
              loading="lazy"
              title="WT vs TG spatial maps by age"
            />
            <div className="min-h-[300px] flex items-center justify-center rounded-xl">
              <p className="mt-2 text-xs text-center text-muted-foreground">
                Spatial distribution of cells in wild-type and transgenic mice across ages
                (2.5, 5.7, and 13+ months). Interactive: pan and zoom enabled.
              </p>
            </div>
          </div>
        </div>

        {/* ROW 2: Plot - Card 2 */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 items-start mb-12">
          {/* Plot 2 placeholder */}
          <div className="relative p-6 rounded-2xl bg-muted/30 border border-dashed border-border">
            {/* ✅ INSERT PLOT HERE
                Put files in:
                  public/figures/if_plaques.png   (if you have one)
                  public/plots/if_plaques.html    (when available)

                - Static PNG example:
                  <img
                    src="/figures/if_plaques.png"
                    alt="IF staining / Aβ plaque visualization"
                    className="w-full rounded-xl border border-border"
                  />

                - Interactive HTML example:
                  <iframe
                    src="/plots/if_plaques.html"
                    className="w-full h-[420px] rounded-xl border border-border"
                    loading="lazy"
                    title="IF staining / Aβ plaque visualization"
                  />
            */}
            <div className="min-h-[300px] flex items-center justify-center rounded-xl">
              <p className="text-sm text-muted-foreground text-center">
                Plot placeholder: IF staining / Aβ plaque visualization
              </p>
            </div>
          </div>

          {/* Card 2: Mutation + plaques */}
          <div className="relative p-6 rounded-2xl bg-card border border-border">
            <div className="flex items-center gap-3 mb-4">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <Image className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-lg font-semibold text-foreground">
                Mutation & IF Plaques
              </h3>
            </div>

            <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
              <p>
                The induced mutation forces murine cells to express the amyloid precursor
                protein (App) carrying known Alzheimer&apos;s disease familial mutations.
                As a result, the transgenic mice express up to 5 times more endogenous App.
                This leads to early and aggressive cerebral amyloid beta (Aβ) plaque deposition
                as soon as 3 months of age.
              </p>

              <p>
                The Aβ plaques are revealed with immunofluorescence (IF) staining, but only
                in transgenic mice at 17.9 months of age; in microscopy figures, plaques
                appear in red.
              </p>
            </div>
          </div>
        </div>

        {/* ROW 3: Card 3 - Plot */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 items-start">
          {/* Card 3: Alignment / variation */}
          <div className="relative p-6 rounded-2xl bg-card border border-border">
            <div className="flex items-center gap-3 mb-4">
              <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
                <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
              </div>
              <h3 className="text-lg font-semibold text-foreground">
                Morphology Variation
              </h3>
            </div>

            <div className="space-y-4 text-sm text-muted-foreground leading-relaxed">
              <p>
                Since each brain slice comes from a different mouse, inter-mouse variation
                in brain morphology is very significant. Our best attempt to align a pair of
                most similar mice (transgenic at 17.9 and 5.7 months) resulted in a Root Mean
                Square Error (RMSE) of 3,390 µm.
              </p>

              <p>
                This is over 50× larger than the median cell-to-plaque distance (61 µm),
                indicating that cross-mouse alignment is not reliable at the spatial scale
                relevant for cell–plaque analyses.
              </p>
            </div>
          </div>

          {/* Plot 3 placeholder */}
          <div className="relative p-6 rounded-2xl bg-muted/30 border border-dashed border-border">
            {/* ✅ INSERT PLOT HERE
                Put files in:
                  public/figures/Tg_17_Tg_5_alignment.png
                  public/plots/Tg_17_Tg_5_alignment.html

                - Static PNG:
                  <img
                    src="/figures/Tg_17_Tg_5_alignment.png"
                    alt="Attempted alignment of Tg 5.7 months onto Tg 17.9 months"
                    className="w-full rounded-xl border border-border"
                  />
                  <p className="mt-2 text-xs text-center text-muted-foreground">
                    Attempted alignment of Tg 5.7 months old mouse onto Tg 17.9 months old mouse
                  </p>

                - Interactive HTML (iframe):
                  <iframe
                    src="/plots/Tg_17_Tg_5_alignment.html"
                    className="w-full h-[420px] rounded-xl border border-border"
                    loading="lazy"
                    title="Tg 5.7 → Tg 17.9 alignment attempt"
                  />
            */}
            <div className="min-h-[300px] flex items-center justify-center rounded-xl">
              <p className="text-sm text-muted-foreground text-center">
                Plot placeholder: Tg 5.7 → Tg 17.9 alignment attempt
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default MiscroscopySection;
