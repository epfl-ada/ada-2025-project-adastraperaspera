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
            CARD 1 — Study design + spatial overview
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
              In this project, we analyze the Xenium dataset from 10X Genomics which contains trascriptomic data accompanied by morphology images. The data comes from sagittal brain slices of 6 mice stained with 4′,6-diamidino-2-phenylindole (DAPI) fluorescent DNA-binding nucleus dye. Three mice constitute healthy controls (wild type, no induced mutations) at 2.5, 5.7, and 13.4 months of age. The remaining mice are mutated (i.e., transgenic) at 2.5, 5.7, and 17.9 months of age.
            </p>
            <p>
              Since each brain slice comes from a different mouse, the inter-mouse variation in brain morphology is very significant. Our best attempt to align a pair of most similar mice in terms of age and disease status(transgenic at 17.9 and 5.7 months) reveals significant divergences in the brain geometry, especially the dentate gyrus. Overall, the Root Mean Square Error (RMSE) for the 8 key point pairs reached 3,390 µm, which is over 50 times larger than the median cell to plaque distance.
            </p>
          </div>

          <PlotFrame
            src={`${base}plots/wt_tg_age_grid.html`}
            title="WT vs TG spatial maps by age"
            size="lg"
            caption="Attempted alignment of Tg 5.7 months old mouse onto Tg 17.9 months old mouse. Interactive: pan and zoom enabled."
          />
        </div>

        {/* =========================
            CARD 2 — Mutation & IF plaques
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
              The induced mutation forces the murine cells to express the amyloid precursor protein (App) carrying known Alzheimer's disease familial mutations. As a result of mutations, the transgenic mice express up to 5 times more of the endogenous App. This leads to early and aggressive cerebral amyloid beta (Aβ) plaque deposition as soon as 3 months of age. The Aβ plaques are revealed with immunofluorescence (IF) staining, but only in transgenic mice at 17.9 months of age; in the figure below, the plaques appear in red.
            </p>
          </div>

          <PlotFrame
            src={`${base}plots/wt_tg_age_grid.html`}
            title="WT vs TG spatial maps by age"
            size="lg"
            caption="Spatial distribution of cells in wild-type and transgenic mice across
              2.5, 5.7, and 13+ months of age. Interactive: pan and zoom enabled."
          />
        </div>

        {/* =========================
            CARD 3 — Morphology variation & alignment
           ========================= */}
        <div className="relative p-8 rounded-2xl bg-card border border-border">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center">
              <AlignHorizontalJustifyCenter className="w-5 h-5 text-primary" />
            </div>
            <h3 className="text-xl font-semibold text-foreground">
              Morphological Variability & Alignment Limits
            </h3>
          </div>

          <div className="space-y-4 text-sm text-muted-foreground leading-relaxed mb-8 max-w-4xl">
            <p>
              Each brain slice originates from a different animal, leading to substantial
              inter-mouse morphological variability.
            </p>
            <p>
              The best alignment achieved yielded an RMSE of approximately 3,390 µm,
              far larger than the biologically relevant scale.
            </p>
          </div>

          <PlotFrame
            title="Tg 5.7 → Tg 17.9 alignment attempt"
            size="md"
            placeholder="Tg 5.7 → Tg 17.9 alignment attempt (plot forthcoming)"
          />
        </div>

      </div>
    </section>
  );
};

export default MicroscopySection;
