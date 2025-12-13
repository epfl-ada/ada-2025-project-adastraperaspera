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
              In this project, we analyze the Xenium dataset from 10X Genomics which
              combines high-resolution spatial transcriptomics with morphology imaging.
              The data consists of sagittal brain slices stained with DAPI, enabling
              precise localization of individual cells.
            </p>
            <p>
              Three mice are wild-type controls at 2.5, 5.7, and 13.4 months of age.
              Three additional mice carry Alzheimer’s disease–associated mutations
              (transgenic) at matched or advanced ages.
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
              The induced mutation forces murine cells to express the amyloid precursor
              protein (App) carrying Alzheimer’s disease–associated familial mutations.
              This results in early and aggressive amyloid beta (Aβ) plaque deposition.
            </p>
            <p>
              Aβ plaques are revealed using immunofluorescence staining and appear only
              in transgenic mice at advanced age.
            </p>
          </div>

          <PlotFrame
            title="IF staining / Aβ plaque visualization"
            size="md"
            placeholder="IF staining / Aβ plaque visualization (plot forthcoming)"
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
