import PlotFrame from "@/components/PlotFrame";
import { Grid2x2Check, CircleChevronRight  } from 'lucide-react'
import ImageGridPlot from "./ImageGridFrame";
const base = import.meta.env.BASE_URL;

const IntroductionSection2 = () => {
  return (
    <section id="introduction" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-3xl">

        {/* Title */}
        <div className="mb-10">
          <h1 className="text-4xl md:text-5xl font-bold text-foreground mb-4">
            The Story of a Sick Mouse
          </h1>
          <p className="text-lg md:text-xl text-muted-foreground italic">
            Prologue: A brain that looks normal until you know where to look
          </p>
        </div>

        {/* Narrative text */}
        <div className="space-y-6 text-lg leading-relaxed text-foreground">
          <p>
            From far away, the slice looks like a brain slice is supposed to look:
            layers, curves, familiar anatomy.
          </p>

          <p>
            But in the transgenic mouse, the one engineered to overproduce amyloid
            precursor protein, something has been accumulating for months.
          </p>

          <p className="text-muted-foreground">
            Tiny lesions that don’t announce themselves in the morphology image
            alone.
          </p>
        </div>

        <ImageGridPlot
                          files={{
                            "wt-2": `${base}images_grid_bw/full/wt-2.webp`,
                            "wt-5": `${base}images_grid_bw/full/wt-5.webp`,
                            "wt-13": `${base}images_grid_bw/full/wt-13.webp`,
                            "tg-2": `${base}images_grid_bw/full/tg-2.webp`,
                            "tg-5": `${base}images_grid_bw/full/tg-5.webp`,
                            "tg-17": `${base}images_grid_bw/full/tg-17.webp`,
                          }}
                          keyToPos={{
                            "wt-2": [0, 0],
                            "wt-5":  [0, 1],
                            "wt-13": [0, 2],
                            "tg-2": [1, 0],
                            "tg-5": [1, 1],
                            "tg-17": [1, 2],
                          }}
                          colTicks={["2", "5.7", "13+"]}
                          rowTicks={["Wild Type", "Transgenic"]}
                          caption="Morphology images of Wt and transgenic mice across ages."
                          cellSize={180}
                        />

      </div>
    </section>
  );
};

export default IntroductionSection2;
