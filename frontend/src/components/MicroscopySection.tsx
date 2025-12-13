import React from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Badge } from "@/components/ui/badge";

type Plot = {
  id: string;
  title: string;
  caption: string;
  pngSrc?: string;   // shown if htmlSrc not available
  htmlSrc?: string;  // preferred when available (interactive)
  widthClass?: string;
};

function Figure({ plot }: { plot: Plot }) {
  const { title, caption, pngSrc, htmlSrc, widthClass } = plot;

  return (
    <Card className="overflow-hidden">
      <CardHeader className="space-y-2">
        <div className="flex items-start justify-between gap-3">
          <CardTitle className="text-base md:text-lg">{title}</CardTitle>
          <Badge variant={htmlSrc ? "default" : "secondary"}>
            {htmlSrc ? "Interactive" : "Static"}
          </Badge>
        </div>
        <p className="text-sm text-muted-foreground">{caption}</p>
      </CardHeader>

      <CardContent>
        <div className="flex justify-center">
          <div className={widthClass ?? "w-full"}>
            {htmlSrc ? (
              <iframe
                title={title}
                src={htmlSrc}
                className="w-full h-[520px] rounded-xl border bg-background"
                loading="lazy"
              />
            ) : (
              <img
                src={pngSrc}
                alt={title}
                className="w-full rounded-xl border bg-background"
                loading="lazy"
              />
            )}
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

const plots: Plot[] = [
  {
    id: "six-mice",
    title: "Microscopy images of 6 mice",
    caption:
      "Sagittal brain slices (DAPI). Three wild-type controls (2.5, 5.7, 13.4 months) and three transgenic mice (2.5, 5.7, 17.9 months).",
    // Put this file in: frontend/public/figures/microscopy_6_mice.png
    pngSrc: "@/figures/microscopy_6_mice.png",
    // Later, when available, add: htmlSrc: "/plots/microscopy_6_mice.html"
    widthClass: "w-full max-w-[860px]",
  },
  {
    id: "alignment",
    title: "Inter-mouse morphology alignment (Tg 5.7 → Tg 17.9)",
    caption:
      "Best alignment attempt between two transgenic mice. RMSE = 3,390 µm — over 50× larger than the median cell-to-plaque distance (61 µm).",
    // Put this file in: frontend/public/figures/Tg_17_Tg_5_alignment.png
    pngSrc: "@/figures/Tg_17_Tg_5_alignment.png",
    // Later, when available, add: htmlSrc: "/plots/Tg_17_Tg_5_alignment.html"
    widthClass: "w-full max-w-[520px]",
  },
];

const MiscroscopySection: React.FC = () => {
  return (
    <section id="microscopy" className="py-24">
      <div className="container mx-auto px-6">
        <div className="text-center mb-14">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Microscopy Data
          </h2>
          <p className="text-lg text-muted-foreground max-w-3xl mx-auto">
            Xenium spatial transcriptomics (10x Genomics) with morphology images
            from sagittal brain slices of 6 mice stained with DAPI.
          </p>
        </div>

        {/* Narrative / Study design */}
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-10">
          <Card className="lg:col-span-2">
            <CardHeader>
              <CardTitle>Study background</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4 text-sm leading-relaxed text-muted-foreground">
              <p>
                We analyze the Xenium dataset from 10X Genomics containing
                transcriptomic measurements accompanied by morphology images.
                The data comes from sagittal brain slices of 6 mice with DAPI
                (4′,6-diamidino-2-phenylindole), a fluorescent DNA-binding
                nucleus dye.
              </p>
              <p>
                Three mice are healthy controls (wild type; no induced mutations)
                at <span className="text-foreground">2.5</span>,{" "}
                <span className="text-foreground">5.7</span>, and{" "}
                <span className="text-foreground">13.4</span> months. The
                remaining mice are transgenic at{" "}
                <span className="text-foreground">2.5</span>,{" "}
                <span className="text-foreground">5.7</span>, and{" "}
                <span className="text-foreground">17.9</span> months.
              </p>
              <p>
                The induced mutation forces cells to express the amyloid
                precursor protein (<span className="text-foreground">App</span>)
                carrying known familial Alzheimer’s mutations. Transgenic mice
                express up to 5× more endogenous App, leading to early and
                aggressive cerebral amyloid beta (Aβ) plaque deposition as soon
                as ~3 months of age.
              </p>
              <p>
                Aβ plaques are revealed with immunofluorescence (IF) staining,
                but only in transgenic mice at{" "}
                <span className="text-foreground">17.9</span> months; in the
                microscopy figure, plaques appear in red.
              </p>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Key takeaways</CardTitle>
            </CardHeader>
            <CardContent className="space-y-3 text-sm text-muted-foreground">
              <div className="flex items-start gap-2">
                <span className="mt-1 h-2 w-2 rounded-full bg-primary/70" />
                <p>
                  Strong inter-mouse morphological variation makes direct spatial
                  alignment unreliable.
                </p>
              </div>
              <div className="flex items-start gap-2">
                <span className="mt-1 h-2 w-2 rounded-full bg-primary/70" />
                <p>
                  Alignment attempt between Tg 5.7 and Tg 17.9 yields RMSE{" "}
                  <span className="text-foreground font-medium">3,390 µm</span>.
                </p>
              </div>
              <div className="flex items-start gap-2">
                <span className="mt-1 h-2 w-2 rounded-full bg-primary/70" />
                <p>
                  RMSE is &gt;50× the median cell-to-plaque distance{" "}
                  <span className="text-foreground font-medium">61 µm</span>.
                </p>
              </div>
            </CardContent>
          </Card>
        </div>

        {/* Plots */}
        <Tabs defaultValue={plots[0].id} className="w-full">
          <div className="flex justify-center mb-6">
            <TabsList className="flex flex-wrap">
              {plots.map((p) => (
                <TabsTrigger key={p.id} value={p.id}>
                  {p.title}
                </TabsTrigger>
              ))}
            </TabsList>
          </div>

          {plots.map((p) => (
            <TabsContent key={p.id} value={p.id} className="mt-0">
              <Figure plot={p} />
            </TabsContent>
          ))}
        </Tabs>
      </div>
    </section>
  );
};

export default MiscroscopySection;
