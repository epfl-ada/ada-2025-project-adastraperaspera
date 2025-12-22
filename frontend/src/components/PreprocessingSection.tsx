import PlotFrame from "@/components/PlotFrame";
import { Grid2x2Check, CircleChevronRight } from 'lucide-react'
import ImageGridPlot from "./ImageGridFrame";

const base = import.meta.env.BASE_URL;
const bw_grid_path = `${base}images_grid_bw/full/`;
const _grid_path = `${base}images_grid_gold/full/`;

const WEIRD_GENE_ROWS = [
  { gene: "Cxcl10", zero_frac: "0.99", weird_score: "9.35" },
  { gene: "Cd74", zero_frac: "0.98", weird_score: "8.13" },
  { gene: "Serpina3n", zero_frac: "0.88", weird_score: "0.79" },
  { gene: "C4b", zero_frac: "0.90", weird_score: "0.19" },
  { gene: "Gfap", zero_frac: "0.69", weird_score: "0.05" },
];

function GeneWeirdnessTable({
  title = "Gene weirdness summary",
  rows = WEIRD_GENE_ROWS,
}: {
  title?: string;
  rows?: { gene: string; zero_frac: string; weird_score: string }[];
}) {
  return (
    <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
      <div className="flex items-center justify-between gap-3 mb-4">
        <h4 className="text-base font-semibold text-foreground">{title}</h4>
      </div>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-muted-foreground">
              <th className="py-2 pr-4 text-left font-medium">Gene</th>
              <th className="py-2 pr-4 text-right font-medium">Zero fraction</th>
              <th className="py-2 text-right font-medium">Weird score</th>
            </tr>
          </thead>

          <tbody>
            {rows.map((r) => (
              <tr key={r.gene} className="border-b border-border/60 last:border-b-0">
                <td className="py-2 pr-4 text-foreground">{r.gene}</td>
                <td className="py-2 pr-4 text-right tabular-nums">{r.zero_frac}</td>
                <td className="py-2 text-right tabular-nums">{r.weird_score}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}


const PreprocessingSection = () => {
  return (
    <section id="preprocessing" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + scrollable content right */}
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              {/* Header (same pattern as GeneExpression / Microscopy) */}
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Grid2x2Check className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Dataset Preprocessing pipeline
                  </div>

                </div>
              </div>

              {/* Navigation (identical hover / spacing / typography) */}
              <nav className="mt-6 space-y-2 text-sm">
                

                <a
                  href="#plaque-detection"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Plaque detection & coordinates
                </a>

                <a
                  href="#plaque-dist"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Distances
                </a>

                <a
                  href="#gene-comp"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Gene composition & sparsity
                </a>
                <a
                  href="#cell-clustering"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Cell clustering
                </a>

              </nav>

              {/* Footer hint */}
              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Analysis content */}
          <div className="space-y-12">
            {/* Intro / method */}

            {/* Plaque detection */}
            <div id="plaque-detection" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">Plaque detection and coordinate alignment</h3>
              

              <p className="text-lg text-muted-foreground leading-relaxed">
                In the oldest Tg mouse (17.9 months), plaques are directly visible via immunofluorescence staining. We treat plaques as spatial objects rather than vague regions. We annotate plaque and plaque-free areas, train a classifier to segment plaques, and then refine the plaque set to improve robustness and biological plausibility (merging overlaps, removing outside-brain detections, filtering very small plaques)
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                This yields a final set of <span className="font-medium text-foreground">1,736</span> plaques that define the "lesion map" of the sick brain.
              </p>

              <PlotFrame
                src={`${base}plots/plaque_geometries.html`}
                title="Detected Aβ plaques after transformation"
                size="md"
              />

            </div>

            {/* Distance to plaque analysis */}
            <div id="plaque-dist" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Distance to plaque</h3>

              

              <p className="text-lg text-muted-foreground leading-relaxed">
                Plaque pathology is expected to be strongest close to plaques, so the key variable is distance to plaque. For each cell, we compute distance to the nearest plaque boundary (not the plaque centroid). This makes proximity a direct geometric measure of the cell’s relationship to the plaque surface.
              </p>

              <PlotFrame
                src={`${base}plots/cell_to_plaque_distance_map.html`}
                title="Brain regions by distance to nearest plaque"
                size="md"
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                A critical technical step is coordinate alignment. Plaque coordinates are transformed from IF space to morphology space using a RANSAC-based transformation trained on <span className="font-medium text-foreground">26 landmarks</span>, achieving <span className="font-medium text-foreground">RMSE = 3.2 µm</span>. This is small relative to the <span className="font-medium text-foreground">median cell–plaque distance (61 µm)</span>, supporting that distance-based trends reflect biology rather than misregistration.
              </p>

              {/*<div className="rounded-2xl border border-border bg-card p-6">*/}
              <p className="text-lg text-muted-foreground leading-relaxed">
                Plaque proximity was visualized as a continuous tissue-wide distance field, with each cell assigned its distance to the nearest plaque boundary.
                Distances are strongly concentrated near plaques, with a <span className="font-medium text-foreground">median of 61 µm </span> (SD 44.4 µm), and <span className="font-medium text-foreground">more than 99% of cells lie within 200 µm</span> of a plaque.
                The maximum observed distance is <span className="font-medium text-foreground">457 µm</span>, indicating that even the most plaque-distant regions remain relatively close in absolute terms.
                The distribution is <span className="font-medium text-foreground">right-skewed with a long tail</span>, reflecting a small subset of anatomically distinct, plaque-sparse regions.
              </p>
              

              <PlotFrame
                src={`${base}plots/cell_to_plaque_distance_distribution.html`}
                title="Cell-to-plaque distance distribution (linear scale)"
                size="md"
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                This confirms that plaque influence is widespread at the tissue scale.
              </p>

            </div>

            {/* Gene panel composition and sparsity */}
            <div id="gene-comp" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Gene panel composition and sparsity</h3>

              

              <p className="text-lg text-muted-foreground leading-relaxed">
                The Xenium panel measures <span className="font-medium text-foreground">347 genes</span>, but expression is sparse. In the most AD-advanced sample (Tg 17.9 months),
                <span className="font-medium text-foreground">302/347 genes have zero median transcript count</span>.
                Across genes, transcript counts are right-skewed and the fraction of nonzero cells varies widely.
                This means signal exists, but it is unevenly detectable and often concentrated in subsets of cells.
              </p>

              <ImageGridPlot
                files={{
                  "wt-2": `${base}images_grid_gold/full/wt-2.webp`,
                  "wt-5": `${base}images_grid_gold/full/wt-5.webp`,
                  "wt-13": `${base}images_grid_gold/full/wt-13.webp`,
                  "tg-2": `${base}images_grid_gold/full/tg-2.webp`,
                  "tg-5": `${base}images_grid_gold/full/tg-5.webp`,
                  "tg-17": `${base}images_grid_gold/full/tg-17.webp`,
                }}
                keyToPos={{
                  "wt-2": [0, 0],
                  "wt-5": [0, 1],
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

              {

              <p className="text-lg text-muted-foreground leading-relaxed">
                We focus on <span className="font-medium text-foreground">16 plaque-induced genes (PIGs)</span> curated from prior literature. Their expression distributions show a consistent pattern: <span className="font-medium text-foreground">13/16 PIGs have a mode at zero</span>, followed by long tails at higher expression. This indicates plaque-linked activation is strong in some microenvironments but absent in most cells.
              </p>

              <PlotFrame
                src={`${base}plots/expression_distribution.html`}
                title="Distribution of distances from each cell centroid to the nearest plaque boundary, showing strong right skew and a long tail of plaque-distant cells."
                size="lg"
                fit="contain"
              />

          

              <p className="text-lg text-muted-foreground leading-relaxed">
                To formalize which genes are most challenging, we compute a composite "weirdness score" combining zero inflation and distributional features.
              </p>


              <GeneWeirdnessTable title="Weird genes (QC)" />

              <p className="text-lg text-muted-foreground leading-relaxed">
                Two genes <em>Cxcl10</em> and <em>Cd74</em> are particularly sparse and show minimal dynamic range. This sparsity constrains downstream modeling and interpretation.
              </p>

            </div>

            {/* Cell clustering workflow (PCA → kNN → Leiden → UMAP) */}
            <div id="cell-clustering" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Cell clustering workflow</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Cells were clustered using expression alone (PCA → kNN → Leiden), yielding <span className="font-medium text-foreground">19 clusters</span>.:
              </p>


              <PlotFrame
                src={`${base}plots/joint_clustering_umap.html`}
                title="Joint clustering UMAP"
                size="lg"
                caption="UMAP embedding colored by Leiden clusters (interactive)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                The key observation is that clusters are not arbitrary. When projected back onto the tissue, they align with anatomical regions, indicating that gene expression encodes anatomical structure and cell identity.
              </p>

              <PlotFrame
                src={`${base}plots/joint_clustering_overlayed.html`}
                title="Spatial overlay of clusters"
                size="lg"
                caption="Spatial coordinates colored by Leiden clusters (interactive)."
              />

              <p className="text-lg text-muted-foreground leading-relaxed">
                This matters because plaques do not occur uniformly across the brain. Any plaque-distance trend can partly reflect <span className="font-medium text-foreground">which cell types and brain regions are near plaques</span>, so later analyses must consider composition and anatomy explicitly.
              </p>

              <p className="text-lg text-muted-foreground leading-relaxed">
                Several clusters align with well-defined anatomical structures, supporting biologically meaningful cell-type assignments.
                <span className="font-medium text-foreground">Cluster 6</span> is enriched along the outer tissue boundary, consistent with vascular cells localized to the epidural space.
                <span className="font-medium text-foreground">Cluster 14</span> traces the hippocampal formation and follows the dentate gyrus morphology, consistent with immature dentate gyrus glutamatergic neurons, while <span className="font-medium text-foreground">cluster 10</span> is concentrated in the amygdala and hypothalamic regions,
                consistent with medial mammillary glutamatergic neurons. In addition, <span className="font-medium text-foreground">cluster 15</span> forms a thin lining along ventricular cavities, consistent with GnRH1-expressing hypothalamic glutamatergic neurons associated with the walls of the third ventricle.
              </p>

            </div>



          </div>
        </div>
      </div>
    </section>
  );
};

export default PreprocessingSection;
