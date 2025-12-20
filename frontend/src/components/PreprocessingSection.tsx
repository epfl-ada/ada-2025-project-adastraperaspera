import PlotFrame from "@/components/PlotFrame";
import { Grid2x2Check, CircleChevronRight  } from 'lucide-react'
import ImageGridPlot from "./ImageGridFrame";

const base = import.meta.env.BASE_URL;
const bw_grid_path = `${base}images_grid_bw/full/`;
const gold_grid_path = `${base}images_grid_gold/full/`

const WEIRD_GENE_ROWS = [
  { gene: "Cxcl10", zero_frac: "1.00", weird_score: "9.35" },
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
                        <Grid2x2Check  className="w-5 h-5 text-primary" />
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
                        href="#xenium-data"
                        className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                    >
                        Xenium AD dataset
                    </a>

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
            <div id="xenium-data" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">Xenium AD dataset</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                We analyze a 10X Genomics Xenium spatial transcriptomics dataset consisting of:
              </p>
              <div className="grid gap-4 md:grid-cols-2">
                    <div className="rounded-2xl border border-border bg-card p-5">
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">6</span> mice</li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">6</span> morphology images</li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">1</span> immunofluorescence (IF) image </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">347</span> genes</li>
                        
                    </ul>
                    </div>

                    <div className="rounded-2xl border border-border bg-card p-5">
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">351,714</span> cells</li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">78,885,074</span> transcripts </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">34.7</span> GB of data</li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">1,736</span> Aβ plaques</li>
                    </ul>
                    </div>
                </div>
              <p className="text-lg text-muted-foreground leading-relaxed">
                The data comes from sagittal brain slices of 6 mice stained with DAPI, a fluorescent DNA-binding nucleus dye. Three mice are healthy controls (wild type or Wt) at <span className="font-medium text-foreground">2.5, 5.7, and 13.4 months</span>. The remaining three are transgenic (Tg) at <span className="font-medium text-foreground">2.5, 5.7, and 17.9 months</span>.
                The remaining three are transgenic (Tg) at <span className="font-medium text-foreground">2.5, 5.7, and 17.9 months</span>.
              </p>
              <p className="text-lg text-muted-foreground leading-relaxed">
                The mutations induced in Tg mice forces the expression of amyloid precursor protein (App) with known familial AD abnormalities. Transgenic mice express up to 5× more App.
                This leads to early and aggressive Aβ plaque deposition.
                The plaques in the <span className="font-medium text-foreground">17.9 month transgenic mouse</span> are revealed with the IF staining and appear in red.
              </p>
              {/*<figure className="space-y-2">
                  <img
                    src={`${base}figures/microscopy_.png`}
                    alt="Microscopy cells unified"
                    className="w-full h-auto rounded-xl border border-border"
                  />
                  <figcaption className="mt-2 text-xs text-muted-foreground text-center">
                    Morphology images of Wt and transgenic mice across ages
                  </figcaption>
                </figure>*/}
                <ImageGridPlot
                  files={{
                    "wt-2": `wt-2.webp`,
                    "wt-5": `${base}images_grid_bw/full/wt-5.webp`,
                    "wt-13": `${base}images_grid_bw/full/wt-13.webp`,
                    "tg-2": `${base}images_grid_bw/full/tg-2.webp`,
                    "tg-5": `${base}images_grid_bw/full/tg-5.webp`,
                    "tg-17": `${base}images_grid_bw/full/tg-17.webp`,

                    // ...
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
                />


            </div>

            {/* Plaque detection */}
            <div id="plaque-detection" className="space-y-4">
              <h3 className="text-2xl font-bold text-foreground">Plaque detection and coordinate alignment</h3>
              <p className="text-lg text-muted-foreground leading-relaxed">
                To obtain plaque coordinates, we hand-labeled <span className="font-medium text-foreground">11 plaque-free regions</span> and <span className="font-medium text-foreground">9 plaques</span> across a range of sizes, then trained a <span className="font-medium text-foreground">random forest classifier</span> to segment the remaining plaques in the IF image.
              </p>
              <p className="text-lg text-muted-foreground leading-relaxed">
                Next, plaque coordinates were transformed from IF-image space to morphology-image space using a <span className="font-medium text-foreground">RANSAC-based transformation</span> trained on <span className="font-medium text-foreground">26 visually aligned landmark pairs</span>. After alignment, we achieved <span className="font-medium text-foreground">RMSE = 3.2 µm</span>, which is small relative to the <span className="font-medium text-foreground">median cell-to-plaque distance (61 µm)</span>—supporting that alignment error is unlikely to dominate distance-based trends.
              </p>
              <p className="text-lg text-muted-foreground leading-relaxed">
                After alignment, we performed plaque post-processing to improve biological plausibility and robustness: we merged intersecting plaques, removed plaques outside the brain boundary, and filtered plaques below the 5th percentile in area. This produced <span className="font-medium text-foreground">1,736 Aβ plaques</span>, visualized below.
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
                For each cell, we computed distance to the nearest plaque as the <span className="font-medium text-foreground">Euclidean distance between the cell centroid and the nearest plaque boundary</span> (not the plaque centroid). This yields a direct geometric measure of proximity to plaque surfaces.
              </p>
              <p className="text-lg text-muted-foreground leading-relaxed">
                To build intuition, we overlay plaque polygons onto the morphology image and color tissue by distance to the nearest plaque, highlighting regions that are consistently near plaques versus regions that are relatively plaque-free.
              </p>

        
              <PlotFrame
                      src={`${base}plots/cell_to_plaque_distance_map.html`}
                      title="Brain regions by distance to nearest plaque"
                      size="md"
                    />

                {/*<div className="rounded-2xl border border-border bg-card p-6">*/}
                    <p className="text-lg text-muted-foreground leading-relaxed">
                    Distances are spatially heterogeneous but concentrated near plaques:
                    </p>
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" />Maximum distance from any plaque: <span className="font-medium text-foreground">457 µm</span></li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">99%</span> of cells are within <span className="font-medium text-foreground">200 µm</span> of a plaque </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> Median distance: <span className="font-medium text-foreground">61 µm</span></li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> Standard deviation: <span className="font-medium text-foreground">44.4 µm</span></li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> Distribution is <span className="font-medium text-foreground">right-skewed</span> with a long tail of cells far from plaques</li>
                    </ul>
                {/*</div>*/}

                <PlotFrame
                      src={`${base}plots/cell_to_plaque_distance_distribution.html`}
                      title="Cell-to-plaque distance distribution (linear scale)"
                      size="md"
                    />

            </div>

            {/* Gene panel composition and sparsity */}
            <div id="gene-comp" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Gene panel composition and sparsity</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                This Xenium dataset provides single-cell expression for <span className="font-medium text-foreground">347 genes</span>, but expression is sparse. In the most AD-advanced sample (the 17.9-month transgenic mouse), <span className="font-medium text-foreground">302/347 genes</span> have <span className="font-medium text-foreground">zero median transcript count</span>. Across genes, transcript counts are <span className="font-medium text-foreground">right-skewed</span>, and the fraction of cells with nonzero counts varies widely (<span className="font-medium text-foreground">0.002 to 0.989</span>), underscoring substantial gene-dependent detection and expression variability.
              </p>
              {/*<div className="rounded-2xl border border-border bg-card p-6">*/}
                    <p className="text-lg text-muted-foreground leading-relaxed">
                The gene panel includes:
                </p>
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">248</span> markers for 8 main cell types, neuronal cortical layer markers, and non-neuronal markers  </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">83</span> genes related to activated microglia and astrocytes </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">16</span> plaque-induced genes (PIGs) curated from primary literature</li>
                    </ul>
                {/*</div>*/}

              <figure className="space-y-2">
                  <img
                    src={`${base}figures/microscopy_cells_unified_no_bg.png`}
                    alt="Microscopy cells unified"
                    className="w-full h-auto rounded-xl border border-border"
                  />
                  <figcaption className="mt-2 text-xs text-muted-foreground text-center">
                    Spatial distribution of individual cells (points) across all mice, colored by a clustering-derived component to illustrate spatial organization and local heterogeneity.
                  </figcaption>
                </figure>
            
                <p className="text-lg text-muted-foreground leading-relaxed">
                    We next compare the distribution of <span className="font-medium text-foreground">log1p-transformed</span> transcript counts for the <span className="font-medium text-foreground">16 PIGs</span> against the typical distribution across all 347 genes, focusing on the <span className="font-medium text-foreground">17.9-month transgenic mouse</span> to study advanced pathology. We plot histograms (B = 50 bins) and (for PIGs) Gaussian KDE curves.
                </p>

                <PlotFrame
                    src={`${base}figures/expression_distribution.png`}
                    title="Distribution of distances from each cell centroid to the nearest plaque boundary, showing strong right skew and a long tail of plaque-distant cells."
                    size="md"
                    fit="contain"
                />
                
                <p className="text-lg text-muted-foreground leading-relaxed">
                    A consistent pattern emerges: for <span className="font-medium text-foreground">13/16 PIGs</span>, the distribution has a <span className="font-medium text-foreground">mode at zero</span>, followed by gradual density decay at higher counts. Different PIGs decay at different rates, indicating heterogeneous activation intensity and/or cell-state specificity.
                </p>
                <p className="text-lg text-muted-foreground leading-relaxed">
                To identify PIGs with particularly unusual expression structure, we combine diagnostics
                (including <span className="font-medium text-foreground">zero inflation</span>, dispersion, and shape)
                into a composite <span className="font-medium text-foreground">weirdness score</span>. Among PIGs,
                <span className="font-medium text-foreground"> 13/16 </span> show strong zero inflation. The five most unusual PIGs are:
                </p>


                <GeneWeirdnessTable title="Weird genes (QC)" />

                <p className="text-lg text-muted-foreground leading-relaxed">
                Cxcl10 and Cd74 clearly stand out: both have weirdness scores &gt; 8 and extremely high zero proportions (&gt;0.97), foreshadowing downstream modeling challenges (e.g., very small effective signal range after log transforms).
                </p>
            </div>

            {/* Cell clustering workflow (PCA → kNN → Leiden → UMAP) */}
            <div id="cell-clustering" className="space-y-4">
              <h3 className="text-xl font-semibold text-foreground">Cell clustering workflow (PCA → kNN → Leiden → UMAP)</h3>

              <p className="text-lg text-muted-foreground leading-relaxed">
                To characterize cell types and anatomical structure, we cluster cells using their 347-dimensional expression vectors:
              </p>
              <div className="rounded-2xl border border-border bg-card p-6">
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> PCA on expression space  </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> kNN graph with  <span className="font-medium text-foreground">15 nearest neighbors</span></li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">Leiden clustering</span>, yielding <span className="font-medium text-foreground">K = 19</span> clusters </li>
                        <li className="flex items-start gap-2"><CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> 2D embedding via <span className="font-medium text-foreground">UMAP</span> for visualization </li>
                    </ul>
                </div>

              <PlotFrame
              src={`${base}plots/joint_clustering_umap.html`}
              title="Joint clustering UMAP"
              size="lg"
              caption="UMAP embedding colored by Leiden clusters (interactive)."
                />
            
                <p className="text-lg text-muted-foreground leading-relaxed">
                    Overlaying clusters on tissue reveals strong correspondence with brain morphology:
                </p>

                <PlotFrame
              src={`${base}plots/joint_clustering_overlayed.html`}
              title="Spatial overlay of clusters"
              size="lg"
              caption="Spatial coordinates colored by Leiden clusters (interactive)."
                />

                {/*<div className="rounded-2xl border border-border bg-card p-6">*/}
                    <p className="text-lg text-muted-foreground leading-relaxed">
                Representative cluster interpretations include:
                </p>
                    <ul className="space-y-3 text-muted-foreground">
                        <li className="flex items-start gap-2">
                        <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">Cluster 6 (vascular cells)</span>
                        is concentrated near the outer rim, consistent with epidural space localization.
                        </li>

                        <li className="flex items-start gap-2">
                        <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">Cluster 14</span>
                        traces hippocampal formation and follows the dentate gyrus shape; inferred as <span className="font-medium text-foreground">dentate gyrus immature glutamatergic neurons</span>.
                        
                            
                        
                        </li>

                        <li className="flex items-start gap-2">
                        <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> <span className="font-medium text-foreground">Cluster 10</span>
                        dominates the amygdala/hypothalamus region; inferred as <span className="font-medium text-foreground">hypothalamic medial mammillary glutamatergic neurons</span>.
                        </li>

                        <li className="flex items-start gap-2">
                        <CircleChevronRight className="w-4 h-4 text-primary mt-1 flex-shrink-0" /> Ventricular cavities show a distinct lining cluster
                        (<span className="font-medium text-foreground">cluster 15</span>), inferred as <span className="font-medium text-foreground">hypothalamic GnRH1-expressing glutamatergic neurons</span>,
                        consistent with hypothalamic contributions to the third ventricle walls.
                        </li>

                    </ul>
                {/*</div>*/}
                
                <p className="text-lg text-muted-foreground leading-relaxed">
                    This anatomical concordance is central for later interpretation: spatial plaque proximity effects can reflect genuine plaque biology, but also the fact that plaques and cell types are unevenly distributed across brain regions.
                </p>
                
            </div>

                

            </div>
        </div>
      </div>
    </section>
  );
};

export default PreprocessingSection;
