const DiscussionSection2 = () => {
  return (
    <section id="discussion" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-5xl space-y-12">
        {/* Title */}
        <div className="space-y-3">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground">
            Discussion and Limitations
          </h2>
        </div>

        {/* 6.a */}
        <div className="space-y-6">
          <h3 className="text-2xl font-semibold text-foreground">
            What we can conclude robustly (and what we cannot)
          </h3>

          <p className="text-lg text-muted-foreground leading-relaxed">
            Combining our analysis of cell type proportions, gene expression
            patterns, and plaque distance modeling, the results indicate that
            plaques play a central role in AD development; however, the presence
            of plaques and their impact significantly by cell type and brain
            region. Based on our findings, we can draw the following conclusions:
          </p>

          <ol className="list-decimal pl-6 space-y-4 text-muted-foreground">
            <li>
              <span className="font-semibold text-foreground">
                Plaque proximity changes cell-type composition
              </span>
              . Immune and vascular cells as well as astrocytes are enriched in
              the plaque viscinity whereas most neurons are depleted (RQ1).
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Gene expression patterns are strongly linked to cell type changes
              </span>
              . Immune enrichment and neuronal depletion explain a significant
              proportion of variation PIG expression. However, some
              distance-related effects still persist even after the cell type is
              taken into account (RQ2).
            </li>
            <li>
              <span className="font-semibold text-foreground">
                All 16 PIGs are elevated around plaques
              </span>
              . Further, <em className="text-foreground">Gfap</em> exhibits the
              strongest spatial changes, whereas{" "}
              <em className="text-foreground">Cxcl10</em> expression remains
              mostly flat due to zero inflation (RQ3).
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Gene expression somewhat encodes distance from plaque
              </span>
              . When attempting to infer plaque distance from gene expression
              alone, we get non-trivial R²; however, the residuals are strongly
              colocalized with plaques and appear heteroscedastic (RQ4).
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Coordinate-based distance models can be deceptively strong
              </span>
              . Applying the coordinate-based models to unseen mice reveals
              little cross-mouse variation; instead, the models output identical
              predictions grounded in brain anatomy as opposed to bone fide
              plaque pathology (RQ4). Even within a single mouse, spatial
              cross-validation (CV) with held-out slices of brain tissue shows a
              significant drop in R² when compared to random CV covering the
              entire brain. Some genes show the generalization gap of hundreds
              (Ctsd) to thousands (Nrep) of percent.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Per age, per genotype average PIG expression reveals AD-specific,
                age-progressive transcriptional programs
              </span>
              . Two cell types (Hypothalamic Gnrh1, Glutaergic and vascular cells)
              show the strongest AD-specific activation (RQ5).
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Plaque occurrence is brain region specific
              </span>
              . This is demonstrated by the brain segmentation inferred by the
              decision-tree. The model outputs a coarse-grained segmentation in
              homogeneous regions around the diencephalon. There, the plaques
              appear evenly spaced, encouraging the model to approximate the
              whole region with a simple average. Meanwile, the tree learned a
              fine-grained segmentation in heterogeneous regions, namely
              hippocampus and isocortex. This indicates a significant variability
              in plaque occurrence in these regions. The denser
              "super-plaque-clusters" in these heteregoneous parts of the brain
              can have an especially debilitating effect on AD patients.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Gene expression in neighboring cells is a strong predictor of PIG levels
              </span>
              . When comparing multiple regression models, those using even a
              single "other" PIG level averaged across neighbors achieve strong
              performance which generalizes to unseen brain slices. This
              highlights the importance of the local microenvironment which can
              dictate the transcriptional program of its cells.
            </li>
          </ol>
        </div>

        {/* Limitations */}
        <div className="space-y-6">
          <h3 className="text-2xl font-semibold text-foreground">
            Limitations
          </h3>

          <ol className="list-decimal pl-6 space-y-4 text-muted-foreground">
            <li>
              <span className="font-semibold text-foreground">
                Without additional experiements, we cannot infer causal relationships among PIG expression, plaques, and cell type composition
              </span>
              . Although we have uncovered strong associations between cell type
              prevalence, plaques, and gene expression, we cannot make causal
              claims on which of these three came first. Further, we have not
              ruled out a confounder which can drive changes in all 3 phenomena.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                The resolution of our analysis is limited by the quality of our alignment
              </span>
              . During the course of this study, we performed 2 types of
              alingment: IF/morphology and inter-mouse morphology image alingment.
              Both come with significant errors that limit the resolution at
              which we operate.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Homoscedasticity of advanced models is countered by sub-optimal R²
              </span>
              . Despite the uniform distribution of residual size at all plaque
              distances, our predictive models fail to achieve a perfect R².
              This means that much of the gene expression variability is
              distance-agnostic and can be explained by unrelated pathology
              mechanisms,
            </li>
          </ol>
        </div>

        <hr className="border-border" />

        {/* 6.b */}
        <div className="space-y-6">
          <h3 className="text-2xl font-semibold text-foreground">
            Confounding and interpretation risks
          </h3>

          <ul className="list-disc pl-6 space-y-4 text-muted-foreground">
            <li>
              <span className="font-semibold text-foreground">
                Anatomical confounding
              </span>
              . Plaque density varies by brain region. Thus, any model using
              coordinates can learn spurious correlation with brain anatomy.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Spatial leakage and evaluation dependence
              </span>
              . Random cross-validation (CV) can inflate performance because
              nearby cells share context. The tile-based spatial block
              cross-validation (391 tiles with 78 held out) provides a more
              strict estimate of generalization capability since sizeable chunks
              of the brain are hidden at train time. The observed gene-to-gene
              variability in the spatial block/random CV performance gap implies
              that we may over-estimate our current generalization capabilities.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Neighborhood features can encode location
              </span>
              . Neighborhood gene expression levels can act as proxies for local
              tissue identity and thus the location in the brain. Our ablation
              study (Section 4.c) partially addresses this. Within-tile
              permutation results in a 62.9% drop in spatial-block OOF
              performance. Further, substituting 100 nearest neighbors with 100
              farthest neighbors causes an 86.4% drop. These results support a
              real dependence on local structure. However, we have not fully
              eliminated the possibility of neighborhood gene expression covertly
              revealing coordinates.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Boundary artifacts and break-away cells
              </span>
              . The cells which appear detached from the brain and peripheral
              cells are often trivially far from plaques and can
              disproportionately influence segmentation, residual maps, and
              permutation behavior. This is especially apparent when tiles mix
              "continental" (main brain mass) and "island" (detached) cells.
            </li>
            <li>
              <span className="font-semibold text-foreground">
                Zero inflation
              </span>
              . Genes like <em className="text-foreground">Cxcl10</em> illustrate
              that statistically significant findings can fail to be practical
              due to flat near-all-zero expression patterns.
            </li>
          </ul>
        </div>

        <hr className="border-border" />

        {/* 6.c */}
        <div className="space-y-6">
          <h3 className="text-2xl font-semibold text-foreground">
            Conclusions
          </h3>

          <p className="text-lg text-muted-foreground leading-relaxed">
            Our findings reveal that plaques are accumulated in a
            brain-region-specific way and radiate a multi-modal influence into
            the surrounding tissue. Plaques re-model cell types and result in
            niches dominated by vascular and immune cells while being devoid of
            neurons. Further, plaques activate distinct transcriptional
            programs, with different genes affected at different distances.
          </p>

          <p className="text-lg text-muted-foreground leading-relaxed">
            When considering genes as AD therapy targets, it is crucial to
            ensure that plaque/gene association is specific to transgenic mice,
            is robust to spatial block cross validation, and persists across
            realistic distances from vasculature. In other words, the drug
            interventions should focus on near-plaque, non-sparse,
            cell-type-specific transcriptional programs that show a clear
            progression with age.
          </p>
        </div>
      </div>
    </section>
  );
};

export default DiscussionSection2;
