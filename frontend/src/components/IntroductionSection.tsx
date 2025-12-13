import { Target, Brain, FlaskConical, TrendingUp } from "lucide-react";

const IntroductionSection = () => {
  return (
    <section id="introduction" className="py-24 bg-background">
      <div className="container mx-auto px-6">

        <div className="max-w-4xl mx-auto space-y-8">
          <div className="p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4 mb-6">
              <div className="w-12 h-12 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                <Target className="w-6 h-6 text-primary" />
              </div>
              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">Research Objective</h3>
                <p className="text-muted-foreground leading-relaxed">
                  This project investigates how Aβ plaques impact the surrounding microenvironment. 
                  We aim to develop a quantitative model to precisely describe the influence of the 
                  plaques on the surrounding tissue. Knowing which genes and cells are impacted at 
                  different plaque distances can refine our understanding of Alzheimer's development.
                </p>
              </div>
            </div>
          </div>

          <div className="p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4 mb-6">
              <div className="w-12 h-12 rounded-xl bg-accent/10 flex items-center justify-center flex-shrink-0">
                <FlaskConical className="w-6 h-6 text-accent" />
              </div>
              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">Drug Development Application</h3>
                <p className="text-muted-foreground leading-relaxed">
                  Developers of new drugs can use our model to select realistic targets within 
                  the plaque regions accessible from the vasculature.
                </p>
              </div>
            </div>
          </div>

          <div className="p-8 rounded-2xl bg-card border border-border">
            <div className="flex items-start gap-4">
              <div className="w-12 h-12 rounded-xl bg-chart-3/10 flex items-center justify-center flex-shrink-0">
                <Brain className="w-6 h-6 text-chart-3" />
              </div>
              <div>
                <h3 className="text-xl font-semibold text-foreground mb-3">Analysis Approach</h3>
                <p className="text-muted-foreground leading-relaxed mb-4">
                  Our story explores how plaque proximity impacts the cellular, molecular, and tissue 
                  environments. We aim to analyze the cell-to-plaque distances and apply rigorous 
                  statistical tests to describe the spatial trends in gene expression and cell composition.
                </p>
                <ul className="space-y-3 text-muted-foreground">
                  <li className="flex items-start gap-2">
                    <TrendingUp className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>Analyze interplay between gene expression and cell composition</span>
                  </li>
                  <li className="flex items-start gap-2">
                    <TrendingUp className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>Determine which phenomenon is the underlying cause</span>
                  </li>
                  <li className="flex items-start gap-2">
                    <TrendingUp className="w-4 h-4 text-primary mt-1 flex-shrink-0" />
                    <span>Benchmark predictive models to infer plaque distance from multigene expression</span>
                  </li>
                </ul>
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default IntroductionSection;
