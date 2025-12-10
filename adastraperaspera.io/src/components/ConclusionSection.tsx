import { CheckCircle2, ExternalLink, Download, FileText } from "lucide-react";

const findings = [
  "High transcript detection efficiency (94.2%) with 2.3M total transcripts mapped",
  "24 distinct cell clusters identified through unsupervised analysis",
  "Clear spatial segregation of epithelial, immune, and stromal compartments",
  "Evidence of immune cell infiltration at epithelial-stromal interface",
  "Neighborhood analysis reveals significant cell-cell interaction patterns"
];

const ConclusionSection = () => {
  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">
        <div className="max-w-4xl mx-auto">
          <div className="text-center mb-12">
            <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
              Key Findings
            </h2>
            <p className="text-lg text-muted-foreground">
              Summary of major discoveries from the Xenium spatial analysis
            </p>
          </div>

          <div className="bg-card rounded-2xl border border-border p-8 shadow-lg mb-12">
            <div className="space-y-4">
              {findings.map((finding, index) => (
                <div 
                  key={index} 
                  className="flex items-start gap-4 p-4 rounded-xl bg-primary/5 border border-primary/10"
                >
                  <CheckCircle2 className="w-6 h-6 text-primary flex-shrink-0 mt-0.5" />
                  <p className="text-foreground">{finding}</p>
                </div>
              ))}
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-6">
            <a 
              href="#" 
              className="group p-6 rounded-2xl bg-card border border-border hover:border-primary/30 transition-all duration-300 hover:shadow-lg"
            >
              <FileText className="w-8 h-8 text-primary mb-4 group-hover:scale-110 transition-transform" />
              <h3 className="font-semibold text-foreground mb-2">Full Report</h3>
              <p className="text-sm text-muted-foreground mb-4">
                Download the complete analysis report with detailed methodology
              </p>
              <span className="text-sm text-primary flex items-center gap-1">
                View PDF <ExternalLink className="w-3 h-3" />
              </span>
            </a>

            <a 
              href="#" 
              className="group p-6 rounded-2xl bg-card border border-border hover:border-primary/30 transition-all duration-300 hover:shadow-lg"
            >
              <Download className="w-8 h-8 text-chart-2 mb-4 group-hover:scale-110 transition-transform" />
              <h3 className="font-semibold text-foreground mb-2">Raw Data</h3>
              <p className="text-sm text-muted-foreground mb-4">
                Access transcript coordinates and cell segmentation files
              </p>
              <span className="text-sm text-primary flex items-center gap-1">
                Download <ExternalLink className="w-3 h-3" />
              </span>
            </a>

            <a 
              href="#" 
              className="group p-6 rounded-2xl bg-card border border-border hover:border-primary/30 transition-all duration-300 hover:shadow-lg"
            >
              <ExternalLink className="w-8 h-8 text-chart-3 mb-4 group-hover:scale-110 transition-transform" />
              <h3 className="font-semibold text-foreground mb-2">Interactive Viewer</h3>
              <p className="text-sm text-muted-foreground mb-4">
                Explore spatial data with the Xenium Explorer application
              </p>
              <span className="text-sm text-primary flex items-center gap-1">
                Launch <ExternalLink className="w-3 h-3" />
              </span>
            </a>
          </div>
        </div>
      </div>
    </section>
  );
};

export default ConclusionSection;
