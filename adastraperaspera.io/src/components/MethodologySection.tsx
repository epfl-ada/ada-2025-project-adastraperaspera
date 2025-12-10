import { FlaskConical, Cpu, Database, BarChart } from "lucide-react";

const steps = [
  {
    icon: FlaskConical,
    step: "01",
    title: "Sample Preparation",
    description: "Fresh-frozen tissue section (10μm) mounted on Xenium slide with optimized permeabilization protocol.",
    details: ["FFPE compatible", "10μm sections", "Quality control passed"]
  },
  {
    icon: Cpu,
    step: "02",
    title: "In Situ Hybridization",
    description: "Targeted gene panel with 480 probes hybridized for high-specificity transcript detection.",
    details: ["480 gene panel", "Padlock probes", "Rolling circle amplification"]
  },
  {
    icon: Database,
    step: "03",
    title: "Image Acquisition",
    description: "Multi-cycle fluorescent imaging at 200nm resolution with automated field-of-view scanning.",
    details: ["200nm resolution", "Z-stack imaging", "12 imaging cycles"]
  },
  {
    icon: BarChart,
    step: "04",
    title: "Data Analysis",
    description: "Cell segmentation, transcript decoding, and unsupervised clustering for cell type annotation.",
    details: ["DAPI segmentation", "Leiden clustering", "Marker-based annotation"]
  }
];

const MethodologySection = () => {
  return (
    <section className="py-24 bg-background">
      <div className="container mx-auto px-6">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Analysis Pipeline
          </h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Standardized workflow for high-quality spatial transcriptomics data generation
          </p>
        </div>

        <div className="relative">
          {/* Connection line */}
          <div className="hidden lg:block absolute top-1/2 left-0 right-0 h-0.5 bg-gradient-to-r from-transparent via-border to-transparent" />

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-8">
            {steps.map((step, index) => (
              <div key={step.step} className="relative group">
                <div className="bg-card rounded-2xl border border-border p-6 shadow-md hover:shadow-lg transition-all duration-300 hover:border-primary/30">
                  {/* Step number badge */}
                  <div className="absolute -top-3 left-6 px-3 py-1 rounded-full bg-primary text-primary-foreground text-sm font-bold">
                    {step.step}
                  </div>

                  <div className="mt-4">
                    <div className="w-14 h-14 rounded-xl bg-primary/10 flex items-center justify-center mb-4 group-hover:scale-110 transition-transform">
                      <step.icon className="w-7 h-7 text-primary" />
                    </div>

                    <h3 className="text-xl font-semibold text-foreground mb-3">{step.title}</h3>
                    <p className="text-muted-foreground text-sm mb-4">{step.description}</p>

                    <div className="space-y-2">
                      {step.details.map((detail) => (
                        <div key={detail} className="flex items-center gap-2 text-sm">
                          <div className="w-1.5 h-1.5 rounded-full bg-primary" />
                          <span className="text-foreground/80">{detail}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>
    </section>
  );
};

export default MethodologySection;
