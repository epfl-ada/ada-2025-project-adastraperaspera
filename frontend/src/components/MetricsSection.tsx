import { TrendingUp, Target, Layers, Zap } from "lucide-react";

const metrics = [
  {
    icon: TrendingUp,
    label: "Median Genes/Cell",
    value: "245",
    change: "+12%",
    description: "Average gene detection per cell"
  },
  {
    icon: Target,
    label: "Detection Rate",
    value: "94.2%",
    change: "+3.5%",
    description: "Transcript assignment accuracy"
  },
  {
    icon: Layers,
    label: "Cell Clusters",
    value: "24",
    change: "Identified",
    description: "Distinct cell populations"
  },
  {
    icon: Zap,
    label: "Spatial Resolution",
    value: "200nm",
    change: "Subcellular",
    description: "Transcript localization precision"
  }
];

const MetricsSection = () => {
  return (
    <section id="overview" className="py-24 bg-background">
      <div className="container mx-auto px-6">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Dataset Overview
          </h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Key performance metrics from the Xenium spatial transcriptomics analysis
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
          {metrics.map((metric, index) => (
            <div
              key={metric.label}
              className="group relative p-6 rounded-2xl bg-card border border-border hover:border-primary/30 transition-all duration-300 hover:shadow-lg"
              style={{ animationDelay: `${index * 100}ms` }}
            >
              <div className="absolute inset-0 rounded-2xl bg-gradient-to-br from-primary/5 to-transparent opacity-0 group-hover:opacity-100 transition-opacity" />
              
              <div className="relative">
                <div className="w-12 h-12 rounded-xl bg-primary/10 flex items-center justify-center mb-4 group-hover:scale-110 transition-transform">
                  <metric.icon className="w-6 h-6 text-primary" />
                </div>
                
                <p className="text-sm text-muted-foreground mb-1">{metric.label}</p>
                <div className="flex items-baseline gap-2 mb-2">
                  <span className="text-3xl font-bold text-foreground">{metric.value}</span>
                  <span className="text-sm font-medium text-primary">{metric.change}</span>
                </div>
                <p className="text-xs text-muted-foreground">{metric.description}</p>
              </div>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
};

export default MetricsSection;
