import { MapPin, Layers, Grid3X3 } from "lucide-react";

const SpatialMapSection = () => {
  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Spatial Distribution
          </h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Visualization of transcript localization and cellular neighborhoods 
            across the tissue section
          </p>
        </div>

        <div className="grid lg:grid-cols-3 gap-8">
          {/* Spatial Map Placeholder */}
          <div className="lg:col-span-2 bg-card rounded-2xl border border-border overflow-hidden shadow-md">
            <div className="p-4 border-b border-border flex items-center justify-between">
              <h3 className="font-semibold text-foreground">Transcript Density Map</h3>
              <div className="flex items-center gap-2 text-sm text-muted-foreground">
                <Grid3X3 className="w-4 h-4" />
                <span>4096 × 4096 px</span>
              </div>
            </div>
            <div className="relative aspect-[4/3] bg-gradient-to-br from-background to-muted flex items-center justify-center">
              {/* Simulated spatial heatmap */}
              <div className="absolute inset-4 rounded-xl overflow-hidden">
                <div 
                  className="w-full h-full"
                  style={{
                    background: `
                      radial-gradient(ellipse 30% 40% at 30% 40%, hsl(180, 60%, 35%, 0.7), transparent),
                      radial-gradient(ellipse 25% 30% at 70% 30%, hsl(340, 75%, 55%, 0.6), transparent),
                      radial-gradient(ellipse 35% 25% at 50% 70%, hsl(45, 90%, 55%, 0.5), transparent),
                      radial-gradient(ellipse 20% 20% at 80% 60%, hsl(260, 60%, 55%, 0.5), transparent),
                      radial-gradient(ellipse 15% 15% at 20% 80%, hsl(120, 50%, 45%, 0.4), transparent),
                      linear-gradient(135deg, hsl(210, 25%, 97%) 0%, hsl(210, 30%, 92%) 100%)
                    `
                  }}
                />
                <div 
                  className="absolute inset-0 opacity-30"
                  style={{
                    backgroundImage: `
                      radial-gradient(circle at random, hsl(180, 60%, 35%, 0.3) 1px, transparent 1px)
                    `,
                    backgroundSize: '12px 12px'
                  }}
                />
              </div>
              
              {/* Annotation markers */}
              <div className="absolute top-1/4 left-1/4 flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-chart-1 animate-pulse" />
                <span className="text-xs font-medium text-foreground bg-card/90 px-2 py-1 rounded">Region A</span>
              </div>
              <div className="absolute top-1/3 right-1/3 flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-chart-2 animate-pulse" />
                <span className="text-xs font-medium text-foreground bg-card/90 px-2 py-1 rounded">Region B</span>
              </div>
              <div className="absolute bottom-1/3 left-1/2 flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-chart-3 animate-pulse" />
                <span className="text-xs font-medium text-foreground bg-card/90 px-2 py-1 rounded">Region C</span>
              </div>
            </div>
          </div>

          {/* Legend and Info */}
          <div className="space-y-6">
            <div className="bg-card rounded-2xl border border-border p-6 shadow-md">
              <h3 className="font-semibold text-foreground mb-4 flex items-center gap-2">
                <Layers className="w-5 h-5 text-primary" />
                Cell Type Legend
              </h3>
              <div className="space-y-3">
                {[
                  { color: "bg-chart-1", label: "Epithelial", density: "High" },
                  { color: "bg-chart-2", label: "Immune", density: "Medium" },
                  { color: "bg-chart-3", label: "Stromal", density: "Medium" },
                  { color: "bg-chart-4", label: "Endothelial", density: "Low" },
                  { color: "bg-chart-5", label: "Other", density: "Sparse" },
                ].map((item) => (
                  <div key={item.label} className="flex items-center justify-between">
                    <div className="flex items-center gap-3">
                      <div className={`w-3 h-3 rounded-full ${item.color}`} />
                      <span className="text-sm text-foreground">{item.label}</span>
                    </div>
                    <span className="text-xs text-muted-foreground">{item.density}</span>
                  </div>
                ))}
              </div>
            </div>

            <div className="bg-card rounded-2xl border border-border p-6 shadow-md">
              <h3 className="font-semibold text-foreground mb-4 flex items-center gap-2">
                <MapPin className="w-5 h-5 text-primary" />
                Spatial Metrics
              </h3>
              <div className="space-y-4">
                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-muted-foreground">Tissue Coverage</span>
                    <span className="font-medium text-foreground">89.3%</span>
                  </div>
                  <div className="h-2 bg-muted rounded-full overflow-hidden">
                    <div className="h-full bg-primary rounded-full" style={{ width: '89.3%' }} />
                  </div>
                </div>
                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-muted-foreground">Cell Density</span>
                    <span className="font-medium text-foreground">2,340/mm²</span>
                  </div>
                  <div className="h-2 bg-muted rounded-full overflow-hidden">
                    <div className="h-full bg-chart-2 rounded-full" style={{ width: '72%' }} />
                  </div>
                </div>
                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-muted-foreground">Neighborhood Score</span>
                    <span className="font-medium text-foreground">0.847</span>
                  </div>
                  <div className="h-2 bg-muted rounded-full overflow-hidden">
                    <div className="h-full bg-chart-3 rounded-full" style={{ width: '84.7%' }} />
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default SpatialMapSection;
