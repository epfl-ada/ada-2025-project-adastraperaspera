import { PieChart, Pie, Cell, ResponsiveContainer, Tooltip, Legend } from "recharts";

const clusterData = [
  { name: "Epithelial", value: 32, color: "hsl(180, 60%, 35%)" },
  { name: "T Cells", value: 18, color: "hsl(340, 75%, 55%)" },
  { name: "Fibroblasts", value: 15, color: "hsl(45, 90%, 55%)" },
  { name: "Macrophages", value: 12, color: "hsl(260, 60%, 55%)" },
  { name: "Endothelial", value: 10, color: "hsl(120, 50%, 45%)" },
  { name: "B Cells", value: 8, color: "hsl(200, 70%, 50%)" },
  { name: "Other", value: 5, color: "hsl(210, 20%, 70%)" },
];

const CellClusteringSection = () => {
  return (
    <section className="py-24 bg-background">
      <div className="container mx-auto px-6">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
            Cell Type Composition
          </h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Unsupervised clustering analysis reveals distinct cell populations 
            within the tissue microenvironment
          </p>
        </div>

        <div className="grid lg:grid-cols-2 gap-12 items-center">
          <div className="bg-card rounded-2xl border border-border p-8 shadow-md">
            <div className="h-[400px]">
              <ResponsiveContainer width="100%" height="100%">
                <PieChart>
                  <Pie
                    data={clusterData}
                    cx="50%"
                    cy="50%"
                    innerRadius={80}
                    outerRadius={140}
                    paddingAngle={3}
                    dataKey="value"
                    stroke="none"
                  >
                    {clusterData.map((entry, index) => (
                      <Cell key={`cell-${index}`} fill={entry.color} />
                    ))}
                  </Pie>
                  <Tooltip
                    contentStyle={{
                      backgroundColor: 'hsl(0, 0%, 100%)',
                      border: '1px solid hsl(210, 20%, 88%)',
                      borderRadius: '12px',
                      boxShadow: '0 8px 24px -4px hsl(210, 50%, 10%, 0.12)'
                    }}
                    formatter={(value: number) => [`${value}%`, 'Proportion']}
                  />
                  <Legend 
                    verticalAlign="bottom"
                    iconType="circle"
                    iconSize={10}
                    formatter={(value) => (
                      <span style={{ color: 'hsl(210, 50%, 10%)', fontSize: '13px' }}>{value}</span>
                    )}
                  />
                </PieChart>
              </ResponsiveContainer>
            </div>
          </div>

          <div className="space-y-6">
            <div className="p-6 rounded-2xl bg-card border border-border">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-4 h-4 rounded-full bg-chart-1" />
                <h3 className="text-xl font-semibold text-foreground">Epithelial Cells (32%)</h3>
              </div>
              <p className="text-muted-foreground">
                Dominant population expressing EPCAM, KRT8, KRT18. Shows distinct 
                spatial organization with clear tissue boundary delineation.
              </p>
            </div>

            <div className="p-6 rounded-2xl bg-card border border-border">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-4 h-4 rounded-full bg-chart-2" />
                <h3 className="text-xl font-semibold text-foreground">T Lymphocytes (18%)</h3>
              </div>
              <p className="text-muted-foreground">
                CD3D+/CD4+/CD8+ populations infiltrating the tissue. Higher density 
                observed in stromal regions and near epithelial interfaces.
              </p>
            </div>

            <div className="p-6 rounded-2xl bg-card border border-border">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-4 h-4 rounded-full bg-chart-3" />
                <h3 className="text-xl font-semibold text-foreground">Stromal Fibroblasts (15%)</h3>
              </div>
              <p className="text-muted-foreground">
                COL1A1+/VIM+ cells forming the supportive tissue framework. 
                Concentrated in interstitial spaces between epithelial structures.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default CellClusteringSection;
