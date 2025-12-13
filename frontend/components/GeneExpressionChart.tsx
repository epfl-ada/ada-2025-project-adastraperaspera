import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell } from "recharts";

const geneData = [
  { gene: "EPCAM", expression: 8.4, color: "hsl(180, 60%, 35%)" },
  { gene: "CD3D", expression: 7.2, color: "hsl(340, 75%, 55%)" },
  { gene: "COL1A1", expression: 6.8, color: "hsl(45, 90%, 55%)" },
  { gene: "PECAM1", expression: 5.9, color: "hsl(260, 60%, 55%)" },
  { gene: "PTPRC", expression: 5.4, color: "hsl(120, 50%, 45%)" },
  { gene: "ACTA2", expression: 4.8, color: "hsl(200, 70%, 50%)" },
  { gene: "CD68", expression: 4.2, color: "hsl(15, 80%, 55%)" },
  { gene: "KRT19", expression: 3.9, color: "hsl(280, 60%, 50%)" },
];

const GeneExpressionChart = () => {
  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">
        <div className="grid lg:grid-cols-2 gap-12 items-center">
          <div>
            <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
              Top Expressed Genes
            </h2>
            <p className="text-lg text-muted-foreground mb-8">
              Marker genes showing highest expression levels across the tissue section, 
              indicating major cell type populations and tissue architecture.
            </p>
            
            <div className="space-y-4">
              <div className="flex items-center gap-3 p-4 rounded-xl bg-card border border-border">
                <div className="w-3 h-3 rounded-full bg-chart-1" />
                <div>
                  <p className="font-medium text-foreground">EPCAM</p>
                  <p className="text-sm text-muted-foreground">Epithelial cell adhesion molecule - marks epithelial cells</p>
                </div>
              </div>
              <div className="flex items-center gap-3 p-4 rounded-xl bg-card border border-border">
                <div className="w-3 h-3 rounded-full bg-chart-2" />
                <div>
                  <p className="font-medium text-foreground">CD3D</p>
                  <p className="text-sm text-muted-foreground">T-cell surface glycoprotein - marks T lymphocytes</p>
                </div>
              </div>
              <div className="flex items-center gap-3 p-4 rounded-xl bg-card border border-border">
                <div className="w-3 h-3 rounded-full bg-chart-3" />
                <div>
                  <p className="font-medium text-foreground">COL1A1</p>
                  <p className="text-sm text-muted-foreground">Collagen type I - marks fibroblasts and stroma</p>
                </div>
              </div>
            </div>
          </div>

          <div className="bg-card rounded-2xl border border-border p-6 shadow-md">
            <h3 className="text-lg font-semibold text-foreground mb-6">
              Mean Log Expression (CPM)
            </h3>
            <div className="h-[400px]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={geneData} layout="vertical" margin={{ left: 20, right: 20 }}>
                  <XAxis 
                    type="number" 
                    domain={[0, 10]}
                    tick={{ fill: 'hsl(210, 20%, 45%)', fontSize: 12 }}
                    axisLine={{ stroke: 'hsl(210, 20%, 88%)' }}
                  />
                  <YAxis 
                    type="category" 
                    dataKey="gene" 
                    tick={{ fill: 'hsl(210, 50%, 10%)', fontSize: 13, fontWeight: 500 }}
                    axisLine={false}
                    tickLine={false}
                    width={70}
                  />
                  <Tooltip
                    contentStyle={{
                      backgroundColor: 'hsl(0, 0%, 100%)',
                      border: '1px solid hsl(210, 20%, 88%)',
                      borderRadius: '12px',
                      boxShadow: '0 8px 24px -4px hsl(210, 50%, 10%, 0.12)'
                    }}
                    formatter={(value: number) => [`${value.toFixed(2)}`, 'Expression']}
                  />
                  <Bar 
                    dataKey="expression" 
                    radius={[0, 6, 6, 0]}
                  >
                    {geneData.map((entry, index) => (
                      <Cell key={`cell-${index}`} fill={entry.color} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default GeneExpressionChart;
