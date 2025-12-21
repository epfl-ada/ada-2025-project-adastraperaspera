import { Dna } from "lucide-react";
import PlotFrame from "@/components/PlotFrame";
const base = import.meta.env.BASE_URL;

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

const Chapter3 = () => {
  return (
    <section id="ch-3" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        {/* Layout: sticky left panel + content right */}
        <div className="grid gap-10 lg:grid-cols-[320px_1fr] items-start">
          {/* Sticky panel (LEFT) */}
          <aside className="lg:sticky lg:top-24">
            <div className="rounded-2xl border border-border bg-card p-6 shadow-sm">
              <div className="flex items-start gap-3">
                <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0">
                  <Dna className="w-5 h-5 text-primary" />
                </div>
                <div>
                  <div className="text-sm font-semibold text-foreground">
                    Chapter 3
                  </div>
                  <div className="text-sm text-muted-foreground">
                    Sparse signals
                  </div>
                </div>
              </div>

              <nav className="mt-6 space-y-2 text-sm">
                <a
                  href="#ch-3-sparsity"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  347 genes, sparse voice
                </a>
                <a
                  href="#ch-3-neighborhoods"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Plaque induced neighborhoods
                </a>
                <a
                  href="#ch-3-troublemakers"
                  className="block rounded-xl px-3 py-2 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition"
                >
                  Two troublemakers
                </a>
              </nav>

              <div className="mt-6 text-xs text-muted-foreground">
                Tip: scroll or use the navigation above.
              </div>
            </div>
          </aside>

          {/* Content (RIGHT) */}
          <div className="space-y-12">
            <div className="space-y-3">
              <h2 className="text-3xl md:text-4xl font-bold text-foreground">
                Chapter 3: The mouse speaks in sparse words
              </h2>
            </div>

            <div id="ch-3-sparsity" className="space-y-4">
              <p className="text-lg text-muted-foreground leading-relaxed">
                This mouse only speaks{" "}
                <span className="font-semibold text-foreground">
                  347 measured genes
                </span>{" "}
                in this dataset, and it speaks them sparsely. In the sickest
                sample, most genes have zero median counts. Many plaque induced
                genes are off in almost all cells and light up only in
                particular neighborhoods.
              </p>
            </div>

            <PlotFrame
                    src={`${base}plots/expression_distribution.html`}
                    title="Distribution of distances from each cell centroid to the nearest plaque boundary, showing strong right skew and a long tail of plaque-distant cells."
                    size="lg"
                    fit="contain"
                />


            <div id="ch-3-neighborhoods" className="space-y-4">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Many plaque induced genes are off in almost all cells and light up
                only in particular neighborhoods.
              </p>
            </div>

            <GeneWeirdnessTable title="Weird genes (QC)" />

            <div id="ch-3-troublemakers" className="space-y-4">
              <p className="text-lg text-muted-foreground leading-relaxed">
                Two genes (<span className="italic text-foreground">Cxcl10</span>,{" "}
                <span className="italic text-foreground">Cd74</span>) are almost
                entirely silent across cells, so silent that even when they
                matter statistically, they barely move in practice. Your{" "}
                <span className="font-semibold text-foreground">
                  weirdness score
                </span>{" "}
                flags them as special troublemakers.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
};

export default Chapter3;
