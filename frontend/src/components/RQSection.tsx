import { Brain, Microscope, BarChart3, ScanLine, Activity, MessageSquareText } from "lucide-react";

type RQCard = {
  id: string;
  title: string;
  question: string;
  href: string;
  icon: React.ElementType;
  bgImage: string; // put a path under /public (e.g. /images/rq/rq1.jpg) or `${import.meta.env.BASE_URL}...`
};

const base = import.meta.env.BASE_URL;

const RQS: RQCard[] = [
  {
    id: "rq-1",
    title: "RQ 1",
    question: "How does cell-type composition vary with plaque proximity?",
    icon: Microscope,
    bgImage: `${base}images/row-1-column-1.webp`,
  },
  {
    id: "rq-2",
    title: "RQ 2",
    question: "How are the cell type composition, PIG expression, and plaque distance related?",
    href: "#rq-2",
    icon: BarChart3,
    bgImage: `${base}images/row-1-column-2.webp`,
  },
  {
    id: "rq-3",
    title: "RQ 3",
    question: "How does the Plaque Induced Gene (PIG) expression change in plaque proximity?  ",
    href: "#rq-3",
    icon: Activity,
    bgImage: `${base}images/row-1-column-3.webp`,
  },
  {
    id: "rq-4",
    title: "RQ 4",
    question: "When modeling plaque distance, which feature modalities are most important?",
    href: "#rq-4",
    icon: ScanLine,
    bgImage: `${base}images/row-1-column-4.webp`,
  },
  {
    id: "rq-5",
    title: "RQ 5",
    question: "How does the gene expression change with age for each cell type and mouse group?",
    href: "#rq-5",
    icon: Brain,
    bgImage: `${base}images/row-1-column-5.webp`,
  },
  {
    id: "discussion",
    title: "Discussion",
    question: "Conclusions, caveats, and implications for targeting.",
    href: "#discussion",
    icon: MessageSquareText,
    bgImage: `${base}images/rq/discussion.jpg`,
  },
];

const RQSection = () => {
  return (
    <section id="research-questions" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-6xl">
        <div className="text-center mb-10">
          <h2 className="text-3xl md:text-4xl font-bold text-foreground">Research questions</h2>
        </div>

        {/* One-row presentation layout on lg+, responsive grid below */}
        <div className="grid gap-6 sm:grid-cols-2 lg:grid-cols-6 items-stretch">
          {RQS.map((rq) => {
            const Icon = rq.icon;
            return (
              <a
                key={rq.id}
                href={rq.href}
                className="group relative overflow-hidden rounded-2xl border border-border bg-card shadow-sm hover:shadow-md transition"
              >
                {/* Background image */}
                <div
                  className="absolute inset-0 bg-cover bg-center"
                  style={{ backgroundImage: `url(${rq.bgImage})` }}
                  aria-hidden="true"
                />
                {/* Soft overlay to keep text readable */}
                <div className="absolute inset-0 bg-gradient-to-b from-background/40 via-background/70 to-background/95" />

                <div className="relative p-5 h-full flex flex-col">
                  {/*<div className="flex items-start gap-3">
                    <div className="w-10 h-10 rounded-xl bg-primary/10 flex items-center justify-center flex-shrink-0 border border-primary/10">
                      <Icon className="w-5 h-5 text-primary" />
                    </div>
                    <div>
                      <div className="text-sm font-semibold text-foreground">{rq.title}</div>
                      {/*<div className="text-xs text-muted-foreground mt-0.5">Click to navigate</div>
                    </div>
                  </div>*/}

                  <p className="mt-4 text-sm text-muted-foreground leading-relaxed">
                    {rq.question}
                  </p>

                  {/*<div className="mt-auto pt-4">
                    <span className="inline-flex items-center rounded-xl px-3 py-2 text-sm text-foreground bg-muted/40 border border-border/60 group-hover:bg-muted/60 transition">
                      Open section
                    </span>
                  </div>*/}
                </div>
              </a>
            );
          })}
        </div>

      </div>
    </section>
  );
};

export default RQSection;
