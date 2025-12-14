import React, { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";
const base = import.meta.env.BASE_URL;

const researchQuestions = [
  {
    id: "rq-1",
    title: "Predicting Cell Composition",
    description: "How well can we infer cell-type composition from spatial transcriptomics features?",
    dotClass: "bg-chart-1",
  },
  {
    id: "rq-2",
    title: "Predicting Gene Expression",
    description: "Can plaque proximity explain changes in PIG gene expression patterns?",
    dotClass: "bg-chart-2",
  },
  {
    id: "rq-3",
    title: "Predicting Gene Expression",
    description: "Do regression models capture PIG expression as a function of plaque distance?",
    dotClass: "bg-chart-3",
  },
  {
    id: "rq-4",
    title: "Predicting Gene Expression",
    description: "How does plaque distance vary across annotated cell types?",
    dotClass: "bg-chart-4",
  },
  {
    id: "rq-5",
    title: "Predicting Plaque Distance",
    description: "Can we predict plaque distance using cell state / expression signatures?",
    dotClass: "bg-chart-5",
  },
  {
    id: "rq-6",
    title: "Predicting Plaque Distance",
    description: "Which models generalize best across samples and conditions?",
    dotClass: "bg-chart-6",
  },
  {
    id: "rq-7",
    title: "Predicting Plaque Distance",
    description: "Where do predictions fail, and what biological factors drive errors?",
    dotClass: "bg-chart-7",
  },
];

function scrollToHash(hash: string) {
  const id = hash.replace("#", "");
  if (!id) return;

  const el = document.getElementById(id);
  if (!el) return;

  el.scrollIntoView({ behavior: "smooth", block: "start" });
}

const RQSection = () => {
  const navigate = useNavigate();
  const location = useLocation();

  useEffect(() => {
    if (location.hash) {
      const t = window.setTimeout(() => scrollToHash(location.hash), 0);
      return () => window.clearTimeout(t);
    }
  }, [location.hash]);

  const onSelect = (id: string) => {
    navigate({ pathname: location.pathname, hash: `#${id}` }, { replace: false });
    scrollToHash(`#${id}`);
  };

  return (
    <section className="py-24 bg-muted/30">
      <div className="container mx-auto px-6">
        <div className="grid lg:grid-cols-2 gap-12 items-start">
          {/* Left: Big box */}
          <div>
            <h2 className="text-3xl md:text-4xl font-bold text-foreground mb-4">
              Research Questions
            </h2>
            <p className="text-lg text-muted-foreground mb-8">
              Seven questions guide the analysis. Click any item to jump to the corresponding
              section on this page.
            </p>

            <div className="mb-8">
                <img
                    src=`${base}rat-brain.svg`
                    className="w-full max-w-md rounded-2xl border border-border shadow-md"
                />
            </div>
          </div>

          {/* Right: RQ list styled like your “dot + text” rows */}
          <div className="space-y-4">
            {researchQuestions.map((rq, idx) => (
              <button
                key={rq.id}
                onClick={() => onSelect(rq.id)}
                className="w-full text-left flex items-center gap-3 p-4 rounded-xl bg-card border border-border hover:bg-muted/40 transition focus:outline-none focus:ring-2 focus:ring-ring"
              >
                {/* the dot like in your previous section (NOT custom colored by data) */}
                <div className={`w-3 h-3 rounded-full ${rq.dotClass}`} />

                <div className="min-w-0">
                  <div className="flex items-baseline gap-2">
                    <p className="font-medium text-foreground">
                      RQ {idx + 1} — {rq.title}
                    </p>
                  </div>
                  <p className="text-sm text-muted-foreground truncate">
                    {rq.description}
                  </p>
                </div>
              </button>
            ))}
          </div>
        </div>
      </div>
    </section>
  );
};

export default RQSection;
