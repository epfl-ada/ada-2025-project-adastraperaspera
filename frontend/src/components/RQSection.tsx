import React, { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";

type RQ = {
  id: string;    // "rq-1"
  title: string; // main label
};

const RQS: RQ[] = [
  { id: "rq-1", title: "Predicting Cell Composition" },
  { id: "rq-2", title: "Predicting Gene Expression" },
  { id: "rq-3", title: "Predicting Gene Expression" },
  { id: "rq-4", title: "Predicting Gene Expression" },
  { id: "rq-5", title: "Predicting Plaque Distance" },
  { id: "rq-6", title: "Predicting Plaque Distance" },
  { id: "rq-7", title: "Predicting Plaque Distance" },
];

function scrollToHash(hash: string) {
  const id = hash.replace("#", "");
  if (!id) return;

  const el = document.getElementById(id);
  if (!el) return;

  el.scrollIntoView({ behavior: "smooth", block: "start" });
}

export default function RQSection() {
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
    <section className="w-full">
      <div className="grid grid-cols-1 gap-4 lg:grid-cols-[340px_1fr]">
        {/* Big box */}
        <div className="rounded-2xl border border-white/10 bg-white/5 p-5 shadow-sm">
          <div className="text-xs uppercase tracking-wide text-white/70">Navigation</div>
          <h2 className="mt-2 text-2xl font-semibold leading-tight">Research Questions</h2>
          <p className="mt-2 text-sm text-white/70">
            Click a card to jump to the corresponding section below.
          </p>
        </div>

        {/* Horizontal cards */}
        <div className="rounded-2xl border border-white/10 bg-white/5 p-4 shadow-sm">
          <div className="flex items-center justify-between">
            <div className="text-sm font-medium text-white/80">Jump to</div>
            <div className="text-xs text-white/60">Scroll →</div>
          </div>

          <div className="mt-3 flex gap-3 overflow-x-auto pb-2 pr-1 [-ms-overflow-style:none] [scrollbar-width:none] [&::-webkit-scrollbar]:hidden">
            {RQS.map((rq, idx) => (
              <button
                key={rq.id}
                onClick={() => onSelect(rq.id)}
                className="group min-w-[260px] max-w-[340px] flex-shrink-0 rounded-2xl border border-white/10 bg-white/5 p-4 text-left transition hover:bg-white/10 hover:border-white/20 focus:outline-none focus:ring-2 focus:ring-white/25"
              >
                <div className="text-xs font-semibold text-white/60">RQ {idx + 1}</div>
                <div className="mt-2 text-base font-semibold leading-snug">{rq.title}</div>
                <div className="mt-3 text-xs text-white/60 group-hover:text-white/75">
                  Go to section →
                </div>
              </button>
            ))}
          </div>
        </div>
      </div>
    </section>
  );
}
