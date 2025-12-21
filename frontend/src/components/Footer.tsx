import { Github } from "lucide-react";

const PARTICIPANTS = [
  { name: "Alexander Sharipov", github: "https://github.com/Alex-T-Sharipov" },
  { name: "Sogand Salehi", github: "https://github.com/sogandstormesalehi" },
  { name: "Zayed Kriem", github: "https://github.com/ZayedK1" },
  { name: "Walid Sofiane", github: "https://github.com/Walsof-14" },
  { name: "Rosa Mayila", github: "https://github.com/rosbotmay" },
];

const Footer = () => {
  return (
    <footer className="py-12 bg-foreground text-background">
      <div className="container mx-auto px-6">
        {/* Top row */}
        <div className="flex flex-col md:flex-row items-center justify-between gap-8">
          <div>
            <h3 className="text-xl font-bold mb-2">Xenium Spatial Analysis</h3>
            <p className="text-background/60 text-sm">
              High-resolution spatial transcriptomics platform
            </p>
          </div>

          <div className="flex items-center gap-8">
            <a
              href="#introduction"
              className="text-background/70 hover:text-background transition-colors text-sm"
            >
              Overview
            </a>
            <a
              href="#rq-1"
              className="text-background/70 hover:text-background transition-colors text-sm"
            >
              Methods
            </a>
            <a
              href="https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=Xenium%20In%20Situ%20Analysis%20of%20Alzheimer%27s%20Disease%20Mouse%20Model%20Brain%20Coronal%20Sections%20from%20One%20Hemisphere%20Over%20a%20Time%20Course&lastTouchOfferType=Dataset&product=chromium&redirectUrl=%2Fdatasets%2Fxenium-in-situ-analysis-of-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard"
              target="_blank"
              rel="noopener noreferrer"
              className="text-background/70 hover:text-background transition-colors text-sm"
            >
              Data Access
            </a>
          </div>
        </div>

        {/* Participants */}
        <div className="mt-8 pt-6 border-t border-background/10">
          <div className="flex flex-wrap justify-center gap-6">
            {PARTICIPANTS.map((p) => (
              <a
                key={p.name}
                href={p.github}
                target="_blank"
                rel="noopener noreferrer"
                className="flex items-center gap-2 text-background/70 hover:text-background transition-colors text-sm"
              >
                <Github className="w-4 h-4" />
                <span>{p.name}</span>
              </a>
            ))}
          </div>
        </div>

        {/* Copyright */}
        <div className="mt-8 text-center">
          <p className="text-background/50 text-sm">
            (c) 2025 Spatial Genomics Research Team Adastraperaspera.
          </p>
        </div>
      </div>
    </footer>
  );
};

export default Footer;
