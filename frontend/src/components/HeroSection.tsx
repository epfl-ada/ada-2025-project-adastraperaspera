import { Microscope, Database, BarChart3, Dna, HardDrive } from "lucide-react";

const HeroSection = () => {
  return (
    <section className="relative min-h-[85vh] flex items-center justify-center overflow-hidden">
      {/* Background gradient */}
      <div 
        className="absolute inset-0"
        style={{ background: 'var(--gradient-hero)' }}
      />
      
      {/* Grid pattern overlay */}
      <div className="absolute inset-0 opacity-10">
        <div 
          className="h-full w-full"
          style={{
            backgroundImage: `
              linear-gradient(hsl(0 0% 100% / 0.1) 1px, transparent 1px),
              linear-gradient(90deg, hsl(0 0% 100% / 0.1) 1px, transparent 1px)
            `,
            backgroundSize: '60px 60px'
          }}
        />
      </div>

      {/* Floating elements */}
      <div className="absolute top-20 left-[15%] w-32 h-32 bg-accent/20 rounded-full blur-3xl animate-pulse-glow" />
      <div className="absolute bottom-32 right-[20%] w-48 h-48 bg-primary-foreground/10 rounded-full blur-3xl" />

      <div className="relative z-10 container mx-auto px-6 text-center">
        <div className="animate-fade-up">
          <div className="inline-flex items-center gap-2 px-4 py-2 rounded-full bg-primary-foreground/10 border border-primary-foreground/20 text-primary-foreground/90 text-sm mb-8">
            <Microscope className="w-4 h-4" />
            <span>Spatial Transcriptomics Analysis</span>
          </div>
        </div>

        <h1 className="animate-fade-up delay-100 text-5xl md:text-7xl font-bold text-primary-foreground mb-6 leading-tight">
          Aβ Plaque Impact on
          <br />
          <span className="text-primary-foreground/80">Brain Microenvironment</span>
        </h1>

        <p className="animate-fade-up delay-200 text-xl md:text-2xl text-primary-foreground/70 max-w-3xl mx-auto mb-12 leading-relaxed">
          Investigating how amyloid-beta plaques influence the surrounding tissue 
          through quantitative spatial transcriptomics modeling
        </p>

        <div className="animate-fade-up delay-300 flex flex-wrap justify-center gap-4 mb-16">
          <div className="flex items-center gap-2 px-5 py-3 rounded-xl bg-primary-foreground/10 backdrop-blur-sm border border-primary-foreground/20">
            <Database className="w-5 h-5 text-accent" />
            <span className="text-primary-foreground font-medium">78.9M Transcripts</span>
          </div>
          <div className="flex items-center gap-2 px-5 py-3 rounded-xl bg-primary-foreground/10 backdrop-blur-sm border border-primary-foreground/20">
            <BarChart3 className="w-5 h-5 text-chart-3" />
            <span className="text-primary-foreground font-medium">347 Genes</span>
          </div>
          <div className="flex items-center gap-2 px-5 py-3 rounded-xl bg-primary-foreground/10 backdrop-blur-sm border border-primary-foreground/20">
            <Dna className="w-5 h-5 text-chart-5" />
            <span className="text-primary-foreground font-medium">351,714 Cells</span>
          </div>
          <div className="flex items-center gap-2 px-5 py-3 rounded-xl bg-primary-foreground/10 backdrop-blur-sm border border-primary-foreground/20">
            <HardDrive className="w-5 h-5 text-chart-2" />
            <span className="text-primary-foreground font-medium">34.7 GB</span>
          </div>
        </div>

        <div className="animate-fade-up delay-400">
          <a 
            href="#introduction" 
            className="inline-flex items-center gap-2 text-primary-foreground/70 hover:text-primary-foreground transition-colors"
          >
            <span>Explore Analysis</span>
            <svg className="w-5 h-5 animate-bounce" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 14l-7 7m0 0l-7-7m7 7V3" />
            </svg>
          </a>
        </div>
      </div>
    </section>
  );
};

export default HeroSection;
