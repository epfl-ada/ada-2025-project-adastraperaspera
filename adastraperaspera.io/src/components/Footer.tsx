const Footer = () => {
  return (
    <footer className="py-12 bg-foreground text-background">
      <div className="container mx-auto px-6">
        <div className="flex flex-col md:flex-row items-center justify-between gap-6">
          <div>
            <h3 className="text-xl font-bold mb-2">Xenium Spatial Analysis</h3>
            <p className="text-background/60 text-sm">
              High-resolution spatial transcriptomics platform
            </p>
          </div>
          
          <div className="flex items-center gap-8">
            <a href="#overview" className="text-background/70 hover:text-background transition-colors text-sm">
              Overview
            </a>
            <a href="#" className="text-background/70 hover:text-background transition-colors text-sm">
              Methods
            </a>
            <a href="#" className="text-background/70 hover:text-background transition-colors text-sm">
              Data Access
            </a>
            <a href="#" className="text-background/70 hover:text-background transition-colors text-sm">
              Contact
            </a>
          </div>
        </div>
        
        <div className="mt-8 pt-8 border-t border-background/10 text-center">
          <p className="text-background/50 text-sm">
            © 2024 Spatial Genomics Research. Analysis performed with Xenium In Situ platform.
          </p>
        </div>
      </div>
    </footer>
  );
};

export default Footer;
