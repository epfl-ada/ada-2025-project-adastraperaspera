const IntroductionSection2 = () => {
  return (
    <section id="introduction" className="py-24 bg-background">
      <div className="container mx-auto px-6 max-w-3xl">

        {/* Title */}
        <div className="mb-10">
          <h1 className="text-4xl md:text-5xl font-bold text-foreground mb-4">
            The Story of a Sick Mouse
          </h1>
          <p className="text-lg md:text-xl text-muted-foreground italic">
            Prologue: A brain that looks normal until you know where to look
          </p>
        </div>

        {/* Narrative text */}
        <div className="space-y-6 text-lg leading-relaxed text-foreground">
          <p>
            From far away, the slice looks like a brain slice is supposed to look:
            layers, curves, familiar anatomy.
          </p>

          <p>
            But in the transgenic mouse, the one engineered to overproduce amyloid
            precursor protein, something has been accumulating for months.
          </p>

          <p className="text-muted-foreground">
            Tiny lesions that don’t announce themselves in the morphology image
            alone.
          </p>
        </div>

      </div>
    </section>
  );
};

export default IntroductionSection2;
