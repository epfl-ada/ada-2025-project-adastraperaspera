type PlotSize = "sm" | "md" | "lg";

const SIZE_CLASSES: Record<PlotSize, string> = {
  sm: "h-[320px] md:h-[380px]",
  md: "h-[420px] md:h-[520px]",
  lg: "h-[520px] md:h-[680px]",
};

interface PlotFrameProps {
  src?: string;
  title: string;
  size?: PlotSize;
  caption?: React.ReactNode;
  placeholder?: string;
}

const PlotFrame = ({
  src,
  title,
  size = "md",
  caption,
  placeholder,
}: PlotFrameProps) => {
  return (
    <div className="w-full">
      {/* Title */}
      <h4 className="mb-2 text-sm font-medium text-foreground">
        {title}
      </h4>

      {/* Frame */}
      <div
        className={`relative w-full ${SIZE_CLASSES[size]} overflow-hidden rounded-xl border border-border`}
      >
        {src ? (
          <iframe
            src={src}
            title={title}
            loading="lazy"
            className="absolute inset-0 h-full w-full border-0"
            scrolling="no"
          />
        ) : (
          <div className="flex h-full items-center justify-center text-sm text-muted-foreground bg-muted/30">
            {placeholder ?? "Plot forthcoming"}
          </div>
        )}
      </div>

      {/* Caption */}
      {caption && (
        <p className="mt-2 text-xs text-muted-foreground leading-snug">
          {caption}
        </p>
      )}
    </div>
  );
};

export default PlotFrame;
