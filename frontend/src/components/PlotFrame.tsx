type PlotSize = "sm" | "md" | "lg" | "xl";

const SIZE_CLASSES: Record<PlotSize, string> = {
  sm: "h-[360px]",
  md: "h-[520px]",
  lg: "h-[700px]",
};

interface PlotFrameProps {
  src?: string;
  title: string;
  size?: PlotSize;
  caption?: string;
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
      {src ? (
        <iframe
          src={src}
          title={title}
          loading="lazy"
          className={`w-full ${SIZE_CLASSES[size]} rounded-xl border border-border`}
        />
      ) : (
        <div
          className={`w-full ${SIZE_CLASSES[size]} flex items-center justify-center rounded-xl border border-dashed border-border bg-muted/30`}
        >
          <p className="text-sm text-muted-foreground text-center">
            {placeholder ?? "Plot forthcoming"}
          </p>
        </div>
      )}

      {caption && (
        <p className="mt-3 text-xs text-center text-muted-foreground">
          {caption}
        </p>
      )}
    </div>
  );
};

export default PlotFrame;
