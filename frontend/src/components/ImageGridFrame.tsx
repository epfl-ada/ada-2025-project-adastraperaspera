import React, { useEffect, useMemo, useRef } from "react";
import Panzoom from "@panzoom/panzoom";
const base = import.meta.env.BASE_URL;

type CellKey = string;

type Props = {
  /** map key -> image URL (from /public). Example: "/images/full/wt-2.webp" */
  files: Record<CellKey, string>;
  /** map key -> [row, col] */
  keyToPos: Record<CellKey, [number, number]>;
  colTicks: string[];
  rowTicks: string[];
  caption?: string;

  /** Optional: flip horizontally for some keys */
  flipHKeys?: CellKey[];

  /** Visual tuning */
  cellSize?: number; // px
  gap?: number; // px
};

export default function ImageGridPlot({
  files,
  keyToPos,
  colTicks,
  rowTicks,
  caption,
  flipHKeys = [],
  cellSize = 220,
  gap = 10,
}: Props) {
  const wrapRef = useRef<HTMLDivElement | null>(null);
  const stageRef = useRef<HTMLDivElement | null>(null);

  const flipSet = useMemo(() => new Set(flipHKeys), [flipHKeys]);

  const { nRows, nCols } = useMemo(() => {
    const positions = Object.values(keyToPos);
    const nRows = Math.max(...positions.map(([r]) => r)) + 1;
    const nCols = Math.max(...positions.map(([, c]) => c)) + 1;
    return { nRows, nCols };
  }, [keyToPos]);

  // Setup pan/zoom
  useEffect(() => {
    if (!stageRef.current || !wrapRef.current) return;

    const panzoom = Panzoom(stageRef.current, {
      maxScale: 8,
      minScale: 1,
      contain: "outside",
      cursor: "grab",
    });

    const onWheel = (e: WheelEvent) => {
      // allow page scroll when user is not hovering the plot
      e.preventDefault();
      panzoom.zoomWithWheel(e);
    };

    wrapRef.current.addEventListener("wheel", onWheel, { passive: false });

    return () => {
      wrapRef.current?.removeEventListener("wheel", onWheel as any);
      panzoom.destroy();
    };
  }, []);

  // Layout sizes
  const plotW = nCols * cellSize + (nCols - 1) * gap;
  const plotH = nRows * cellSize + (nRows - 1) * gap;

  return (
    <figure className="w-full">
      <div
        ref={wrapRef}
        className="relative w-full overflow-hidden rounded-2xl border border-border bg-transparent"
        style={{ height: Math.min(640, plotH + 90) }} // viewport height (zoom/pan inside)
      >
        {/* Stage = the thing we pan/zoom */}
        <div
          ref={stageRef}
          className="absolute left-0 top-0"
          style={{
            width: plotW + 120, // room for y labels
            height: plotH + 80, // room for x labels
            transformOrigin: "0 0",
          }}
        >
          {/* Y tick labels */}
          <div
            className="absolute left-0 top-0"
            style={{ width: 110, height: plotH }}
          >
            {rowTicks.map((lab, i) => {
              // rowTicks top->bottom; map to y center
              const y = i * (cellSize + gap) + cellSize / 2;
              return (
                <div
                  key={lab}
                  className="absolute right-3 -translate-y-1/2 text-xs text-muted-foreground whitespace-nowrap"
                  style={{ top: y }}
                >
                  {lab}
                </div>
              );
            })}
            <div className="absolute left-2 top-1/2 -translate-y-1/2 -rotate-90 text-xs text-muted-foreground">
              Type
            </div>
          </div>

          {/* Grid area */}
          <div
            className="absolute"
            style={{ left: 110, top: 0, width: plotW, height: plotH }}
          >
            {Object.entries(keyToPos).map(([key, [r, c]]) => {
              const src = files[key];
              if (!src) return null;

              const x = c * (cellSize + gap);
              const y = r * (cellSize + gap);

              return (
                <div
                  key={key}
                  className="absolute overflow-hidden rounded-xl bg-transparent"
                  style={{
                    left: x,
                    top: y,
                    width: cellSize,
                    height: cellSize,
                  }}
                >
                  <img
                    src={src}
                    alt={key}
                    className="h-full w-full object-cover"
                    style={{
                      transform: flipSet.has(key) ? "scaleX(-1)" : undefined,
                    }}
                    draggable={false}
                    loading="lazy"
                  />
                </div>
              );
            })}
          </div>

          {/* X tick labels */}
          <div
            className="absolute"
            style={{ left: 110, top: plotH + 8, width: plotW, height: 60 }}
          >
            {colTicks.map((lab, i) => {
              const x = i * (cellSize + gap) + cellSize / 2;
              return (
                <div
                  key={lab}
                  className="absolute -translate-x-1/2 text-xs text-muted-foreground whitespace-nowrap"
                  style={{ left: x, top: 0 }}
                >
                  {lab}
                </div>
              );
            })}
            <div className="absolute left-1/2 -translate-x-1/2 top-6 text-xs text-muted-foreground">
              Age (months)
            </div>
          </div>
        </div>

        {/* Small hint overlay */}
        <div className="pointer-events-none absolute right-3 top-3 rounded-xl bg-background/70 px-3 py-2 text-xs text-muted-foreground backdrop-blur">
          Scroll to zoom · Drag to pan
        </div>
      </div>

      {caption ? (
        <figcaption className="mt-3 text-sm text-muted-foreground">
          {caption}
        </figcaption>
      ) : null}
    </figure>
  );
}
