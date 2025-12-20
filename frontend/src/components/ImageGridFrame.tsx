import React, { useEffect, useMemo, useRef } from "react";

type CellKey = string;

type Props = {
  files: Record<CellKey, string>;
  keyToPos: Record<CellKey, [number, number]>;
  colTicks: string[];
  rowTicks: string[];
  caption?: string;
  flipHKeys?: CellKey[];
  cellSize?: number;
  gap?: number;
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

  const plotW = nCols * cellSize + (nCols - 1) * gap;
  const plotH = nRows * cellSize + (nRows - 1) * gap;

  useEffect(() => {
    const wrap = wrapRef.current;
    const stage = stageRef.current;
    if (!wrap || !stage) return;

    let scale = 1;
    let x = 0;
    let y = 0;

    let isDragging = false;
    let lastX = 0;
    let lastY = 0;

    const clamp = (v: number, min: number, max: number) => Math.min(max, Math.max(min, v));

    const apply = () => {
      stage.style.transform = `translate(${x}px, ${y}px) scale(${scale})`;
    };

    apply();

    const onWheel = (e: WheelEvent) => {
      // Zoom on wheel, centered around cursor
      e.preventDefault();

      const rect = wrap.getBoundingClientRect();
      const cx = e.clientX - rect.left;
      const cy = e.clientY - rect.top;

      const prevScale = scale;
      const delta = -e.deltaY;

      // smooth zoom factor
      const factor = delta > 0 ? 1.08 : 1 / 1.08;
      scale = clamp(scale * factor, 1, 8);

      // adjust translation so zoom centers on cursor
      const k = scale / prevScale;
      x = cx - (cx - x) * k;
      y = cy - (cy - y) * k;

      apply();
    };

    const onPointerDown = (e: PointerEvent) => {
      isDragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      wrap.setPointerCapture(e.pointerId);
      wrap.style.cursor = "grabbing";
    };

    const onPointerMove = (e: PointerEvent) => {
      if (!isDragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;
      x += dx;
      y += dy;
      apply();
    };

    const onPointerUp = (e: PointerEvent) => {
      isDragging = false;
      wrap.releasePointerCapture(e.pointerId);
      wrap.style.cursor = "grab";
    };

    wrap.addEventListener("wheel", onWheel, { passive: false });
    wrap.addEventListener("pointerdown", onPointerDown);
    wrap.addEventListener("pointermove", onPointerMove);
    wrap.addEventListener("pointerup", onPointerUp);
    wrap.addEventListener("pointercancel", onPointerUp);

    wrap.style.cursor = "grab";

    return () => {
      wrap.removeEventListener("wheel", onWheel as any);
      wrap.removeEventListener("pointerdown", onPointerDown as any);
      wrap.removeEventListener("pointermove", onPointerMove as any);
      wrap.removeEventListener("pointerup", onPointerUp as any);
      wrap.removeEventListener("pointercancel", onPointerUp as any);
    };
  }, []);

  return (
    <figure className="w-full">
      <div
        ref={wrapRef}
        className="relative w-full overflow-hidden rounded-2xl border border-border bg-transparent"
        style={{ height: Math.min(640, plotH + 90) }}
      >
        <div
          ref={stageRef}
          className="absolute left-0 top-0"
          style={{
            width: plotW + 120,
            height: plotH + 80,
            transformOrigin: "0 0",
            willChange: "transform",
          }}
        >
          {/* Y labels */}
          <div className="absolute left-0 top-0" style={{ width: 110, height: plotH }}>
            {rowTicks.map((lab, i) => {
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

          {/* Grid */}
          <div className="absolute" style={{ left: 110, top: 0, width: plotW, height: plotH }}>
            {Object.entries(keyToPos).map(([key, [r, c]]) => {
              const src = files[key];
              if (!src) return null;

              const left = c * (cellSize + gap);
              const top = r * (cellSize + gap);

              return (
                <div
                  key={key}
                  className="absolute overflow-hidden rounded-xl bg-transparent"
                  style={{ left, top, width: cellSize, height: cellSize }}
                >
                  <img
                    src={src}
                    alt={key}
                    className="h-full w-full object-cover"
                    style={{ transform: flipSet.has(key) ? "scaleX(-1)" : undefined }}
                    draggable={false}
                    loading="lazy"
                  />
                </div>
              );
            })}
          </div>

          {/* X labels */}
          <div className="absolute" style={{ left: 110, top: plotH + 8, width: plotW, height: 60 }}>
            {colTicks.map((lab, i) => {
              const left = i * (cellSize + gap) + cellSize / 2;
              return (
                <div
                  key={lab}
                  className="absolute -translate-x-1/2 text-xs text-muted-foreground whitespace-nowrap"
                  style={{ left, top: 0 }}
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

        <div className="pointer-events-none absolute right-3 top-3 rounded-xl bg-background/70 px-3 py-2 text-xs text-muted-foreground backdrop-blur">
          Scroll to zoom · Drag to pan
        </div>
      </div>

      {caption ? (
        <figcaption className="mt-3 text-sm text-muted-foreground">{caption}</figcaption>
      ) : null}
    </figure>
  );
}
