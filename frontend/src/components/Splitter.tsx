import { useCallback, useState } from "react";
import type { CSSProperties, KeyboardEvent, PointerEvent } from "react";

/** A panel size (px) the user can drag, remembered for the browser session. */
export function useDragSize(key: string, initial: number, min: number, max: number) {
  const [size, setSize] = useState(() => {
    try {
      const saved = Number(sessionStorage.getItem(`hdv:size:${key}`));
      return saved >= min && saved <= max ? saved : initial;
    } catch {
      return initial;
    }
  });
  const update = useCallback((value: number) => {
    const clamped = Math.round(Math.min(max, Math.max(min, value)));
    setSize(clamped);
    try {
      sessionStorage.setItem(`hdv:size:${key}`, String(clamped));
    } catch {
      /* storage unavailable: the size just isn't remembered */
    }
  }, [key, min, max]);
  return [size, update] as const;
}

/** A drag handle (the old interface's split.js gutters). `axis` "x" resizes a
 * width, "y" a height; `sign` -1 when dragging towards the panel grows it
 * (e.g. the right panel, whose handle is on its left edge). Arrow keys also work. */
export function Splitter({ axis, size, onResize, sign = 1, label, style }: {
  axis: "x" | "y";
  size: number;
  onResize(size: number): void;
  sign?: 1 | -1;
  label: string;
  style?: CSSProperties;
}) {
  const start = (event: PointerEvent<HTMLDivElement>) => {
    event.preventDefault();
    const origin = axis === "x" ? event.clientX : event.clientY;
    const initial = size;
    const target = event.currentTarget;
    target.setPointerCapture(event.pointerId);
    const move = (e: globalThis.PointerEvent) =>
      onResize(initial + sign * ((axis === "x" ? e.clientX : e.clientY) - origin));
    const stop = () => {
      target.removeEventListener("pointermove", move);
      target.removeEventListener("pointerup", stop);
      target.removeEventListener("pointercancel", stop);
    };
    target.addEventListener("pointermove", move);
    target.addEventListener("pointerup", stop);
    target.addEventListener("pointercancel", stop);
  };
  const keys = (event: KeyboardEvent<HTMLDivElement>) => {
    const step = event.shiftKey ? 50 : 10;
    const grow = axis === "x" ? { ArrowRight: 1, ArrowLeft: -1 } : { ArrowDown: 1, ArrowUp: -1 };
    const direction = grow[event.key as keyof typeof grow];
    if (direction) {
      event.preventDefault();
      onResize(size + sign * direction * step);
    }
  };
  return (
    <div role="separator" aria-orientation={axis === "x" ? "vertical" : "horizontal"} aria-label={label}
         aria-valuenow={size} tabIndex={0} onPointerDown={start} onKeyDown={keys}
         className={`splitter splitter-${axis}`} style={style} />
  );
}
