// ROI drawing on the map. One drawing at a time across layers (as the Bokeh
// UI's drawing controller): arming a tool cancels the previous owner's drawing.
// Shapes are in map coordinates; the server turns them into ROIs.
import { bridge } from "../bokeh/bridge";
import { mapPlot } from "./layers";

export type DrawMode = "box" | "lasso" | "radius";
export type Drawn =
  | { kind: "box"; box: [number, number, number, number] }
  | { kind: "lasso"; xs: number[]; ys: number[] }
  | { kind: "point"; centre: [number, number] };

// Polygons, each a list of rings (exterior first), as in Bokeh's multi_polygons
export interface Outline { xs: number[][][]; ys: number[][][] }
export interface SavedOutline extends Outline { colour: string }

// One owner (a layer's ROI panel) holds the drawing tool and the draft; a new
// owner takes them over and the previous one is told it lost them.
let owner: { id: string; lost(): void } | null = null;
let stopDrawing: (() => void) | null = null;

function mapModel(name: string): any {
  return bridge.docs.get("map")?.get_model_by_name(name) ?? null;
}

/** Become the owner of the drawing and draft (e.g. before a seed preview). */
export function claim(id: string, lost: () => void): void {
  if (owner && owner.id !== id) {
    const previous = owner;
    disown(previous.id);
    previous.lost();
  }
  owner = { id, lost };
}

/** Drop `id`'s drawing tool and draft, if it still owns them. */
export function disown(id: string): void {
  if (owner?.id !== id) return;
  stopDrawing?.();
  showDraft(id, null);
  owner = null;
}

/** Arm the map tool for `mode`; `done` gets the drawn shape. */
export function arm(id: string, mode: DrawMode, done: (shape: Drawn) => void, lost: () => void): boolean {
  claim(id, lost);
  stopDrawing?.();
  const tool = mapModel({ box: "roi_box", lasso: "roi_lasso", radius: "roi_tap" }[mode]);
  const plot = mapPlot();
  if (!tool || !plot) return false;
  const gesture = mode === "radius" ? "tap" : "pan";
  const previous = plot.toolbar.tools.filter((t: any) => t.active && [t.event_type].flat().includes(gesture));
  const finish = (shape: Drawn) => {
    stopDrawing?.();
    done(shape);
  };
  const offGeometry = bridge.on("geometry", ({ geometry, final }) => {
    if (!final || !tool.active || !geometry) return;
    if (mode === "box" && geometry.type === "rect") {
      finish({ kind: "box", box: [geometry.x0, geometry.y0, geometry.x1, geometry.y1] });
    } else if (mode === "lasso" && geometry.type === "poly") {
      finish({ kind: "lasso", xs: Array.from(geometry.x as ArrayLike<number>), ys: Array.from(geometry.y as ArrayLike<number>) });
    }
  });
  const offTap = bridge.on("tap", ({ x, y }) => {
    if (tool.active) finish({ kind: "point", centre: [x, y] });
  });
  // Give the map its previous tools back once the shape is drawn or dropped
  stopDrawing = () => {
    stopDrawing = null;
    offGeometry();
    offTap();
    tool.active = false;
    previous.forEach((t: any) => (t.active = true));
  };
  tool.active = true;
  return true;
}

/** The dashed outline of `id`'s ROI being drawn (null clears it). */
export function showDraft(id: string, outline: Outline | null): void {
  const draft = mapModel("draft");
  if (draft && owner?.id === id) draft.data = outline ? { xs: [outline.xs], ys: [outline.ys] } : { xs: [], ys: [] };
}

/** A layer's saved ROI outlines, in their colours. */
export function showOutlines(layer: string, outlines: SavedOutline[]): void {
  const source = bridge.model(layer, "rois");
  if (source) {
    source.data = { xs: outlines.map((o) => o.xs), ys: outlines.map((o) => o.ys), color: outlines.map((o) => o.colour) };
  }
}

export function focus(x0: number, x1: number, y0: number, y1: number): void {
  const plot = mapPlot();
  if (!plot) return;
  plot.x_range.setv({ start: x0, end: x1 });
  plot.y_range.setv({ start: y0, end: y1 });
}
