// The only place that touches BokehJS. Figures come from the API as json_item
// snapshots with named models; this module embeds them, applies data frames to
// the named models and serves the callbacks the figures make to `window.hdv`.
import { Frame } from "../api/frames";
import { sendJson } from "../api/client";

/* Minimal typing of the BokehJS objects we use (public API only, except the
   hover view redraw noted below). */
interface Model {
  [key: string]: any;
  setv(attrs: Record<string, unknown>): void;
}
interface Document {
  get_model_by_name(name: string): Model | null;
}
declare global {
  interface Window {
    Bokeh: any;
    hdv?: Bridge;
  }
}

type Listener = (payload: any) => void;

class Bridge {
  private listeners = new Map<string, Set<Listener>>();
  private names = new Map<string, Map<number, string>>(); // layer -> row -> barcode
  private rowVersions = new Map<string, number>();
  private pending = new Map<string, Set<number>>();
  private asked = new Set<string>();
  private scheduled = false;
  readonly docs = new Map<string, Document>();

  on(event: string, listener: Listener): () => void {
    if (!this.listeners.has(event)) this.listeners.set(event, new Set());
    this.listeners.get(event)!.add(listener);
    return () => this.listeners.get(event)?.delete(listener);
  }

  emit(event: string, payload: unknown = {}): void {
    if (this.quiet && event.endsWith("selection")) return;
    this.listeners.get(event)?.forEach((listener) => listener(payload));
  }

  // -- selections set from the server: not echoed back --------------------------
  private quiet = false;

  /** Set a source's selection without reporting it (it came from the server). */
  selectQuietly(source: Model | null, rows: number[]): void {
    if (!source) return;
    this.quiet = true;
    try {
      source.selected.indices = rows;
    } finally {
      this.quiet = false;
    }
  }

  /** Show `layer`'s selection (e.g. an ROI's members) on the map and tell the
   * other views (the embedding) that it changed. */
  showSelection(layer: string, rows: number[]): void {
    this.selectQuietly(this.model(layer, "source"), rows);
    this.emit("selection-changed", { layer });
  }

  // -- hover barcodes: fetched once per row position, batched per frame ------
  setRowVersion(layer: string, version: number): void {
    if (this.rowVersions.get(layer) !== version) {
      this.rowVersions.set(layer, version);
      this.names.set(layer, new Map());
      this.pending.delete(layer);
    }
  }

  rowVersion(layer: string): number | undefined {
    return this.rowVersions.get(layer);
  }

  hoverName(layer: string, row: number): string {
    const cached = this.names.get(layer)?.get(row);
    if (cached !== undefined) return cached;
    const key = `${layer}:${this.rowVersions.get(layer)}:${row}`;
    if (!this.asked.has(key)) {
      this.asked.add(key);
      if (!this.pending.has(layer)) this.pending.set(layer, new Set());
      this.pending.get(layer)!.add(row);
      if (!this.scheduled) {
        this.scheduled = true;
        requestAnimationFrame(() => this.flushNames());
      }
    }
    return "…";
  }

  private async flushNames(): Promise<void> {
    this.scheduled = false;
    const batches = [...this.pending.entries()];
    this.pending.clear();
    for (const [layer, rows] of batches) {
      const version = this.rowVersions.get(layer);
      if (version === undefined) continue;
      const reply = await sendJson<{ row_version: number; rows: [number, string][] }>(
        "POST", this.namesPath(layer), { row_version: version, rows: [...rows] });
      if (reply.row_version !== this.rowVersions.get(layer)) continue;
      const cache = this.names.get(layer)!;
      reply.rows.forEach(([row, name]) => cache.set(row, name));
      this.redrawTooltip(layer);
    }
  }

  /** Names come from the map layer, or from a layer's embedding ("embedding:<layer>"). */
  private namesPath(key: string): string {
    const [kind, layer] = key.startsWith("embedding:") ? ["embedding", key.slice(10)] : ["map", key];
    return `layers/${encodeURIComponent(layer)}/${kind === "embedding" ? "embedding/names" : "names"}`;
  }

  private redrawTooltip(layer: string): void {
    // Re-run the hover at the last mouse position so the name appears without
    // moving the mouse (BokehJS internals; harmless if they change).
    const hover = layer.startsWith("embedding:")
      ? this.docs.get("embedding")?.get_model_by_name("embedding:hover") ?? null
      : this.model(layer, "hover");
    try {
      const view = hover && window.Bokeh.index.find_one(hover);
      if (view?._current_sxy) view._inspect(...view._current_sxy);
    } catch {
      /* ignore */
    }
  }

  // -- figures ------------------------------------------------------------------
  private roots = new Map<string, { view: any; element: HTMLElement }>();

  async embed(item: any, element: HTMLElement, key: string): Promise<Document> {
    const previous = this.roots.get(key);
    try {
      if (previous?.element === element) previous.view?.remove();   // re-embedded in place
    } catch {
      /* already gone */
    }
    element.innerHTML = "";
    const views = await window.Bokeh.embed.embed_item(item, element);
    const view = Array.isArray(views) ? views[0] : views;
    this.roots.set(key, { view, element });
    const doc: Document = view?.model?.document ?? window.Bokeh.documents?.[window.Bokeh.documents.length - 1];
    this.docs.set(key, doc);
    return doc;
  }

  /** Opacity of the tissue image (tile or single-image renderer of the map). */
  setImageAlpha(alpha: number): void {
    const doc = this.docs.get("map");
    const tiles = doc?.get_model_by_name("image:tiles");
    if (tiles) tiles.alpha = alpha;
    const single = doc?.get_model_by_name("image:url");
    if (single) single.glyph.global_alpha = alpha;
  }

  model(layer: string, part: string): Model | null {
    return this.docs.get("map")?.get_model_by_name(`${layer}:${part}`) ?? null;
  }

  /** Apply an API frame to the layer's named models. */
  applyFrame(frame: Frame): void {
    const { header, ops, extra } = frame;
    const layer = header.layer;
    const source = this.model(layer, "source");
    const renderer = this.model(layer, "renderer");
    if (!source || !renderer) return;
    for (const op of ops) {
      if (op.op === "append") source.stream(op.columns);
      else if (op.op === "replace") source.data = { ...op.columns };
      else if (op.op === "set") source.data = { ...source.data, [op.name]: op.column };
    }
    const constants = header.constants ?? {};
    if (Object.keys(constants).length) {
      for (const glyph of [renderer.glyph, renderer.hover_glyph, renderer.nonselection_glyph,
                           renderer.selection_glyph, renderer.muted_glyph]) {
        if (glyph && typeof glyph === "object" && "setv" in glyph) {
          const attrs = Object.fromEntries(Object.entries(constants).filter(([name]) => name in glyph));
          if (Object.keys(attrs).length) glyph.setv(attrs);
        }
      }
    }
    const mapper = this.model(layer, "mapper");
    if (mapper && header.palette?.length) {
      mapper.setv({ palette: header.palette, low: -0.5, high: header.palette.length - 0.5 });
    }
    for (const [attr, palette] of Object.entries(header.palettes ?? {})) {
      const extra = this.model(layer, `mapper:${attr}`);
      if (extra && palette.length) extra.setv({ palette, low: -0.5, high: palette.length - 0.5 });
    }
    const filter = this.model(layer, "filter");
    if (filter) {
      const mask = extra.visible as ArrayLike<number> | undefined;
      filter.booleans = header.filtered && mask ? Array.from(mask, (value) => value !== 0) : null;
    }
    const hover = this.model(layer, "hover");
    if (hover && header.tooltips) hover.tooltips = header.tooltips;
    renderer.visible = header.visible;
    this.setRowVersion(layer, header.row_version);
  }
}

export const bridge = new Bridge();
window.hdv = bridge;
