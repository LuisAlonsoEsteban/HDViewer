// Keeps each layer's BokehJS data source in step with the server: frames are
// applied in order; a frame built on another data version triggers a resync.
import { ApiError, frame, getJson, sendJson } from "../api/client";
import { Frame } from "../api/frames";
import { bridge } from "../bokeh/bridge";
import { useApp } from "../store/app";

const versions = new Map<string, number>();
const inFlight = new Set<string>();
const queued = new Set<string>();
let seq = 0;
let host: HTMLElement | null = null;

export function setMapHost(element: HTMLElement | null) {
  host = element;
}

export function mapPlot(): any {
  return bridge.docs.get("map")?.get_model_by_name("map") ?? null;
}

export async function applyLayerFrame(result: Frame): Promise<void> {
  const layer = result.header.layer;
  if (result.ops.length && result.header.base !== (versions.get(layer) ?? 0)) {
    const state = await frame("GET", `layers/${encodeURIComponent(layer)}/state`);
    bridge.applyFrame(state);
    versions.set(layer, state.header.version);
    return;
  }
  bridge.applyFrame(result);
  versions.set(layer, result.header.version);
}

/** Load what the current view needs for `layer` (one request in flight per layer). */
export async function requestViewport(layer: string): Promise<void> {
  const plot = mapPlot();
  if (!plot) return;
  if (inFlight.has(layer)) {
    queued.add(layer);
    return;
  }
  inFlight.add(layer);
  try {
    const result = await frame("POST", `layers/${encodeURIComponent(layer)}/viewport`, {
      seq: ++seq,
      x0: plot.x_range.start, x1: plot.x_range.end, y0: plot.y_range.start, y1: plot.y_range.end,
      width: host?.clientWidth || 800, height: host?.clientHeight || 800,
    });
    await applyLayerFrame(result);
  } catch (exc) {
    if (!(exc instanceof ApiError && exc.status === 409)) console.error(exc);
  } finally {
    inFlight.delete(layer);
    if (queued.delete(layer)) requestViewport(layer);
  }
}

export async function sendSelection(layer: string, rows: number[]): Promise<void> {
  const rowVersion = bridge.rowVersion(layer);
  if (rowVersion === undefined) return;
  try {
    await sendJson("POST", `layers/${encodeURIComponent(layer)}/selection`, { row_version: rowVersion, rows });
    bridge.emit("selection-changed", { layer });          // e.g. the embedding highlights it
  } catch (exc) {
    if (!(exc instanceof ApiError && exc.status === 409)) console.error(exc);
  }
}

/** Re-read a layer after an analysis changed it on the server, then reload the view. */
export async function refreshLayer(layer: string): Promise<void> {
  const state = await frame("GET", `layers/${encodeURIComponent(layer)}/state`);
  bridge.applyFrame(state);
  versions.set(layer, state.header.version);
  await requestViewport(layer);
}

export async function setStyle(layer: string, changes: Record<string, unknown>): Promise<void> {
  await applyLayerFrame(await frame("PATCH", `layers/${encodeURIComponent(layer)}/style`, changes));
}

/** PATCH a layer resource that answers with a frame (clusters, groups). */
export async function patchFrame(path: string, changes: Record<string, unknown>): Promise<void> {
  await applyLayerFrame(await frame("PATCH", path, changes));
}

/** Rebuild the map figure at the current view (its glyphs changed, e.g. from the
 * Parameters tab), then reload every layer's data into the new figure. */
export async function rebuildMap(): Promise<void> {
  const plot = mapPlot();
  if (!plot || !host) return;
  const view = new URLSearchParams({
    x0: String(plot.x_range.start), x1: String(plot.x_range.end),
    y0: String(plot.y_range.start), y1: String(plot.y_range.end),
  });
  const item = await getJson(`figures/map?${view}`);
  await bridge.embed(item, host, "map");
  bridge.setImageAlpha(imageAlpha);
  versions.clear();
  await Promise.all((useApp.getState().workspace?.layers ?? []).map(({ id: layer }) => refreshLayer(layer)));
}

let imageAlpha = 1;

export const imageOpacity = () => imageAlpha;

/** Opacity of the tissue image; kept across figure rebuilds. */
export function setImageAlpha(alpha: number): void {
  imageAlpha = alpha;
  bridge.setImageAlpha(alpha);
}
