import { useEffect, useRef, useState } from "react";
import { Autocomplete, Button, Checkbox, Group, NumberInput, Select, Stack, Text, TextInput } from "@mantine/core";
import { ApiError, getJson, sendJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { refreshLayer } from "../map/layers";
import { arm, claim, disown, DrawMode, Drawn, focus, Outline, showDraft, showOutlines } from "../map/roi";
import { Layer, useApp } from "../store/app";

interface ROIInfo extends Outline {
  id: string;
  name: string;
  kind: string;
  colour: string;
  count: number | null;
  visible: boolean;
  error: string | null;
  group: string | null;
  seed_group: string | null;
  seed_count: number | null;
  seed_barcode: string | null;
  group_excluded: boolean;
}
interface ROIList { rois: ROIInfo[]; unit: string; seed_groups: string[] }
interface Preview extends Outline { bounds: [number, number, number, number]; seeds: number | null }
type Shape = Drawn
  | { kind: "barcode"; barcode: string }
  | { kind: "group"; group: string; exclude_group: boolean };
type Message = { text: string; error: boolean } | null;

const KIND_LABELS: Record<string, string> = {
  box: "bounding box", lasso: "lasso", point: "seed + radius", barcode: "barcode + radius", group: "group + radius",
};

/** ROIs tool: draw regions, then group, select or focus them. */
export function ROIPanel({ layer }: { layer: Layer }) {
  const config = useApp((s) => s.config);
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const analysesChanged = useApp((s) => s.analysesChanged);
  const [list, setList] = useState<ROIList>();
  const [name, setName] = useState("");
  const [method, setMethod] = useState<DrawMode>("box");
  const [radius, setRadius] = useState<number | string>(config?.roi.default_radius ?? 50);
  const [seedSource, setSeedSource] = useState("map");
  const [barcode, setBarcode] = useState("");
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [seedGroup, setSeedGroup] = useState<string | null>(null);
  const [region, setRegion] = useState("full");
  const [draft, setDraft] = useState<{ shape: Shape; seeds: number | null } | null>(null);
  const [drawing, setDrawing] = useState(false);
  const [previewing, setPreviewing] = useState(false);
  const [saving, setSaving] = useState(false);
  const [message, setMessage] = useState<Message>({ text: "Choose a method and draw an ROI.", error: false });
  const preview = useRef<AbortController | null>(null);
  const base = `layers/${encodeURIComponent(layer.id)}`;
  const seeded = layer.aggregation_factors.length > 0;          // barcode and group seeds: bins layers
  const mode = method !== "radius" || !seeded ? "draw" : seedSource;   // draw | map | barcode | group

  useEffect(() => { getJson<ROIList>(`${base}/rois`).then(setList); }, [version, base]);
  // Saved outlines on the map (hidden with the layer)
  useEffect(() => {
    if (list) showOutlines(layer.id, layer.visible ? list.rois.filter((r) => r.visible) : []);
  }, [list, layer.visible, layer.id]);
  // A draft belongs to this panel: leaving the tab or hiding the layer drops it
  useEffect(() => () => cancel(false), []);
  useEffect(() => { if (!layer.visible) cancel(); }, [layer.visible]);

  // Local state only: the map's tool and draft are dropped by disown(), or
  // were taken over by another layer's panel (claim)
  const reset = (text = "Choose a method and draw an ROI.") => {
    preview.current?.abort();
    setDraft(null);
    setDrawing(false);
    setPreviewing(false);
    if (text) setMessage({ text, error: false });
  };
  const cancel = (say = true) => {
    disown(layer.id);
    reset(say ? undefined : "");
  };
  const lost = () => reset("Drawing moved to another layer.");
  const radiusValue = () => (typeof radius === "number" && radius > 0 ? radius : null);

  const showPreview = async (shape: Shape, zoom = false) => {
    preview.current?.abort();
    const controller = (preview.current = new AbortController());
    const r = radiusValue();
    if (shape.kind !== "box" && shape.kind !== "lasso" && r === null) {
      setMessage({ text: "Enter a positive radius.", error: true });
      return;
    }
    claim(layer.id, lost);
    setPreviewing(true);
    setMessage({ text: shape.kind === "group" ? "Combining the group's seed neighbourhoods…" : "Previewing…", error: false });
    try {
      const result = await sendJson<Preview>("POST", `${base}/rois/preview`, { ...shape, radius: r }, controller.signal);
      if (controller.signal.aborted) return;
      showDraft(layer.id, result);
      setDraft({ shape, seeds: result.seeds });
      setMessage({ text: "Draft ready. Save ROI or cancel.", error: false });
      if (zoom) {
        const [x0, y0, x1, y1] = result.bounds;
        const px = 0.05 * (x1 - x0), py = 0.05 * (y1 - y0);
        focus(x0 - px, x1 + px, y0 - py, y1 + py);
      }
    } catch (exc: any) {
      if (controller.signal.aborted || (exc instanceof ApiError && exc.status === 409)) return;   // superseded
      showDraft(layer.id, null);
      setDraft(null);
      setMessage({ text: exc.message, error: true });
    } finally {
      if (preview.current === controller) setPreviewing(false);
    }
  };

  const begin = () => {
    cancel(false);
    if (mode === "barcode") {
      if (!barcode.trim()) return setMessage({ text: "Type a native barcode or bin index and choose a suggestion.", error: true });
      return void showPreview({ kind: "barcode", barcode: barcode.trim() }, true);
    }
    if (mode === "group") {
      if (!seedGroup) return setMessage({ text: "Choose a seed group to preview one combined ROI.", error: true });
      return void showPreview({ kind: "group", group: seedGroup, exclude_group: region === "minus_group" });
    }
    if (method === "radius" && radiusValue() === null) return setMessage({ text: "Enter a positive radius.", error: true });
    const armed = arm(layer.id, method, (shape) => { setDrawing(false); showPreview(shape); }, lost);
    if (!armed) return setMessage({ text: "The map is not ready yet.", error: true });
    setDrawing(true);
    setMessage({
      text: `Drawing for ${layer.title}: ${method === "radius" ? "click the seed on the map." : "drag the region on the map."}`,
      error: false,
    });
  };

  // A new radius, group or region redraws the seed-based draft (debounced)
  useEffect(() => {
    if (!draft || draft.shape.kind === "box" || draft.shape.kind === "lasso") return;
    const shape = draft.shape.kind === "group"
      ? { ...draft.shape, group: seedGroup ?? draft.shape.group, exclude_group: region === "minus_group" }
      : draft.shape;
    const timer = setTimeout(() => showPreview(shape), config?.roi.preview_debounce_ms ?? 200);
    return () => clearTimeout(timer);
  }, [radius, seedGroup, region]);

  useEffect(() => {
    if (mode !== "barcode" || !barcode.trim()) return setSuggestions([]);
    const timer = setTimeout(() => getJson<string[]>(`${base}/barcodes?q=${encodeURIComponent(barcode.trim())}`)
      .then(setSuggestions).catch(() => setSuggestions([])), 150);
    return () => clearTimeout(timer);
  }, [barcode, mode, base]);

  const save = async () => {
    if (!draft) return;
    setSaving(true);
    setMessage({ text: "Computing ROI membership…", error: false });
    try {
      const result = await sendJson<ROIList>("POST", `${base}/rois`, { ...draft.shape, name, radius: radiusValue() });
      const roi = result.rois[result.rois.length - 1];
      cancel(false);
      setList(result);
      setName("");
      setMessage({ text: `${roi.name}: ${(roi.count ?? 0).toLocaleString()} observations`, error: false });
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    } finally {
      setSaving(false);
    }
  };

  const act = async (roi: ROIInfo, action: string) => {
    try {
      if (action === "focus") {
        const box = await sendJson<{ x0: number; x1: number; y0: number; y1: number }>("POST", `${base}/rois/${roi.id}/focus`);
        focus(box.x0, box.x1, box.y0, box.y1);
      } else if (action === "select" || action === "unselect") {
        const result = await sendJson<{ row_version: number; rows: number[]; selected: number }>(
          "POST", `${base}/rois/${roi.id}/${action}`);
        if (result.row_version === bridge.rowVersion(layer.id)) bridge.showSelection(layer.id, result.rows);
        setMessage({ text: `${result.selected.toLocaleString()} observations selected.`, error: false });
      } else if (action === "group") {
        const result = await sendJson<ROIList & { group: string }>("POST", `${base}/rois/${roi.id}/group`);
        setList(result);
        setMessage({ text: `Created group: ${result.group}`, error: false });
        await refreshLayer(layer.id);
        analysesChanged(layer.id);
      } else if (action === "delete") {
        const result = await sendJson<ROIList & { recoloured: boolean }>("DELETE", `${base}/rois/${roi.id}`);
        setList(result);
        setMessage({ text: roi.group ? `ROI and group '${roi.group}' deleted.` : "ROI deleted.", error: false });
        if (result.recoloured) await refreshLayer(layer.id);
        analysesChanged(layer.id);
      }
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    }
  };

  const outline = async (roi: ROIInfo, visible: boolean) => {
    setList((l) => l && { ...l, rois: l.rois.map((r) => (r.id === roi.id ? { ...r, visible } : r)) });
    setList(await sendJson<ROIList>("PATCH", `${base}/rois/${roi.id}`, { visible }));
  };

  const unit = list?.unit ?? "µm";
  const groups = list?.seed_groups ?? [];
  const busy = previewing || saving;
  const changeMethod = (update: () => void) => { cancel(); update(); };

  return (
    <Stack gap="sm">
      <TextInput label="ROI name (optional)" placeholder="ROI 1" value={name} onChange={(e) => setName(e.currentTarget.value)} />
      <Select label="Geometry" value={method} allowDeselect={false}
              data={[{ value: "box", label: "Bounding box" }, { value: "lasso", label: "Lasso" }, { value: "radius", label: "Seed + radius" }]}
              onChange={(value) => changeMethod(() => setMethod(value as DrawMode))} />
      {method === "radius" && (
        <NumberInput label={`Radius (${unit})`} min={0} decimalScale={3} value={radius} onChange={setRadius} />
      )}
      {method === "radius" && seeded && (
        <Select label="Seed source" value={seedSource} allowDeselect={false}
                data={[{ value: "map", label: "Click on map" }, { value: "barcode", label: "Find barcode / bin index" },
                       { value: "group", label: "All bins in a group" }]}
                onChange={(value) => changeMethod(() => setSeedSource(value ?? "map"))} />
      )}
      {mode === "barcode" && (
        <Autocomplete label="Barcode / bin index" placeholder="Type a native barcode…" data={suggestions} value={barcode}
                      onChange={(value) => { if (draft) cancel(); setBarcode(value); }}
                      onOptionSubmit={(value) => { setBarcode(value); showPreview({ kind: "barcode", barcode: value }, true); }} />
      )}
      {mode === "group" && (
        <>
          <Select label="Seed group" placeholder={groups.length ? "Choose a group" : "Create a group first"} data={groups}
                  value={seedGroup && groups.includes(seedGroup) ? seedGroup : null}
                  onChange={(value) => { setSeedGroup(value); if (value && !draft) showPreview({ kind: "group", group: value, exclude_group: region === "minus_group" }); }} />
          <Select label="ROI region" value={region} allowDeselect={false}
                  data={[{ value: "full", label: "Full ROI" }, { value: "minus_group", label: "ROI minus group" }]}
                  onChange={(value) => setRegion(value ?? "full")} />
        </>
      )}
      {draft?.shape.kind === "barcode" && <Text size="xs" c="dark.1">Seed: {draft.shape.barcode} · native bin centre.</Text>}
      {draft?.shape.kind === "group" && draft.seeds !== null && (
        <Text size="xs" c="dark.1">
          All {draft.seeds.toLocaleString()} native bins are seeds; their neighbourhoods form one ROI.
          {draft.shape.exclude_group && " Source-group bin shapes and members are excluded."}
        </Text>
      )}
      <Group gap="xs">
        <Button onClick={begin} disabled={!layer.visible || saving} loading={previewing || drawing}>
          {mode === "barcode" ? "Preview barcode" : mode === "group" ? "Preview group" : "Draw"}
        </Button>
        <Button color="green" onClick={save} disabled={!draft || previewing} loading={saving}>Save ROI</Button>
        <Button variant="default" onClick={() => cancel()} disabled={saving}>Cancel</Button>
      </Group>
      {message && <Text size="xs" c={message.error ? "red.4" : "dark.1"}>{message.text}</Text>}

      <Text fw={600} size="sm">Saved ROIs</Text>
      {list && list.rois.length === 0 && <Text size="xs" c="dark.2">No ROIs saved in this layer.</Text>}
      {list?.rois.map((roi) => (
        <Stack key={roi.id} gap={4} data-testid="roi">
          <div>
            <Text size="sm"><Text span fw={700} c={roi.colour}>{roi.name}</Text> · {KIND_LABELS[roi.kind] ?? roi.kind}</Text>
            <Text size="xs" c={roi.error ? "red.4" : "dark.1"}>
              {roi.error ? `Error: ${roi.error}` : roi.count === null ? "Computing membership…" : `${roi.count.toLocaleString()} observations`}
            </Text>
            {roi.seed_group !== null && (
              <Text size="xs" c="dark.2">
                {(roi.seed_count ?? 0).toLocaleString()} seeds from {roi.seed_group}{roi.group_excluded && " · source group excluded"}
              </Text>
            )}
            {roi.seed_barcode !== null && <Text size="xs" c="dark.2">Seed: {roi.seed_barcode}</Text>}
          </div>
          <Group gap={4}>
            <Button size="compact-xs" variant="default" onClick={() => act(roi, "focus")}>Focus</Button>
            <Button size="compact-xs" variant="default" disabled={!roi.count} onClick={() => act(roi, "select")}>Select</Button>
            <Button size="compact-xs" variant="default" disabled={!roi.count} onClick={() => act(roi, "unselect")}>Unselect</Button>
            <Button size="compact-xs" variant="default" disabled={!roi.count || roi.group !== null} onClick={() => act(roi, "group")}
                    title={roi.group ? `Group '${roi.group}'` : undefined}>
              {roi.group ? "Group created" : "Create group"}
            </Button>
            <Button size="compact-xs" color="red" variant="light" onClick={() => act(roi, "delete")}
                    aria-label={`Delete ${roi.name}`}>Delete</Button>
          </Group>
          <Checkbox size="xs" label="Outline" checked={roi.visible} onChange={(e) => outline(roi, e.currentTarget.checked)} />
        </Stack>
      ))}
    </Stack>
  );
}
