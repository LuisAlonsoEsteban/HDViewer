import { useEffect, useRef, useState } from "react";
import { Autocomplete, Button, Group, SegmentedControl, Select, Stack, Text } from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { useApp } from "../store/app";
import { BokehFigure } from "./BokehFigure";

interface State { method: string | null; n: number; version: number; methods: Record<string, string>; default_colour_by: string }

const embeddingSource = () => bridge.docs.get("embedding")?.get_model_by_name("embedding:source") ?? null;

/** The Bokeh UI's Embeddings panel: PCA / UMAP / t-SNE of a layer, coloured by
 * any column; selections are shared with the map. */
export function EmbeddingPanel() {
  const layers = useApp((s) => s.workspace?.layers ?? []);
  const [layer, setLayer] = useState(layers[0]?.id ?? "");
  const job = useApp((s) => s.jobs[layer]?.embedding);
  const clustering = useApp((s) => s.jobs[layer]?.clustering);
  const version = useApp((s) => s.analysisVersion[layer] ?? 0);
  const [state, setState] = useState<State>();
  const [method, setMethod] = useState<string | null>("pca");
  const [colourBy, setColourBy] = useState("");
  const [query, setQuery] = useState("");
  const [options, setOptions] = useState<string[]>([]);
  const [item, setItem] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);
  const shownColour = useRef("");
  const seenClustering = useRef<string | undefined>(undefined);
  const base = `layers/${encodeURIComponent(layer)}/embedding`;
  const key = `embedding:${layer}`;

  // State; a new embedding (version) is drawn with the current colour
  useEffect(() => {
    getJson<State>(base).then((s) => {
      setState(s);
      // After clustering, show the new clusters (as the Bokeh panel)
      const done = clustering?.state === "done" ? clustering.message : undefined;
      if (done !== seenClustering.current) {
        seenClustering.current = done;
        if (done && s.default_colour_by) setColourBy(s.default_colour_by);
      }
      if (!colourBy && s.default_colour_by) setColourBy(s.default_colour_by);
    });
  }, [base, version]);

  useEffect(() => {
    setItem(null);
    if (!state?.method) return;
    getJson<any>(`${base}/figure?colour_by=${encodeURIComponent(colourBy)}`).then((figure) => {
      shownColour.current = colourBy;
      bridge.setRowVersion(key, figure.version);
      setItem(figure);
    }).catch((exc) => setError(exc.message));
  }, [state?.method, state?.n, job?.state === "done" ? job.message : "", base]);

  // Recolour in place (keeps the zoom)
  useEffect(() => {
    if (!item || colourBy === shownColour.current) return;
    getJson<{ version: number; palette: string[]; codes: number[] }>(`${base}/colours?colour_by=${encodeURIComponent(colourBy)}`)
      .then(({ version: v, palette, codes }) => {
        const source = embeddingSource();
        const mapper = bridge.docs.get("embedding")?.get_model_by_name("embedding:mapper");
        if (!source || !mapper) return;
        mapper.setv({ palette, low: -0.5, high: palette.length - 0.5 });
        source.data = { ...source.data, c: Uint16Array.from(codes) };
        bridge.setRowVersion(key, v);
        shownColour.current = colourBy;
      }).catch((exc) => setError(exc.message));
  }, [colourBy, item]);

  useEffect(() => {
    const timer = setTimeout(() => getJson<string[]>(`layers/${encodeURIComponent(layer)}/columns?q=${encodeURIComponent(query)}`)
      .then(setOptions).catch(() => setOptions([])), 150);
    return () => clearTimeout(timer);
  }, [query, layer]);

  // Selections: embedding -> map, and map (or an ROI) -> embedding
  useEffect(() => {
    const offEmbedding = bridge.on("embedding-selection", async ({ layer: from, indices }) => {
      if (from !== layer) return;
      const result = await sendJson<{ row_version: number; rows: number[] }>("POST", `${base}/selection`,
                                                                              { indices: Array.from(indices as number[]) });
      if (result.row_version === bridge.rowVersion(layer)) bridge.selectQuietly(bridge.model(layer, "source"), result.rows);
    });
    const offMap = bridge.on("selection-changed", async ({ layer: changed }) => {
      if (changed !== layer || !embeddingSource()) return;
      const { indices } = await getJson<{ indices: number[] }>(`${base}/selection`);
      bridge.selectQuietly(embeddingSource(), indices);
    });
    return () => { offEmbedding(); offMap(); };
  }, [layer, base]);

  const run = async () => {
    setError(null);
    try {
      await sendJson("POST", base, { method });
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  const running = job?.state === "running";
  const title = layers.find((l) => l.id === layer)?.title ?? layer;
  return (
    <Stack gap="sm">
      {layers.length > 1 && (
        <SegmentedControl size="xs" value={layer} onChange={(id) => { setLayer(id); setColourBy(""); }}
                          data={layers.map((l) => ({ value: l.id, label: l.title }))} />
      )}
      <Group grow align="end">
        <Select label="Method" value={method} onChange={setMethod} allowDeselect={false}
                data={Object.entries(state?.methods ?? { pca: "PCA", umap: "UMAP", tsne: "t-SNE" }).map(([value, label]) => ({ value, label }))} />
        <Button color="green" onClick={run} loading={running}>Run embedding</Button>
      </Group>
      <Autocomplete label="Color by" placeholder="Gene or metadata column" data={options} value={query || colourBy}
                    onChange={setQuery} onOptionSubmit={(name) => { setQuery(""); setColourBy(name); }} />
      <Text size="xs" c={job?.state === "failed" || error ? "red.4" : "dark.1"}>
        {error ?? (job && job.state !== "idle" ? job.message : state?.method ? "" : `Select a method and run an embedding for ${title}.`)}
      </Text>
      <BokehFigure item={item} docKey="embedding"
                   onReady={() => bridge.emit("selection-changed", { layer })} />
    </Stack>
  );
}
