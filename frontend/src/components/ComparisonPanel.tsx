import { useEffect, useMemo, useState } from "react";
import {
  Button, Group, ScrollArea, Select, Slider, Stack, Table, Tabs, Text, TextInput, UnstyledButton,
} from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { requestViewport, setStyle } from "../map/layers";
import { JobProgress } from "./JobProgress";
import { useApp } from "../store/app";
import { BokehFigure } from "./BokehFigure";

interface LayerOptions { id: string; title: string; clusters: string[]; celltypes: { value: string; label: string }[] }
interface Options {
  layers: LayerOptions[];
  default_sample1: string | null;
  default_sample2: string | null;
  predicates: { value: string; label: string }[];
  predicate: string;
  min_coverage: number;
  coverage_step: number;
  ready: boolean;
}
interface Result {
  low: string;
  high: string;
  low_title: string;
  n_local: number;
  global_metrics: { metric: string; value: number | null; help: string }[];
  metrics: string[];
  metric_help: Record<string, string>;
  default_metric: string | null;
  default_x: string | null;
  default_y: string | null;
  gene_columns: { key: string; title: string; help: string }[];
  genes: (string | number | null)[][];
  prefix: string;
}

export const WORKSPACE = "workspace";        // the `layer` of the comparison's job events
const NONE = "";
const GENE_ROWS = 200;

/** Comparison panel: coherence of two samples at different resolutions. */
export function ComparisonPanel() {
  const job = useApp((s) => s.jobs[WORKSPACE]?.comparison);
  const version = useApp((s) => s.analysisVersion[WORKSPACE] ?? 0);
  const updateLayer = useApp((s) => s.updateLayer);
  const [options, setOptions] = useState<Options>();
  const [sample1, setSample1] = useState<string | null>(null);
  const [sample2, setSample2] = useState<string | null>(null);
  const [predicate, setPredicate] = useState<string | null>(null);
  const [coverage, setCoverage] = useState(0.8);
  const [keys, setKeys] = useState<Record<string, string>>({});
  const [result, setResult] = useState<Result>();
  const [tab, setTab] = useState<string | null>("global");
  const [helpRow, setHelpRow] = useState<number | null>(null);
  const [metric, setMetric] = useState<string | null>(null);
  const [x, setX] = useState<string | null>(null);
  const [y, setY] = useState<string | null>(null);
  const [plots, setPlots] = useState<any>(null);
  const [geneQuery, setGeneQuery] = useState("");
  const [geneSort, setGeneSort] = useState<{ column: number; desc: boolean }>({ column: 1, desc: true });
  const [message, setMessage] = useState<{ text: string; error: boolean } | null>(null);

  useEffect(() => {
    getJson<Options>("comparison").then((o) => {
      setOptions(o);
      setSample1((s) => s ?? o.default_sample1);
      setSample2((s) => s ?? o.default_sample2);
      setPredicate((p) => p ?? o.predicate);
      if (!options) setCoverage(o.min_coverage);
      if (o.ready) getJson<Result>("comparison/result").then((r) => {
        setResult(r);
        setMetric(r.default_metric);
        setX(r.default_x ?? r.metrics[0] ?? null);
        setY(r.default_y ?? r.metrics[1] ?? r.metrics[0] ?? null);
        setHelpRow(null);
      });
    });
  }, [version]);

  useEffect(() => {
    setPlots(null);
    if (!result || !metric || !x || !y) return;
    const params = new URLSearchParams({ metric, x, y });
    getJson(`comparison/plots?${params}`).then(setPlots).catch((exc) => setMessage({ text: exc.message, error: true }));
  }, [result, metric, x, y]);

  // Scatter selection -> the low-res layer's selection on the map
  useEffect(() => bridge.on("comparison-selection", async ({ indices }) => {
    if (!x || !y) return;
    const picked = await sendJson<{ layer: string; row_version: number; rows: number[] }>(
      "POST", "comparison/selection", { x, y, indices: Array.from(indices as number[]) });
    if (picked.row_version === bridge.rowVersion(picked.layer)) bridge.showSelection(picked.layer, picked.rows);
  }), [x, y]);

  const layerOptions = (id: string | null) => options?.layers.find((l) => l.id === id);
  const key = (name: string) => keys[name] || null;
  const compare = async () => {
    setMessage(null);
    try {
      await sendJson("POST", "comparison", {
        sample1, sample2, predicate, min_coverage: coverage,
        sample1_cluster_key: key("cluster1"), sample2_cluster_key: key("cluster2"),
        sample1_celltype_key: key("celltype1"), sample2_celltype_key: key("celltype2"),
      });
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    }
  };

  const show = async (layers: string[], column: string) => {
    try {
      for (const layer of layers) {
        await setStyle(layer, { colour_by: column });
        updateLayer(layer, { colour_by: { kind: "metadata", name: column } });
        requestViewport(layer);
      }
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    }
  };

  const genes = useMemo(() => {
    if (!result) return [];
    const q = geneQuery.trim().toLowerCase();
    const rows = q ? result.genes.filter((g) => String(g[0]).toLowerCase().includes(q)) : [...result.genes];
    const { column, desc } = geneSort;
    rows.sort((a, b) => {
      const p = a[column], r = b[column];
      if (p === null) return r === null ? 0 : 1;
      if (r === null) return -1;
      const order = typeof p === "string" ? p.localeCompare(String(r)) : (p as number) - (r as number);
      return desc ? -order : order;
    });
    return rows;
  }, [result, geneQuery, geneSort]);

  const keySelect = (name: string, label: string, data: { value: string; label: string }[]) => (
    <Select label={label} value={keys[name] ?? NONE} allowDeselect={false} size="xs"
            data={[{ value: NONE, label: "(none)" }, ...data]} onChange={(v) => setKeys((k) => ({ ...k, [name]: v ?? NONE }))} />
  );
  const running = job?.state === "running";
  const first = layerOptions(sample1), second = layerOptions(sample2);
  const sampleData = (options?.layers ?? []).map((l) => ({ value: l.id, label: l.title }));
  const format = (value: number | null) => (value === null ? "—" : value.toFixed(3));

  return (
    <Stack gap="sm">
      <Group grow>
        <Select label="Sample 1" data={sampleData} value={sample1} allowDeselect={false} size="xs"
                onChange={(v) => { setSample1(v); setKeys({}); }} />
        <Select label="Sample 2" data={sampleData} value={sample2} allowDeselect={false} size="xs"
                onChange={(v) => { setSample2(v); setKeys({}); }} />
      </Group>
      <Group grow align="end">
        <Select label="Match high-res geometries that are" data={options?.predicates ?? []} value={predicate}
                allowDeselect={false} size="xs" onChange={setPredicate} />
        <div>
          <Text size="xs">Min. coverage: {coverage.toFixed(2)}</Text>
          <Slider min={0} max={1} step={options?.coverage_step ?? 0.05} value={coverage} onChange={setCoverage}
                  label={(v) => v.toFixed(2)} />
        </div>
      </Group>
      <Group grow>
        {keySelect("cluster1", "Sample 1 clusters", (first?.clusters ?? []).map((c) => ({ value: c, label: c })))}
        {keySelect("cluster2", "Sample 2 clusters", (second?.clusters ?? []).map((c) => ({ value: c, label: c })))}
      </Group>
      <Group grow>
        {keySelect("celltype1", "Sample 1 cell types", first?.celltypes ?? [])}
        {keySelect("celltype2", "Sample 2 cell types", second?.celltypes ?? [])}
      </Group>
      <Button color="green" onClick={compare} loading={running} disabled={(options?.layers.length ?? 0) < 2}>
        Compare samples
      </Button>
      <JobProgress job={job} message={false} />
      <Text size="xs" c={job?.state === "failed" || message?.error ? "red.4" : "dark.1"}>
        {message?.text ?? (job && job.state !== "idle" ? job.message
          : (options?.layers.length ?? 0) < 2 ? "Load at least two samples to compare them."
          : "Choose the samples and press Compare samples.")}
      </Text>

      {result && (
        <Tabs value={tab} onChange={setTab}>
          <Tabs.List>
            <Tabs.Tab value="global">Global</Tabs.Tab>
            <Tabs.Tab value="local" disabled={!result.metrics.length}>Local</Tabs.Tab>
            <Tabs.Tab value="genes" disabled={!result.genes.length}>Genes</Tabs.Tab>
          </Tabs.List>
          <Tabs.Panel value="global" pt="xs">
            <ScrollArea.Autosize mah={380}>
              <Table fz="xs" highlightOnHover stickyHeader>
                <Table.Thead><Table.Tr><Table.Th>Metric</Table.Th><Table.Th style={{ textAlign: "right" }}>Value</Table.Th></Table.Tr></Table.Thead>
                <Table.Tbody>
                  {result.global_metrics.map((row, i) => (
                    <Table.Tr key={row.metric} onClick={() => setHelpRow(i)} style={{ cursor: "help" }}
                              bg={helpRow === i ? "dark.5" : undefined}>
                      <Table.Td title={row.help}><Text span c="blue.4">ⓘ</Text> {row.metric}</Table.Td>
                      <Table.Td style={{ textAlign: "right" }}>{format(row.value)}</Table.Td>
                    </Table.Tr>
                  ))}
                </Table.Tbody>
              </Table>
            </ScrollArea.Autosize>
            <Text size="xs" p={6} style={{ borderLeft: "3px solid #4e9af1" }}>
              {helpRow === null ? <i>ⓘ Click a metric to see what it captures.</i>
                : <><b>{result.global_metrics[helpRow].metric}</b>: {result.global_metrics[helpRow].help || "No description."}</>}
            </Text>
          </Tabs.Panel>
          <Tabs.Panel value="local" pt="xs">
            <Stack gap="xs">
              <Group grow align="end">
                <Select label="Local metric" size="xs" data={result.metrics} value={metric} onChange={setMetric} searchable />
                <Button size="xs" variant="default" disabled={!metric}
                        onClick={() => metric && show([result.low], result.prefix + metric)}>Show on map</Button>
              </Group>
              {metric && <Text size="xs" c="dark.1">{result.metric_help[metric] || "No description."}</Text>}
              <Group grow>
                <Select label="x" size="xs" data={result.metrics} value={x} onChange={setX} searchable />
                <Select label="y" size="xs" data={result.metrics} value={y} onChange={setY} searchable />
              </Group>
              <Text size="xs" c="dark.2" fs="italic">Select points to highlight the geometries on the map.</Text>
              <BokehFigure item={plots} docKey="comparison" />
            </Stack>
          </Tabs.Panel>
          <Tabs.Panel value="genes" pt="xs">
            <Stack gap="xs">
              <Text size="xs" c="dark.2" fs="italic">Click a gene to show it on both samples.</Text>
              <TextInput size="xs" placeholder="Search gene..." value={geneQuery} onChange={(e) => setGeneQuery(e.currentTarget.value)} />
              <ScrollArea.Autosize mah={380}>
                <Table fz="xs" highlightOnHover stickyHeader>
                  <Table.Thead>
                    <Table.Tr>
                      {[{ key: "gene", title: "Gene", help: "" }, ...result.gene_columns].map((column, i) => (
                        <Table.Th key={column.key} title={column.help} style={{ whiteSpace: "nowrap" }}>
                          <UnstyledButton onClick={() => setGeneSort((s) => ({ column: i, desc: s.column === i ? !s.desc : i > 0 }))}>
                            <Text size="xs" fw={600}>{column.title}{geneSort.column === i ? (geneSort.desc ? " ▼" : " ▲") : ""}</Text>
                          </UnstyledButton>
                        </Table.Th>
                      ))}
                    </Table.Tr>
                  </Table.Thead>
                  <Table.Tbody>
                    {genes.slice(0, GENE_ROWS).map((row) => (
                      <Table.Tr key={String(row[0])} style={{ cursor: "pointer" }}
                                onClick={() => show([result.low, result.high], String(row[0]))}>
                        <Table.Td>{row[0]}</Table.Td>
                        {row.slice(1).map((value, i) => <Table.Td key={i}>{format(value as number | null)}</Table.Td>)}
                      </Table.Tr>
                    ))}
                  </Table.Tbody>
                </Table>
              </ScrollArea.Autosize>
              <Text size="xs" c="dark.2">
                {genes.length > GENE_ROWS ? `First ${GENE_ROWS} of ${genes.length.toLocaleString()} genes · search to find any gene`
                  : `${genes.length.toLocaleString()} genes`}
              </Text>
            </Stack>
          </Tabs.Panel>
        </Tabs>
      )}
    </Stack>
  );
}
