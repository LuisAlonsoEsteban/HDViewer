import { useEffect, useState } from "react";
import { Button, Checkbox, ColorSwatch, Group, NumberInput, Select, Slider, Stack, Text } from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { patchFrame, requestViewport } from "../map/layers";
import { Layer, useApp } from "../store/app";

interface Options {
  algorithms: string[];
  defaults: Record<string, any>;
  resolution: { min: number; max: number; step: number; value: number };
}
interface Cluster { label: string; colour: string; size: number; visible: boolean; markers: string[] }
interface Clusters { column: string | null; clusters: Cluster[]; summary: string; shown: boolean }

const LABELS: Record<string, string> = {
  leiden: "Leiden (graph)", kmeans: "k-means", dbscan: "DBSCAN", hdbscan: "HDBSCAN",
  gaussian_mixture: "Gaussian mixture", agglomerative: "Agglomerative",
};

/** The Bokeh UI's Clustering tab: run an algorithm, then show/hide clusters. */
export function ClusteringPanel({ layer }: { layer: Layer }) {
  const job = useApp((s) => s.jobs[layer.id]?.clustering);
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const [options, setOptions] = useState<Options>();
  const [params, setParams] = useState<Record<string, any>>({});
  const [result, setResult] = useState<Clusters>();
  const [error, setError] = useState<string | null>(null);
  const path = `layers/${encodeURIComponent(layer.id)}/clusters`;

  useEffect(() => {
    getJson<Options>("clustering/options").then((o) => {
      setOptions(o);
      setParams({ ...o.defaults, resolution: o.resolution.value });
    });
  }, []);
  useEffect(() => { getJson<Clusters>(path).then(setResult); }, [version, path]);

  const run = async () => {
    setError(null);
    const { algorithm, resolution, n_clusters, eps, min_samples, min_cluster_size, covariance_type, linkage } = params;
    try {
      await sendJson("POST", `layers/${encodeURIComponent(layer.id)}/clustering`,
                     { algorithm, resolution, n_clusters, eps, min_samples, min_cluster_size, covariance_type, linkage });
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  const display = async (changes: Record<string, unknown>) => {
    if (Array.isArray(changes.visible)) {
      const shown = new Set(changes.visible as string[]);   // show the change at once
      setResult((r) => r && { ...r, clusters: r.clusters.map((c) => ({ ...c, visible: shown.has(c.label) })) });
    }
    await patchFrame(path, changes);
    setResult(await getJson<Clusters>(path));
  };
  const visible = result?.clusters.filter((c) => c.visible).map((c) => c.label) ?? [];
  const toggle = (label: string, on: boolean) =>
    display({ visible: on ? [...visible, label] : visible.filter((l) => l !== label) });

  const algorithm = params.algorithm ?? "leiden";
  const set = (key: string) => (value: unknown) => setParams((p) => ({ ...p, [key]: value }));
  const running = job?.state === "running";

  return (
    <Stack gap="sm">
      <Select label="Algorithm" value={algorithm} onChange={set("algorithm")}
              data={(options?.algorithms ?? []).map((a) => ({ value: a, label: LABELS[a] ?? a }))} />
      {algorithm === "leiden" && options && (
        <div>
          <Text size="sm">Resolution: {params.resolution}</Text>
          <Slider min={options.resolution.min} max={options.resolution.max} step={options.resolution.step}
                  value={params.resolution} onChange={set("resolution")} />
        </div>
      )}
      {["kmeans", "gaussian_mixture", "agglomerative"].includes(algorithm) &&
        <NumberInput label="Number of clusters" min={2} value={params.n_clusters} onChange={set("n_clusters")} />}
      {algorithm === "dbscan" && <NumberInput label="Neighbourhood radius (eps)" min={0} step={0.1} decimalScale={3}
                                              value={params.eps} onChange={set("eps")} />}
      {algorithm === "hdbscan" && <NumberInput label="Minimum cluster size" min={2} value={params.min_cluster_size}
                                               onChange={set("min_cluster_size")} />}
      {["dbscan", "hdbscan"].includes(algorithm) &&
        <NumberInput label="Minimum samples" min={1} value={params.min_samples} onChange={set("min_samples")} />}
      {algorithm === "gaussian_mixture" &&
        <Select label="Covariance" value={params.covariance_type} onChange={set("covariance_type")}
                data={["full", "tied", "diag", "spherical"]} />}
      {algorithm === "agglomerative" &&
        <Select label="Linkage" value={params.linkage} onChange={set("linkage")} data={["ward", "complete", "average", "single"]} />}
      <Text size="xs" c="dark.2">Clustering uses expression PCA coordinates. Density methods label outliers as Noise.</Text>
      <Button onClick={run} loading={running} disabled={!options}>Run clustering</Button>
      {job && job.state !== "idle" && <Text size="xs" c={job.state === "failed" ? "red.4" : "dark.1"}>{job.message}</Text>}
      {error && <Text size="xs" c="red.4">{error}</Text>}

      {result && result.clusters.length > 0 && (
        <Stack gap={6}>
          <Group justify="space-between">
            <Text fw={600} size="sm">{result.summary}</Text>
            <Button size="compact-xs" variant="default"
                    onClick={() => display({ visible: visible.length === result.clusters.length ? [] : result.clusters.map((c) => c.label) })}>
              Select all / none
            </Button>
          </Group>
          {!result.shown && <Button size="xs" variant="light" onClick={async () => { await display({ show: true }); requestViewport(layer.id); }}>
            Colour the map by clusters</Button>}
          {result.clusters.map((cluster) => (
            <div key={cluster.label}>
              <Checkbox checked={cluster.visible} onChange={(e) => toggle(cluster.label, e.currentTarget.checked)}
                        label={<Group gap={6} wrap="nowrap"><ColorSwatch color={cluster.colour} size={12} />
                          <Text size="sm" c={cluster.colour} fw={600}>
                            {cluster.label === "Noise" ? "Noise" : `Cluster ${cluster.label}`}
                          </Text><Text size="xs" c="dark.2">({cluster.size} {layer.kind})</Text></Group>} />
              <Text size="xs" c="dark.2" pl={30}>{cluster.markers.length ? cluster.markers.join(", ") : running ? "computing marker genes…" : "–"}</Text>
            </div>
          ))}
        </Stack>
      )}
    </Stack>
  );
}
