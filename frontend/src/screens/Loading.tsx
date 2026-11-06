import { useEffect, useRef, useState } from "react";
import {
  Alert, Anchor, Button, Card, Checkbox, Container, FileButton, Group, Progress, Select, Stack, Tabs, Text, Title,
} from "@mantine/core";
import { getJson, sendJson, uploadFile } from "../api/client";
import { PluginParamForm, Param } from "../components/PluginParamForm";
import { useApp } from "../store/app";

interface Dataset { id: string; title: string; bins_size: number; has_reference: boolean }
interface Plugin { name: string; title: string; description: string; url?: string; params: Param[] }
interface Plugins { available: boolean; error: string | null; runner: string | null; segmentation: Plugin[] }
interface UploadState { status: string; filename: string; size: number; received: number; error: string | null; required: boolean }

const FIELDS: { key: string; label: string; accept: string; help: string }[] = [
  { key: "expression", label: "Expression matrix", accept: ".h5", help: "binned_outputs/square_XXXum/filtered_feature_bc_matrix.h5" },
  { key: "positions", label: "Bin positions", accept: ".parquet", help: "binned_outputs/square_XXXum/spatial/tissue_positions.parquet" },
  { key: "image", label: "Tissue image", accept: ".btf,.png,.tif,.tiff", help: "The full-resolution image given to Space Ranger" },
  { key: "reference", label: "Single-cell reference (optional)", accept: ".h5ad,.h5", help: "For deconvolution: AnnData with labels, or a 10x .h5" },
  { key: "reference_annotations", label: "Reference annotations (optional)", accept: ".csv,.parquet", help: "Barcode first, then cell-type labels" },
];

function bytes(size: number) {
  return size >= 1 << 30 ? `${(size / (1 << 30)).toFixed(1)} GB` : `${(size / (1 << 20)).toFixed(1)} MB`;
}

export function Loading() {
  const { config, load, failures, setLoad, setScreen } = useApp();
  const [source, setSource] = useState<string | null>("example");
  const [datasets, setDatasets] = useState<Dataset[]>([]);
  const [dataset, setDataset] = useState<string | null>(null);
  const [plugins, setPlugins] = useState<Plugins>();
  const [selected, setSelected] = useState<Record<string, Record<string, unknown> | undefined>>({});
  const [pseudospots, setPseudospots] = useState(true);
  const [binSize, setBinSize] = useState<string | null>(null);
  const [uploads, setUploads] = useState<Record<string, UploadState>>({});
  const [sending, setSending] = useState<Record<string, number>>({});
  const [error, setError] = useState<string | null>(null);
  const aborts = useRef<Record<string, AbortController>>({});
  const running = load.state === "running";

  useEffect(() => {
    getJson<Dataset[]>("datasets").then((list) => { setDatasets(list); setDataset(list[0]?.id ?? null); });
    getJson<Plugins>("plugins").then(setPlugins);
    getJson<Record<string, UploadState>>("uploads").then(setUploads);
  }, []);
  useEffect(() => { if (config) setBinSize(String(config.upload.default_bin_size)); }, [config]);

  const choose = async (key: string, file: File | null) => {
    if (!file) return;
    aborts.current[key]?.abort();
    const abort = (aborts.current[key] = new AbortController());
    setError(null);
    setSending((s) => ({ ...s, [key]: 0 }));
    try {
      await uploadFile(key, file, (sent) => setSending((s) => ({ ...s, [key]: sent / file.size })), abort.signal);
    } catch (exc: any) {
      if (exc.name !== "AbortError") setError(exc.message);
    } finally {
      setSending((s) => { const { [key]: _, ...rest } = s; return rest; });
      setUploads(await getJson("uploads"));
    }
  };

  const remove = async (key: string) => {
    aborts.current[key]?.abort();
    setUploads(await sendJson("DELETE", `uploads/${key}`));
  };

  const start = async () => {
    setError(null);
    const segmentation = Object.entries(selected).filter(([, v]) => v).map(([plugin, params]) => ({ plugin, params }));
    try {
      setLoad(await sendJson("POST", "load", {
        source, dataset_id: dataset, bins_size: Number(binSize), pseudospots, segmentation,
      }));
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  const required = FIELDS.filter((f) => uploads[f.key]?.required);
  const uploadsReady = required.every((f) => uploads[f.key]?.status === "complete") && !Object.keys(sending).length;
  const canLoad = !running && (source === "example" ? Boolean(dataset) : uploadsReady);

  return (
    <Container size="sm" py="xl">
      <Stack>
        <Group justify="space-between">
          <Title order={2}>Load a dataset</Title>
          <Button variant="subtle" onClick={() => setScreen("landing")}>Back</Button>
        </Group>
        <Card withBorder>
          <Tabs value={source} onChange={setSource}>
            <Tabs.List>
              <Tabs.Tab value="example" disabled={running}>Load example</Tabs.Tab>
              <Tabs.Tab value="upload" disabled={running}>Upload files</Tabs.Tab>
            </Tabs.List>
            <Tabs.Panel value="example" pt="md">
              <Select label="Example dataset" data={datasets.map((d) => ({ value: d.id, label: d.title }))}
                      value={dataset} onChange={setDataset} disabled={running} />
            </Tabs.Panel>
            <Tabs.Panel value="upload" pt="md">
              <Stack>
                {FIELDS.map((field) => {
                  const state = uploads[field.key];
                  const progress = sending[field.key];
                  return (
                    <Stack key={field.key} gap={4}>
                      <Group justify="space-between">
                        <div>
                          <Text fw={600} size="sm">{field.label}</Text>
                          <Text size="xs" c="dark.2">{field.help}</Text>
                        </div>
                        <Group gap="xs">
                          <FileButton accept={field.accept} onChange={(file) => choose(field.key, file)} disabled={running}>
                            {(props) => <Button size="xs" variant="default" {...props}>Choose file…</Button>}
                          </FileButton>
                          {(state?.status !== "empty" || progress !== undefined) &&
                            <Button size="xs" variant="subtle" color="red" onClick={() => remove(field.key)} disabled={running}>
                              {progress !== undefined ? "Cancel" : "Remove"}
                            </Button>}
                        </Group>
                      </Group>
                      {progress !== undefined && <Progress value={progress * 100} animated />}
                      {state?.status === "complete" && <Text size="xs" c="teal.3">Uploaded: {state.filename} ({bytes(state.size)})</Text>}
                      {state?.status === "error" && <Text size="xs" c="red.4">{state.error}</Text>}
                    </Stack>
                  );
                })}
                <Select label="Bin size of the uploaded files" value={binSize} onChange={setBinSize}
                        data={(config?.upload.bin_sizes ?? [2, 8, 16]).map((s) => ({ value: String(s), label: `${s} µm` }))} />
                <Text size="xs" c="dark.2">Up to {bytes(config?.upload.max_bytes ?? 0)} per file. Files are deleted when your session ends.</Text>
              </Stack>
            </Tabs.Panel>
          </Tabs>
        </Card>

        <Card withBorder>
          <Stack>
            <Title order={4}>Preprocessing</Title>
            <Checkbox label="Generate Visium v2 pseudospots (55 µm)" checked={pseudospots} disabled={running}
                      onChange={(e) => setPseudospots(e.currentTarget.checked)} />
            <Text fw={600} size="sm">Cell segmentation</Text>
            {plugins && !plugins.available && <Text size="sm" c="red.4">Segmentation is unavailable: {plugins.error}</Text>}
            {plugins?.segmentation.map((plugin) => (
              <Stack key={plugin.name} gap={6}>
                <Checkbox label={plugin.title} checked={Boolean(selected[plugin.name])} disabled={running}
                          onChange={(e) => setSelected((s) => ({ ...s, [plugin.name]: e.currentTarget.checked ? {} : undefined }))} />
                {selected[plugin.name] && (
                  <Stack pl="lg" gap="xs">
                    <Text size="xs" c="dark.2">{plugin.description} {plugin.url && <Anchor href={plugin.url} target="_blank" size="xs">More</Anchor>}</Text>
                    <PluginParamForm params={plugin.params} values={selected[plugin.name]!}
                                     onChange={(values) => setSelected((s) => ({ ...s, [plugin.name]: values }))} />
                  </Stack>
                )}
              </Stack>
            ))}
          </Stack>
        </Card>

        {error && <Alert color="red">{error}</Alert>}
        {load.state === "failed" && <Alert color="red">{load.error}</Alert>}
        {failures.map((message) => <Alert key={message} color="yellow">{message}</Alert>)}

        <Button size="md" onClick={start} disabled={!canLoad} loading={running}>Load data</Button>
        {load.tasks && (
          <Card withBorder>
            <Stack gap="xs">
              {Object.entries(load.tasks).map(([task, state]) => (
                <div key={task}>
                  <Group justify="space-between"><Text size="sm">{task}</Text><Text size="xs" c="dark.2">{state.status}</Text></Group>
                  <Progress value={state.pct} animated={running && state.pct < 100} />
                </div>
              ))}
            </Stack>
          </Card>
        )}
      </Stack>
    </Container>
  );
}
