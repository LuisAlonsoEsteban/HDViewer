import { useEffect, useRef, useState } from "react";
import { Anchor, Button, Select, Stack, Text } from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { requestViewport, setStyle } from "../map/layers";
import { Layer, useApp } from "../store/app";
import { JobProgress } from "./JobProgress";
import { Param, PluginParamForm } from "./PluginParamForm";

interface Method {
  name: string;
  title: string;
  description: string;
  url: string | null;
  citation: string | null;
  params: Param[];
  unavailable: string;
}
interface Run {
  id: string;
  method: string;
  cell_types: string[];
  fraction_columns: string[];
  abundance_columns: string[];
  abundance_units: string;
}
interface State {
  unavailable: string;
  reference: { name: string; cells: number; genes: number } | null;
  columns: string[];
  label_column: string;
  runner: string | null;
  methods: Method[];
  runs: Run[];
}

/** Deconvolution tool: run a method against the single-cell
 * reference, then colour the layer by a cell type of a saved run. */
export function DeconvolutionPanel({ layer }: { layer: Layer }) {
  const job = useApp((s) => s.jobs[layer.id]?.deconvolution);
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const updateLayer = useApp((s) => s.updateLayer);
  const [state, setState] = useState<State>();
  const [label, setLabel] = useState<string | null>(null);
  const [method, setMethod] = useState<string | null>(null);
  const [values, setValues] = useState<Record<string, Record<string, unknown>>>({});   // per method
  const [reason, setReason] = useState("");
  const [run, setRun] = useState<string | null>(null);
  const [readout, setReadout] = useState("fractions");
  const [display, setDisplay] = useState<string | null>(layer.colour_by.name ?? null);
  const seenRuns = useRef<number | null>(null);
  const [error, setError] = useState<string | null>(null);
  const base = `layers/${encodeURIComponent(layer.id)}/deconvolution`;
  const running = job?.state === "running";

  useEffect(() => {
    getJson<State>(base).then((s) => {
      setState(s);
      setLabel((current) => (current && s.columns.includes(current) ? current : s.label_column || null));
      setMethod((current) => current ?? s.methods.find((m) => m.name === "nnls")?.name ?? s.methods[0]?.name ?? null);
      // A run that just finished is chosen and its first cell type shown (as in the Bokeh tab)
      const fresh = seenRuns.current !== null && s.runs.length > seenRuns.current;
      seenRuns.current = s.runs.length;
      const newest = s.runs[s.runs.length - 1];
      if (fresh) {
        setRun(newest.id);
        setReadout("fractions");
        show(newest.fraction_columns[0] ?? null);
      } else {
        setRun((current) => (s.runs.some((r) => r.id === current) ? current : newest?.id ?? null));
      }
    });
  }, [version, base]);

  const plugin = state?.methods.find((m) => m.name === method);
  const params = (method && values[method]) || {};

  // Live availability, as the Bokeh tab's hint under the form (debounced validation)
  useEffect(() => {
    if (!state) return;
    if (state.unavailable) return setReason(state.unavailable);
    if (!plugin) return setReason("");
    if (plugin.unavailable) return setReason(plugin.unavailable);
    const timer = setTimeout(() => {
      sendJson("POST", `${base}/check`, { method, label_column: label ?? "", params })
        .then(() => setReason("")).catch((exc) => setReason(exc.message));
    }, 250);
    return () => clearTimeout(timer);
  }, [state, method, label, values]);

  const start = async () => {
    setError(null);
    try {
      await sendJson("POST", base, { method, label_column: label, params });
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  const record = state?.runs.find((r) => r.id === run);
  const columns = record ? (readout === "abundance" ? record.abundance_columns : record.fraction_columns) : [];
  const choices = columns.map((column, i) => ({ value: column, label: record!.cell_types[i] ?? column }));

  const show = async (column: string | null) => {
    setDisplay(column);
    if (!column) return;
    setError(null);
    try {
      await setStyle(layer.id, { colour_by: column });
      updateLayer(layer.id, { colour_by: { kind: "metadata", name: column } });
      requestViewport(layer.id);
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  const ready = state?.reference && !state.unavailable;
  return (
    <Stack gap="sm">
      <Text size="xs" c="dark.1">
        Run a method on all native observations of this layer. Use linear counts and matching gene identifiers.
        Empty observations remain unassigned.
      </Text>
      <Select label="Reference cell-type column (.obs)" data={state?.columns ?? []} value={label} onChange={setLabel}
              disabled={!state?.columns.length || running} allowDeselect={false} />
      <Select label="Deconvolution method" value={method} allowDeselect={false} disabled={!state?.methods.length || running}
              data={(state?.methods ?? []).map((m) => ({ value: m.name, label: m.title }))} onChange={setMethod} />
      {plugin && (
        <Text size="xs" c="dark.1">
          {plugin.description}{" "}
          {plugin.url && <Anchor size="xs" href={plugin.url} target="_blank">Documentation</Anchor>}
          {plugin.citation && <><br />{plugin.citation}</>}
        </Text>
      )}
      {plugin && (
        <PluginParamForm params={plugin.params} values={params} columns={state?.columns ?? []} disabled={running}
                         onChange={(next) => setValues((v) => ({ ...v, [plugin.name]: next }))} />
      )}
      {reason && <Text size="xs" c="yellow.5">{reason}</Text>}
      <Button color="green" onClick={start} loading={running} disabled={!!reason || !plugin || !label}>
        Run deconvolution
      </Button>
      <JobProgress job={job} message={false} />
      <Text size="xs" c={job?.state === "failed" || error ? "red.4" : "dark.1"} style={{ whiteSpace: "pre-wrap" }}>
        {error ?? (job && job.state !== "idle" ? job.message
          : ready ? `Reference ready: ${state!.reference!.name} (${state!.reference!.cells.toLocaleString()} cells, `
                    + `${state!.reference!.genes.toLocaleString()} genes). Runner: ${state!.runner ?? "unavailable"}.`
          : "")}
      </Text>

      <Select label="Saved deconvolution run" disabled={!state?.runs.length} value={run} allowDeselect={false}
              data={(state?.runs ?? []).map((r) => ({ value: r.id, label: `Run ${r.id}: ${r.method}` }))}
              onChange={(value) => {
                setRun(value);
                setReadout("fractions");
                show(state?.runs.find((r) => r.id === value)?.fraction_columns[0] ?? null);
              }} />
      <Select label="Readout" disabled={!record} value={readout} allowDeselect={false}
              data={[{ value: "fractions", label: "Relative fractions" },
                     ...(record?.abundance_columns.length ? [{ value: "abundance", label: `Abundance: ${record.abundance_units}` }] : [])]}
              onChange={(value) => {
                setReadout(value ?? "fractions");
                show((value === "abundance" ? record?.abundance_columns : record?.fraction_columns)?.[0] ?? null);
              }} />
      <Select label="Display cell type" placeholder={choices.length ? "Choose a cell type" : undefined}
              disabled={!choices.length} data={choices} allowDeselect={false} searchable
              value={choices.some((c) => c.value === display) ? display : null} onChange={show} />
    </Stack>
  );
}
