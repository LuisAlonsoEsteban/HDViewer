import { useEffect, useRef, useState } from "react";
import { Button, Checkbox, MultiSelect, Select, Stack, Text } from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { Layer, useApp } from "../store/app";
import { DESet, SetList } from "./DEPanel";

interface Preset {
  key: string;
  title: string;
  all_columns: boolean;           // a deconvolution readout: every cell type by default
  columns: { name: string; label: string }[];
}
interface Options {
  groups: DESet[];
  clusters: DESet[];
  presets: Preset[];
  default_preset: string;
  hint: string;
  ready: boolean;
}

/** The Bokeh UI's Metadata comparison tab: summaries of numeric metadata (e.g.
 * deconvolution fractions) per group/cluster; the charts open below the map. */
export function MetadataPanel({ layer }: { layer: Layer }) {
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const showResults = useApp((s) => s.showResults);
  const shown = useApp((s) => s.results?.layer === layer.id && s.results.kind === "metadata");
  const [options, setOptions] = useState<Options>();
  const [groups, setGroups] = useState<Set<string>>(new Set());
  const [clusters, setClusters] = useState<Set<string>>(new Set());
  const [rest, setRest] = useState(false);
  const [preset, setPreset] = useState<string | null>(null);
  const [columns, setColumns] = useState<string[]>([]);
  const [statistic, setStatistic] = useState("mean");
  const [busy, setBusy] = useState(false);
  const [message, setMessage] = useState<{ text: string; error: boolean } | null>(null);
  const base = `layers/${encodeURIComponent(layer.id)}/metadata`;
  const chosenPreset = useRef<string | null>(null);
  const chosenColumns = useRef<string[]>([]);
  useEffect(() => { chosenColumns.current = columns; }, [columns]);

  const choosePreset = (key: string, presets: Preset[], keep: string[] = []) => {
    const spec = presets.find((p) => p.key === key);
    if (!spec) return;
    const names = spec.columns.map((c) => c.name);
    const kept = keep.filter((name) => names.includes(name));
    chosenPreset.current = key;
    setPreset(key);
    setColumns(kept.length ? kept : spec.all_columns ? names : []);
  };

  useEffect(() => {
    getJson<Options>(`${base}/options`).then((o) => {
      setOptions(o);
      const existing = (sets: DESet[]) => (chosen: Set<string>) =>
        new Set(sets.map((s) => s.name).filter((n) => chosen.has(n)));
      setGroups(existing(o.groups));
      setClusters(existing(o.clusters));
      // As the Bokeh tab: the latest run's set replaces Numeric metadata while no column is chosen
      const current = chosenPreset.current;
      const keep = current !== null && o.presets.some((p) => p.key === current)
        && (current !== "numeric" || chosenColumns.current.length > 0);
      if (!keep && o.default_preset !== current) choosePreset(o.default_preset, o.presets);
    });
  }, [version, base]);

  const toggle = (setter: typeof setGroups) => (name: string, on: boolean) =>
    setter((chosen) => {
      const next = new Set(chosen);
      if (on) next.add(name);
      else next.delete(name);
      return next;
    });

  const compare = async () => {
    const chosen = groups.size + clusters.size;
    if (chosen < (rest ? 1 : 2) || !columns.length) {
      const requirement = rest ? "at least one group or cluster with The rest" : "at least two groups or clusters";
      return setMessage({ text: `Select ${requirement} and one or more numeric metadata columns.`, error: true });
    }
    const ordered = (sets: DESet[], picked: Set<string>) => sets.map((s) => s.name).filter((n) => picked.has(n));
    setBusy(true);
    setMessage({ text: "Comparing metadata…", error: false });
    try {
      await sendJson("POST", `${base}/compare`, {
        groups: ordered(options!.groups, groups), clusters: ordered(options!.clusters, clusters),
        include_rest: rest, preset, columns, statistic,
      });
      setMessage({ text: "Metadata comparison ready below the map.", error: false });
      showResults(layer.id, "metadata");
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    } finally {
      setBusy(false);
    }
  };

  const spec = options?.presets.find((p) => p.key === preset);
  return (
    <Stack gap="sm">
      <Text size="xs" c="dark.1">Compare selected groups/clusters, or compare each selected set against The rest.</Text>
      <SetList title="Groups" sets={options?.groups ?? []} chosen={groups} toggle={toggle(setGroups)}
               hint="Create groups in the Groups tab (or from an ROI)." />
      <SetList title="Clusters" sets={options?.clusters ?? []} chosen={clusters} toggle={toggle(setClusters)}
               hint="Run clustering in the Clustering tab." />
      <Checkbox label="The rest (outside all defined custom groups and the compared set)" checked={rest}
                onChange={(e) => setRest(e.currentTarget.checked)} />
      <Select label="Metadata set" value={preset} allowDeselect={false}
              data={(options?.presets ?? []).map((p) => ({ value: p.key, label: p.title }))}
              onChange={(key) => key && options && choosePreset(key, options.presets)} />
      <MultiSelect label="Metadata columns" searchable clearable value={columns} onChange={setColumns}
                   disabled={!spec?.columns.length} placeholder={columns.length ? undefined : "Choose columns"}
                   data={(spec?.columns ?? []).map((c) => ({ value: c.name, label: c.label }))} />
      <Select label="Summary" value={statistic} allowDeselect={false} onChange={(value) => setStatistic(value ?? "mean")}
              data={[{ value: "mean", label: "Mean" }, { value: "median", label: "Median" }]} />
      {options && (
        <Text size="xs" c="dark.2">
          {options.hint} Select The rest for one plot per selected group/cluster against the observations outside
          all defined custom groups and outside that group/cluster.
        </Text>
      )}
      <Button onClick={compare} loading={busy} disabled={!options}>Compare metadata</Button>
      {message && <Text size="xs" c={message.error ? "red.4" : "dark.1"}>{message.text}</Text>}
      {options?.ready && !shown && (
        <Button size="xs" variant="light" onClick={() => showResults(layer.id, "metadata")}>Show metadata comparison</Button>
      )}
    </Stack>
  );
}
