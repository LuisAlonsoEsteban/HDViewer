import { useEffect, useState } from "react";
import { Button, Checkbox, ColorSwatch, Group, SegmentedControl, Select, Stack, Text } from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { Layer, useApp } from "../store/app";
import { JobProgress } from "./JobProgress";

export interface DESet { name: string; colour: string | null }
interface Options {
  groups: DESet[];
  clusters: DESet[];
  factors: number[];
  displayed_factor: number;
  ready: boolean;
  header: string | null;
}

export function SetList({ title, sets, chosen, toggle, hint }:
  { title: string; sets: DESet[]; chosen: Set<string>; toggle(name: string, on: boolean): void; hint: string }) {
  return (
    <Stack gap={4}>
      <Text size="xs" c="dark.2">{title}</Text>
      {sets.length === 0 && <Text size="xs" c="dark.2" fs="italic">{hint}</Text>}
      {sets.map((set) => (
        <Checkbox key={set.name} checked={chosen.has(set.name)} onChange={(e) => toggle(set.name, e.currentTarget.checked)}
                  label={<Group gap={6} wrap="nowrap">
                    {set.colour && <ColorSwatch color={set.colour} size={12} />}
                    <Text size="sm">{set.name}</Text>
                  </Group>} />
      ))}
    </Stack>
  );
}

/** DE tool: compare groups and/or clusters; results open below the map. */
export function DEPanel({ layer }: { layer: Layer }) {
  const job = useApp((s) => s.jobs[layer.id]?.de);
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const showResults = useApp((s) => s.showResults);
  const shown = useApp((s) => s.results?.layer === layer.id && s.results.kind === "de");
  const [options, setOptions] = useState<Options>();
  const [groups, setGroups] = useState<Set<string>>(new Set());
  const [clusters, setClusters] = useState<Set<string>>(new Set());
  const [mode, setMode] = useState("0");
  const [level, setLevel] = useState("displayed");
  const [error, setError] = useState<string | null>(null);
  const base = `layers/${encodeURIComponent(layer.id)}`;

  const load = () => getJson<Options>(`${base}/de/options`).then((o) => {
    setOptions(o);
    // Keep the choices that still exist
    const keep = (sets: DESet[]) => (chosen: Set<string>) => new Set(sets.map((s) => s.name).filter((n) => chosen.has(n)));
    setGroups(keep(o.groups));
    setClusters(keep(o.clusters));
  });
  useEffect(() => { load(); }, [version, base]);

  const toggle = (setter: typeof setGroups) => (name: string, on: boolean) =>
    setter((chosen) => {
      const next = new Set(chosen);
      if (on) next.add(name);
      else next.delete(name);
      return next;
    });

  const run = async () => {
    setError(null);
    const chosen = groups.size + clusters.size;
    if (!chosen) return setError("Select at least one group or cluster to compare.");
    if (mode === "1" && chosen < 2) return setError('"Between selected" needs at least 2 selected groups or clusters.');
    const ordered = (sets: DESet[], picked: Set<string>) => sets.map((s) => s.name).filter((n) => picked.has(n));
    try {
      await sendJson("POST", `${base}/de`, {
        groups: ordered(options!.groups, groups), clusters: ordered(options!.clusters, clusters), mode: Number(mode),
        factor: level === "displayed" ? null : Number(level),
      });
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  const running = job?.state === "running";
  const displayed = options?.displayed_factor ?? 1;
  const factors = options?.factors ?? [];

  return (
    <Stack gap="sm">
      <Text fw={600} size="sm">Compare</Text>
      <SetList title="Groups" sets={options?.groups ?? []} chosen={groups} toggle={toggle(setGroups)}
               hint="Create groups in the Groups tab (or from an ROI)." />
      <SetList title="Clusters" sets={options?.clusters ?? []} chosen={clusters} toggle={toggle(setClusters)}
               hint="Run clustering in the Clustering tab." />
      <Text fw={600} size="sm">Against</Text>
      <SegmentedControl value={mode} onChange={setMode}
                        data={[{ value: "0", label: "vs Rest" }, { value: "1", label: "Between selected" }]} />
      {factors.length > 0 && (
        <Select label="Resolution" value={level} allowDeselect={false} onDropdownOpen={load} onChange={(value) => setLevel(value ?? "displayed")}
                data={[{ value: "displayed", label: `Displayed level (${displayed === 1 ? "native" : `${displayed} × ${displayed}`})` },
                       { value: "1", label: "Native bins" },
                       ...factors.filter((f) => f > 1).map((f) => ({ value: String(f), label: `${f} × ${f} bins (summed counts)` }))]} />
      )}
      <Button onClick={run} loading={running} disabled={!options}>Run DE</Button>
      <JobProgress job={job} />
      {job && job.state !== "idle" && job.state !== "running" && <Text size="xs" c={job.state === "failed" ? "red.4" : "dark.1"}>{job.message}</Text>}
      {error && <Text size="xs" c="red.4">{error}</Text>}
      {options?.ready && !shown && (
        <Button size="xs" variant="light" onClick={() => showResults(layer.id, "de")}>Show results: {options.header}</Button>
      )}
    </Stack>
  );
}
