import { useEffect, useState } from "react";
import { getJson } from "../api/client";
import { useApp } from "../store/app";
import type { Facts, RoiInfo, RunInfo } from "./types";

interface Remote {
  methods: Facts["methods"];
  runs: RunInfo[];
  groups: string[];
  rois: RoiInfo[];
}
const EMPTY: Remote = { methods: [], runs: [], groups: [], rois: [] };

let datasetsRequest: Promise<string[]> | null = null;
/** Example datasets offered by the server (asked once per page). */
export function exampleDatasets(): Promise<string[]> {
  datasetsRequest ??= getJson<{ id: string }[]>("datasets").then((list) => list.map((d) => d.id))
    .catch(() => { datasetsRequest = null; return []; });
  return datasetsRequest;
}

/** The facts the tour's step checks read. The server is asked again only when the
 * bins layer's analyses change (a group, ROI or deconvolution run), never on a timer. */
export function useTourFacts(active: boolean, step: number): Facts {
  const screen = useApp((s) => s.screen);
  const load = useApp((s) => s.load);
  const workspace = useApp((s) => s.workspace);
  const activeLayer = useApp((s) => s.activeLayer);
  const results = useApp((s) => s.results?.kind ?? null);
  const bins = workspace?.layers.find((l) => l.kind === "bins") ?? null;
  const version = useApp((s) => (bins ? s.analysisVersion[bins.id] ?? 0 : 0));
  const deconvJob = useApp((s) => (bins ? s.jobs[bins.id]?.deconvolution : undefined));
  const [datasets, setDatasets] = useState<string[]>([]);
  const [remote, setRemote] = useState<Remote>(EMPTY);
  const [now, setNow] = useState(() => Date.now());

  useEffect(() => { if (active) exampleDatasets().then(setDatasets); }, [active]);

  useEffect(() => {
    if (!active || screen !== "workspace" || !bins) return;
    let stale = false;
    const base = `layers/${encodeURIComponent(bins.id)}`;
    Promise.all([
      getJson<{ methods: Remote["methods"]; runs: RunInfo[] }>(`${base}/deconvolution`).catch(() => null),
      getJson<{ groups: { name: string }[] }>(`${base}/groups`).catch(() => null),
      getJson<{ rois: RoiInfo[] }>(`${base}/rois`).catch(() => null),
    ]).then(([deconv, groups, rois]) => {
      if (stale) return;
      setRemote({
        methods: deconv?.methods.map(({ name, unavailable }) => ({ name, unavailable })) ?? [],
        runs: deconv?.runs ?? [],
        groups: groups?.groups.map((g) => g.name) ?? [],
        rois: rois?.rois ?? [],
      });
    });
    return () => { stale = true; };
  }, [active, screen, bins?.id, version, deconvJob?.state, step]);

  // A clock for "this is taking long" hints, only while a deconvolution runs
  useEffect(() => {
    if (!active || deconvJob?.state !== "running") return;
    const timer = setInterval(() => setNow(Date.now()), 5000);
    return () => clearInterval(timer);
  }, [active, deconvJob?.state]);

  const activeKind = workspace?.layers.find((l) => l.id === activeLayer)?.kind ?? null;
  return {
    screen, loadState: load.state, dataset: load.dataset, datasets,
    hasWorkspace: Boolean(workspace), binsLayer: bins?.id ?? null, activeKind,
    colourBy: bins?.colour_by.name ?? null, deconvJob, ...remote, results, now,
  };
}
