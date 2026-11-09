import { create } from "zustand";
import { LoadStatus } from "../api/client";

export interface Config {
  branding: Record<string, string>;
  bokeh_version: string;
  tutorial_url: string;
  use_case_url: string;
  upload: { max_bytes: number; chunk_bytes: number; bin_sizes: number[]; default_bin_size: number };
  expression_scales: { value: string; label: string }[];
  roi: { default_radius: number; preview_debounce_ms: number };
}

export interface Layer {
  id: string;
  kind: "bins" | "spots" | "cells";
  title: string;
  visible: boolean;
  colour: string;
  alpha: number;
  colour_by: { kind: string; name?: string };
  scale: string;
  show_not_sequenced: boolean;
  aggregation_factors: number[];
  fixed_aggregation: number | null;
}

export interface Workspace {
  layers: Layer[];
  image: { size: number; inverse: number[] | null; um_per_px: number | null };
  tiles: { template: string; single: boolean };
  reference: string | null;
}

type Screen = "landing" | "loading" | "workspace";

/** The displayed bins (frame headers of bins layers; the map's overlay). */
export interface BinInfo {
  factor: number;
  width: number | null;
  height: number | null;
  rotation: number | null;
  bin_size_um: number | null;
  fixed: number | null;
  minimum: number | null;
  factors: number[];
}

export interface Job {
  name: string | null;
  state: "idle" | "running" | "done" | "failed";
  message: string;
  started?: number;      // ms since epoch, when the browser saw the run start
  log?: string[];        // recent step messages of the run (newest last)
}

const LOG_LINES = 6;
export type ResultKind = "de" | "metadata";

interface State {
  screen: Screen;
  config?: Config;
  load: LoadStatus;
  failures: string[];
  workspace?: Workspace;
  jobs: Record<string, Record<string, Job>>;   // layer -> job name -> last state
  analysisVersion: Record<string, number>;   // bumped when a layer's analyses change
  results: { layer: string; kind: ResultKind; run: number } | null;   // shown below the map
  showResults(layer: string | null, kind?: ResultKind): void;
  setJob(layer: string, job: Job): void;
  analysesChanged(layer: string): void;
  setScreen(screen: Screen): void;
  setConfig(config: Config): void;
  setLoad(load: LoadStatus): void;
  progress(task: string, pct: number, status: string): void;
  fail(message: string): void;
  setWorkspace(workspace: Workspace): void;
  updateLayer(id: string, changes: Partial<Layer>): void;
  activeLayer: string | null;               // the layer the tool panel works on
  setActiveLayer(id: string): void;
  binInfo: Record<string, BinInfo>;
  setBinInfo(layer: string, info: BinInfo): void;
}

export const useApp = create<State>((set) => ({
  screen: "landing",
  load: {},
  failures: [],
  jobs: {},
  analysisVersion: {},
  results: null,
  binInfo: {},
  setBinInfo: (layer, info) => set((state) => ({ binInfo: { ...state.binInfo, [layer]: info } })),
  showResults: (layer, kind = "de") =>
    set((state) => ({ results: layer ? { layer, kind, run: (state.results?.run ?? 0) + 1 } : null })),
  setJob: (layer, job) =>
    set((state) => {
      const previous = state.jobs[layer]?.[job.name ?? ""];
      const continuing = previous?.state === "running" && job.state === "running";
      const log = continuing ? [...(previous.log ?? [])] : [];
      if (job.message && log[log.length - 1] !== job.message) log.push(job.message);
      const entry = { ...job, started: continuing ? previous.started : Date.now(), log: log.slice(-LOG_LINES) };
      return { jobs: { ...state.jobs, [layer]: { ...state.jobs[layer], [job.name ?? ""]: entry } } };
    }),
  analysesChanged: (layer) =>
    set((state) => ({ analysisVersion: { ...state.analysisVersion, [layer]: (state.analysisVersion[layer] ?? 0) + 1 } })),
  setScreen: (screen) => set({ screen }),
  setConfig: (config) => set({ config }),
  setLoad: (load) => set({ load, failures: [] }),
  progress: (task, pct, status) =>
    set((state) => ({
      load: { ...state.load, tasks: { ...(state.load.tasks ?? {}), [task]: { pct, status } } },
    })),
  fail: (message) => set((state) => ({ failures: [...state.failures, message] })),
  setWorkspace: (workspace) => set((state) => ({
    workspace,
    activeLayer: workspace.layers.some((l) => l.id === state.activeLayer) ? state.activeLayer : workspace.layers[0]?.id ?? null,
  })),
  activeLayer: null,
  setActiveLayer: (id) => set({ activeLayer: id }),
  updateLayer: (id, changes) =>
    set((state) => state.workspace ? {
      workspace: { ...state.workspace, layers: state.workspace.layers.map((l) => (l.id === id ? { ...l, ...changes } : l)) },
    } : {}),
}));
