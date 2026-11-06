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

export interface Job { name: string | null; state: "idle" | "running" | "done" | "failed"; message: string }
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
}

export const useApp = create<State>((set) => ({
  screen: "landing",
  load: {},
  failures: [],
  jobs: {},
  analysisVersion: {},
  results: null,
  showResults: (layer, kind = "de") =>
    set((state) => ({ results: layer ? { layer, kind, run: (state.results?.run ?? 0) + 1 } : null })),
  setJob: (layer, job) =>
    set((state) => ({ jobs: { ...state.jobs, [layer]: { ...state.jobs[layer], [job.name ?? ""]: job } } })),
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
  setWorkspace: (workspace) => set({ workspace }),
  updateLayer: (id, changes) =>
    set((state) => state.workspace ? {
      workspace: { ...state.workspace, layers: state.workspace.layers.map((l) => (l.id === id ? { ...l, ...changes } : l)) },
    } : {}),
}));
