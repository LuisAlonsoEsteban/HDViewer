import type { ReactNode } from "react";
import type { Job } from "../store/app";

export interface RunInfo {
  id: string;
  method: string;
  cell_types: string[];
  fraction_columns: string[];
}

export interface RoiInfo {
  name: string;
  group: string | null;
  seed_group: string | null;
  group_excluded: boolean;
}

/** What the tour knows about the app, gathered from the store and the API. Step
 * checks are pure functions of these facts (so they can be unit tested). */
export interface Facts {
  screen: "landing" | "loading" | "workspace";
  loadState: "running" | "ready" | "failed" | undefined;
  dataset: string | null | undefined;        // the loaded example's id
  datasets: string[];                        // examples offered by this server
  hasWorkspace: boolean;
  binsLayer: string | null;
  activeKind: string | null;
  colourBy: string | null;                   // what the bins layer is coloured by
  deconvJob: Job | undefined;
  methods: { name: string; unavailable: string }[];
  runs: RunInfo[];
  groups: string[];
  rois: RoiInfo[];
  results: "de" | "metadata" | null;
  now: number;
}

export interface TourStep {
  id: string;
  title: string;
  body(facts: Facts): ReactNode;
  /** `data-tour` keys to highlight; their elements are outlined together. */
  anchor?: string[];
  /** The screen the step happens on (the card offers a way back if elsewhere). */
  screen?: Facts["screen"];
  /** The workspace tool to open when the step starts. */
  tool?: string;
  done?(facts: Facts): boolean;
  /** Section of the biological use-case page with the full explanation. */
  readMore?: string;
}

export interface Tour {
  id: string;
  title: string;
  dataset: string;                           // the example dataset the tour needs
  steps: TourStep[];
}
