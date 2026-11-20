import { useEffect, useState } from "react";
import { useApp } from "../store/app";
import { exampleDatasets } from "./facts";
import { TOURS } from "./p1crc";
import { useTour } from "./store";

/** Starts a tour; with its dataset already loaded, it begins in the workspace. */
export function startTour(id: string) {
  const tour = TOURS[id];
  if (!tour) return;
  const app = useApp.getState();
  const { start, go } = useTour.getState();
  start(id);
  if (app.workspace && app.load.state === "ready" && app.load.dataset === tour.dataset) {
    app.setScreen("workspace");
    go(Math.max(0, tour.steps.findIndex((s) => s.screen === "workspace")));
  }
}

/** Whether this server offers the tour's example dataset. */
export function useTourAvailable(id: string): boolean {
  const [available, setAvailable] = useState(false);
  useEffect(() => {
    const dataset = TOURS[id]?.dataset;
    exampleDatasets().then((ids) => setAvailable(Boolean(dataset && ids.includes(dataset))));
  }, [id]);
  return available;
}

/** `?tour=<id>` (from the use-case page) starts that tour, then leaves the URL clean. */
export function startTourFromUrl() {
  const url = new URL(window.location.href);
  const id = url.searchParams.get("tour");
  if (!id) return;
  url.searchParams.delete("tour");
  window.history.replaceState(null, "", url.toString());
  if (TOURS[id]) startTour(id);
}
