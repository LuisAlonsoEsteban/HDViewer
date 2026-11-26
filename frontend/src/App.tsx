import { useEffect, useState } from "react";
import { ApiError, connectEvents, getJson, LoadStatus, onCapacity } from "./api/client";
import { CapacityNotice } from "./components/CapacityNotice";
import { CookieNotice } from "./components/CookieNotice";
import { refreshLayer } from "./map/layers";
import { ErrorNotices, notifyError, reportError } from "./notify";
import { Config, useApp, Workspace as WorkspaceData } from "./store/app";
import { Landing } from "./screens/Landing";
import { Loading } from "./screens/Loading";
import { Workspace } from "./screens/Workspace";
import { startTourFromUrl } from "./tour/start";
import { TourOverlay } from "./tour/TourOverlay";

export function App() {
  const { screen, setConfig, setLoad, progress, fail, setWorkspace, setScreen } = useApp();
  // Seconds until the next try when every analysis session on the server is in use
  const [capacity, setCapacity] = useState<number | null>(null);
  useEffect(() => onCapacity(setCapacity), []);
  useEffect(startTourFromUrl, []);

  /** A free session? Then start again cleanly (the page sets up its session). */
  const retrySession = async () => {
    try {
      await getJson<LoadStatus>("load");
      window.location.reload();
      return true;
    } catch (exc) {
      if (!(exc instanceof ApiError && exc.reason === "capacity")) setCapacity(null);
      return false;
    }
  };

  useEffect(() => {
    const openWorkspace = async () => {
      try {
        setWorkspace(await getJson<WorkspaceData>("workspace"));
        setScreen("workspace");
      } catch (exc) {
        notifyError(exc, "The workspace could not be opened");
      }
    };
    getJson<Config>("config").then((config) => {
      setConfig(config);
      if (window.Bokeh?.version && window.Bokeh.version !== config.bokeh_version) {
        console.warn(`BokehJS ${window.Bokeh.version} differs from the server's Bokeh ${config.bokeh_version}`);
      }
    }).catch(reportError("The server could not be reached. Reload the page in a moment."));
    let stop = () => {};
    // GET /api/load creates the session cookie the events socket needs
    getJson<LoadStatus>("load").catch((exc) => {
      if (!(exc instanceof ApiError && exc.reason === "capacity")) {
        notifyError(exc, "The server could not be reached. Reload the page in a moment.");
      }
      return null;
    }).then((load) => {
      if (!load) return;          // no session yet: the capacity notice retries
      setLoad(load);
      if (load.state === "ready") openWorkspace();
      stop = connectEvents((event) => {
        if (event.type === "status") setLoad(event.load);
        else if (event.type === "progress") progress(event.task, event.pct, event.status);
        else if (event.type === "failed") {
          fail(event.task ? `${event.task}: ${event.message}` : event.message);
          if (event.fatal) setLoad({ ...useApp.getState().load, state: "failed", error: event.message });
        } else if (event.type === "ready") openWorkspace();
        else if (event.type === "job") {
          const { setJob, analysesChanged } = useApp.getState();
          setJob(event.layer, { name: event.name, state: event.state, message: event.message });
          if (event.state !== "running") analysesChanged(event.layer);
          if (event.name === "de" && event.state === "done") useApp.getState().showResults(event.layer, "de");
        } else if (event.type === "layer_stale") {
          refreshLayer(event.layer).catch((exc) => console.error(exc));
          useApp.getState().analysesChanged(event.layer);
        }
      });
    });
    // Fallback while loading, in case a proxy blocks the events socket
    const poll = setInterval(async () => {
      if (useApp.getState().load.state !== "running") return;
      let load: LoadStatus;
      try {
        load = await getJson<LoadStatus>("load");
      } catch {
        return;                  // a network blip: try again on the next tick
      }
      if (load.state === "running") return;
      setLoad(load);
      if (load.state === "ready" && useApp.getState().screen !== "workspace") openWorkspace();
    }, 3000);
    return () => { stop(); clearInterval(poll); };
  }, []);

  const notice = <>
    {capacity !== null && <CapacityNotice retryAfter={capacity} onRetry={retrySession} />}
    <CookieNotice />
    <ErrorNotices />
    <TourOverlay />
  </>;
  if (screen === "workspace") return <>{notice}<Workspace /></>;
  if (screen === "loading") return <>{notice}<Loading /></>;
  return <>{notice}<Landing /></>;
}
