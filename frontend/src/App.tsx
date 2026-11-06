import { useEffect } from "react";
import { connectEvents, getJson, LoadStatus } from "./api/client";
import { refreshLayer } from "./map/layers";
import { Config, useApp, Workspace as WorkspaceData } from "./store/app";
import { Landing } from "./screens/Landing";
import { Loading } from "./screens/Loading";
import { Workspace } from "./screens/Workspace";

export function App() {
  const { screen, setConfig, setLoad, progress, fail, setWorkspace, setScreen } = useApp();

  useEffect(() => {
    const openWorkspace = async () => {
      setWorkspace(await getJson<WorkspaceData>("workspace"));
      setScreen("workspace");
    };
    getJson<Config>("config").then((config) => {
      setConfig(config);
      if (window.Bokeh?.version && window.Bokeh.version !== config.bokeh_version) {
        console.warn(`BokehJS ${window.Bokeh.version} differs from the server's Bokeh ${config.bokeh_version}`);
      }
    });
    let stop = () => {};
    // GET /api/load creates the session cookie the events socket needs
    getJson<LoadStatus>("load").then((load) => {
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
          refreshLayer(event.layer);
          useApp.getState().analysesChanged(event.layer);
        }
      });
    });
    // Fallback while loading, in case a proxy blocks the events socket
    const poll = setInterval(async () => {
      if (useApp.getState().load.state !== "running") return;
      const load = await getJson<LoadStatus>("load");
      if (load.state === "running") return;
      setLoad(load);
      if (load.state === "ready" && useApp.getState().screen !== "workspace") openWorkspace();
    }, 3000);
    return () => { stop(); clearInterval(poll); };
  }, []);

  if (screen === "workspace") return <Workspace />;
  if (screen === "loading") return <Loading />;
  return <Landing />;
}
