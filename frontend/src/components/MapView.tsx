import { useEffect, useRef, useState } from "react";
import { Loader, Text } from "@mantine/core";
import { getJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { imageOpacity, requestViewport, sendSelection, setImageAlpha, setMapHost } from "../map/layers";
import { useApp } from "../store/app";

/** The tissue map: a Bokeh figure from the server whose layers load lazily. */
export function MapView() {
  const host = useRef<HTMLDivElement>(null);
  const [status, setStatus] = useState<"loading" | "ready" | "error">("loading");

  useEffect(() => {
    let timer: ReturnType<typeof setTimeout> | undefined;
    const visibleLayers = () => (useApp.getState().workspace?.layers ?? []).filter((l) => l.visible);
    const reload = () => {
      clearTimeout(timer);
      timer = setTimeout(() => visibleLayers().forEach((layer) => requestViewport(layer.id)), 120);
    };
    const offRanges = bridge.on("ranges", reload);
    const offSelection = bridge.on("selection", ({ layer, indices }) => sendSelection(layer, indices));
    setMapHost(host.current);
    (async () => {
      try {
        const item = await getJson("figures/map");
        await bridge.embed(item, host.current!, "map");
        setImageAlpha(imageOpacity());
        setStatus("ready");
        reload();
      } catch (exc) {
        console.error(exc);
        setStatus("error");
      }
    })();
    return () => {
      offRanges();
      offSelection();
      clearTimeout(timer);
      setMapHost(null);
    };
  }, []);

  return (
    <div style={{ position: "relative", height: "100%" }}>
      <div ref={host} className="map-host" />
      {status === "loading" && <Loader style={{ position: "absolute", top: "50%", left: "50%" }} />}
      {status === "error" && <Text c="red.4" p="md">The map could not be loaded.</Text>}
    </div>
  );
}
