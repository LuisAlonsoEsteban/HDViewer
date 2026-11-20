// Hosts <Vitessce> with the config at ?config=<url>. The benchmark driver
// changes the view (target, zoom) and the coloured gene by passing a new
// config object with the same uid: Vitessce updates the view in place, as
// its Jupyter widget does when Python edits the config.
import React, { useEffect, useState } from "react";
import { createRoot } from "react-dom/client";
import { Vitessce } from "vitessce";

function setAll(space, type, value) {
  if (!space[type]) return;
  for (const scope of Object.keys(space[type])) space[type][scope] = value;
}

function App() {
  const [config, setConfig] = useState(null);
  useEffect(() => {
    const url = new URLSearchParams(window.location.search).get("config");
    fetch(url).then((r) => r.json()).then((c) => { window.__vitConfigLoaded = performance.now(); setConfig(c); });
    window.__setView = (x, y, zoom) => setConfig((c) => {
      const next = structuredClone(c);
      setAll(next.coordinationSpace, "spatialTargetX", x);
      setAll(next.coordinationSpace, "spatialTargetY", y);
      setAll(next.coordinationSpace, "spatialZoom", zoom);
      return next;
    });
    window.__setGene = (gene) => setConfig((c) => {
      const next = structuredClone(c);
      setAll(next.coordinationSpace, "featureSelection", [gene]);
      setAll(next.coordinationSpace, "obsColorEncoding", "geneSelection");
      return next;
    });
  }, []);
  if (!config) return null;
  return <Vitessce config={config} height={window.innerHeight} theme="light" remountOnUidChange={false} />;
}

createRoot(document.getElementById("root")).render(<App />);
