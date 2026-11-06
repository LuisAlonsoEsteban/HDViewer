import { useEffect, useRef } from "react";
import { bridge } from "../bokeh/bridge";

/** A Bokeh json_item from the API, registered in the bridge under `docKey`. */
export function BokehFigure({ item, docKey, onReady }: { item: any; docKey: string; onReady?(): void }) {
  const host = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const element = host.current!;
    if (!item) return;
    let current = true;
    bridge.embed(item, element, docKey).then(() => current && onReady?.()).catch((exc) => console.error(exc));
    return () => {
      current = false;
      element.innerHTML = "";
    };
  }, [item, docKey]);
  return <div ref={host} style={{ minWidth: 0 }} />;
}
