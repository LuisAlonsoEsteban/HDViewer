import { useEffect, useRef, useState } from "react";
import { Button, Stack } from "@mantine/core";
import { getJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { mapPlot } from "../map/layers";
import { useApp } from "../store/app";

/** Visible region as x0,y0,x1,y1 in full-resolution image pixels. */
function viewportBox(inverse: number[] | null): number[] | null {
  const plot = mapPlot();
  if (!plot) return null;
  const [x0, x1, y0, y1] = [plot.x_range.start, plot.x_range.end, plot.y_range.start, plot.y_range.end];
  if (!inverse) return [x0, y0, x1, y1].map(Math.round);
  const corners = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]].map(([x, y]) =>
    [inverse[0] * x + inverse[1] * y + inverse[2], inverse[3] * x + inverse[4] * y + inverse[5]]);
  const xs = corners.map((p) => p[0]), ys = corners.map((p) => p[1]);
  return [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)].map(Math.round);
}

export function Thumbnail() {
  const host = useRef<HTMLDivElement>(null);
  const inverse = useApp((s) => s.workspace?.image.inverse ?? null);
  const [copied, setCopied] = useState(false);

  useEffect(() => {
    let box: any = null;
    const follow = () => {
      const plot = mapPlot();
      if (box && plot) box.setv({ left: plot.x_range.start, right: plot.x_range.end,
                                  bottom: plot.y_range.start, top: plot.y_range.end });
    };
    const off = bridge.on("ranges", follow);
    getJson("figures/thumbnail").then(async (item) => {
      const doc = await bridge.embed(item, host.current!, "thumbnail");
      box = doc.get_model_by_name("viewport_box");
      follow();
    });
    return off;
  }, []);

  const copy = async () => {
    const text = viewportBox(inverse)?.join(",");
    if (!text) return;
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
      setCopied(true);
      setTimeout(() => setCopied(false), 1200);
    } else {
      window.prompt("Visible region (x0,y0,x1,y1, full-res px):", text);
    }
  };

  return (
    <Stack gap={6}>
      <div ref={host} style={{ maxWidth: 150, margin: "0 auto", width: "100%" }} />
      <Button size="compact-xs" variant="default" onClick={copy} aria-label="Copy viewport coordinates"
              title="Copy the visible region as x0,y0,x1,y1 in full-resolution image pixels">
        {copied ? "Copied" : "Copy coordinates"}
      </Button>
    </Stack>
  );
}
