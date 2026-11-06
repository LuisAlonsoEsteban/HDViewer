import { useEffect, useRef, useState } from "react";
import { CloseButton, Group, Loader, ScrollArea, Text, Title } from "@mantine/core";
import { getJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { useApp } from "../store/app";

/** A layer's last metadata comparison: one bar chart per comparison. */
export function MetadataResults({ layer, run }: { layer: string; run: number }) {
  const version = useApp((s) => s.analysisVersion[layer] ?? 0);
  const showResults = useApp((s) => s.showResults);
  const title = useApp((s) => s.workspace?.layers.find((l) => l.id === layer)?.title ?? layer);
  const [status, setStatus] = useState<"loading" | "ready" | "cleared" | "failed">("loading");
  const plots = useRef<HTMLDivElement>(null);
  const base = `layers/${encodeURIComponent(layer)}/metadata`;

  // Charts of each new comparison; cleared when the compared sets change
  useEffect(() => {
    let current = true;
    const element = plots.current!;
    (async () => {
      const { ready } = await getJson<{ ready: boolean }>(`${base}/options`);
      if (!current) return;
      if (!ready) {
        element.innerHTML = "";
        return setStatus("cleared");
      }
      const item = await getJson(`${base}/plots`);
      if (!current) return;
      await bridge.embed(item, element, "metadata");
      setStatus("ready");
    })().catch((exc) => {
      console.error(exc);
      if (current) setStatus("failed");
    });
    return () => { current = false; };
  }, [base, run, version]);

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%", borderTop: "1px solid #344d45" }}
         data-testid="metadata-results">
      <Group justify="space-between" px="md" py={6} wrap="nowrap">
        <Title order={5}>Metadata comparison <Text span size="xs" c="dark.2">· {title}</Text></Title>
        <CloseButton aria-label="Close metadata comparison" onClick={() => showResults(null)} />
      </Group>
      <ScrollArea style={{ flex: 1 }} px="md" pb="md">
        {status === "loading" && <Loader size="sm" />}
        {status === "cleared" && (
          <Text size="sm" c="dark.2">
            The metadata comparison was cleared because the compared groups or clusters changed. Compare again.
          </Text>
        )}
        {status === "failed" && <Text size="sm" c="red.4">The charts could not be shown.</Text>}
        <div ref={plots} style={{ minWidth: 0 }} />
      </ScrollArea>
    </div>
  );
}
