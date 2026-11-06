import { useState } from "react";
import { Accordion, Anchor, AppShell, Button, Group, ScrollArea, Slider, Stack, Tabs, Text, Title } from "@mantine/core";
import { ClusteringPanel } from "../components/ClusteringPanel";
import { ComparisonPanel } from "../components/ComparisonPanel";
import { DEPanel } from "../components/DEPanel";
import { DEResults } from "../components/DEResults";
import { DeconvolutionPanel } from "../components/DeconvolutionPanel";
import { EmbeddingPanel } from "../components/EmbeddingPanel";
import { GroupsPanel } from "../components/GroupsPanel";
import { LayerPanel } from "../components/LayerPanel";
import { MapView } from "../components/MapView";
import { MetadataPanel } from "../components/MetadataPanel";
import { MetadataResults } from "../components/MetadataResults";
import { ParametersPanel } from "../components/ParametersPanel";
import { ROIPanel } from "../components/ROIPanel";
import { SummaryPanel } from "../components/SummaryPanel";
import { Thumbnail } from "../components/Thumbnail";
import { imageOpacity, setImageAlpha } from "../map/layers";
import { useApp } from "../store/app";

export function Workspace() {
  const { workspace, config, results } = useApp();
  const [analyses, setAnalyses] = useState(true);
  if (!workspace) return null;
  return (
    <AppShell navbar={{ width: 420, breakpoint: "sm" }} padding={0}
              aside={{ width: 480, breakpoint: "md", collapsed: { desktop: !analyses, mobile: !analyses } }}>
      <AppShell.Navbar p="sm">
        <AppShell.Section>
          <Group justify="space-between">
            <Title order={4}>{config?.branding.name ?? "VisiumVisualizer"}</Title>
            <Group gap="xs">
              <Button size="compact-xs" variant={analyses ? "light" : "default"} onClick={() => setAnalyses((a) => !a)}
                      aria-pressed={analyses}>Summary · embeddings · comparison</Button>
              {config && <Anchor size="sm" href={config.tutorial_url} target="_blank">Help</Anchor>}
            </Group>
          </Group>
        </AppShell.Section>
        <AppShell.Section grow component={ScrollArea} my="sm">
          <Accordion multiple defaultValue={[workspace.layers[0]?.id]}>
            {workspace.layers.map((layer) => (
              <Accordion.Item key={layer.id} value={layer.id}>
                <Accordion.Control>{layer.title}</Accordion.Control>
                <Accordion.Panel>
                  <Tabs defaultValue="display" keepMounted={false} styles={{ tab: { paddingInline: 6, fontSize: 13 } }}>
                    <Tabs.List grow>
                      <Tabs.Tab value="display">Display</Tabs.Tab>
                      <Tabs.Tab value="clustering">Clustering</Tabs.Tab>
                      <Tabs.Tab value="groups">Groups</Tabs.Tab>
                      <Tabs.Tab value="rois">ROIs</Tabs.Tab>
                      <Tabs.Tab value="de">DE</Tabs.Tab>
                      <Tabs.Tab value="deconvolution">Deconv.</Tabs.Tab>
                      <Tabs.Tab value="metadata">Metadata</Tabs.Tab>
                      <Tabs.Tab value="parameters">Params</Tabs.Tab>
                    </Tabs.List>
                    <Tabs.Panel value="display" pt="sm"><LayerPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="clustering" pt="sm"><ClusteringPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="groups" pt="sm"><GroupsPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="rois" pt="sm"><ROIPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="de" pt="sm"><DEPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="deconvolution" pt="sm"><DeconvolutionPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="metadata" pt="sm"><MetadataPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="parameters" pt="sm"><ParametersPanel layer={layer} /></Tabs.Panel>
                  </Tabs>
                </Accordion.Panel>
              </Accordion.Item>
            ))}
          </Accordion>
        </AppShell.Section>
        <AppShell.Section>
          <Stack gap="xs">
            <Thumbnail />
            <ImageOpacity />
            {workspace.reference && <Text size="xs" c="dark.2">Single-cell reference: {workspace.reference}</Text>}
          </Stack>
        </AppShell.Section>
      </AppShell.Navbar>
      {/* Workspace-wide analyses (the Bokeh UI's Summary, Embeddings and Comparison panels) */}
      <AppShell.Aside p="sm">
        <Tabs defaultValue="summary" keepMounted={false} style={{ display: "flex", flexDirection: "column", minHeight: 0, flex: 1 }}>
          <Tabs.List grow>
            <Tabs.Tab value="summary">Summary</Tabs.Tab>
            <Tabs.Tab value="embeddings">Embeddings</Tabs.Tab>
            <Tabs.Tab value="comparison">Comparison</Tabs.Tab>
          </Tabs.List>
          <ScrollArea style={{ flex: 1 }} pt="sm">
            <Tabs.Panel value="summary"><SummaryPanel /></Tabs.Panel>
            <Tabs.Panel value="embeddings"><EmbeddingPanel /></Tabs.Panel>
            <Tabs.Panel value="comparison"><ComparisonPanel /></Tabs.Panel>
          </ScrollArea>
        </Tabs>
      </AppShell.Aside>
      <AppShell.Main style={{ height: "100vh", display: "flex", flexDirection: "column" }}>
        <div style={{ flex: 1, minHeight: 0, position: "relative" }}><MapView /></div>
        {/* DE and metadata comparison results below the map, as in the Bokeh UI */}
        {results && (
          <div style={{ height: "45%", minHeight: 240 }}>
            {results.kind === "de"
              ? <DEResults key={`de:${results.layer}`} layer={results.layer} run={results.run} />
              : <MetadataResults key={`metadata:${results.layer}`} layer={results.layer} run={results.run} />}
          </div>
        )}
      </AppShell.Main>
    </AppShell>
  );
}

/** Transparency of the tissue image (changed in the browser only). */
function ImageOpacity() {
  return (
    <div>
      <Text size="xs">Image opacity</Text>
      <Slider size="sm" min={0} max={1} step={0.05} defaultValue={imageOpacity()} label={(v) => `${Math.round(v * 100)}%`}
              onChange={setImageAlpha} thumbLabel="Image opacity" />
    </div>
  );
}
