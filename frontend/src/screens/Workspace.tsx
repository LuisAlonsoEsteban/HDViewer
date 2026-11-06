import { useState } from "react";
import { Accordion, Anchor, AppShell, Button, Group, ScrollArea, Slider, Stack, Tabs, Text, Title, UnstyledButton } from "@mantine/core";
import { ClusteringPanel } from "../components/ClusteringPanel";
import { ComparisonPanel } from "../components/ComparisonPanel";
import { DEPanel } from "../components/DEPanel";
import { DEResults } from "../components/DEResults";
import { DeconvolutionPanel } from "../components/DeconvolutionPanel";
import { EmbeddingPanel } from "../components/EmbeddingPanel";
import { GroupsPanel } from "../components/GroupsPanel";
import { HelpButton } from "../components/HelpButton";
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
import { Splitter, useDragSize } from "../components/Splitter";
import { BrandMark } from "./Landing";

export function Workspace() {
  const { workspace, config, results, setScreen } = useApp();
  const [analyses, setAnalyses] = useState(true);
  const [left, setLeft] = useDragSize("left", 420, 300, 760);
  const [right, setRight] = useDragSize("right", 480, 320, 900);
  const [below, setBelow] = useDragSize("results", Math.round(window.innerHeight * 0.45), 160, Math.round(window.innerHeight * 0.8));
  if (!workspace) return null;
  return (
    <AppShell navbar={{ width: left, breakpoint: "sm" }} padding={0}
              aside={{ width: right, breakpoint: "md", collapsed: { desktop: !analyses, mobile: !analyses } }}>
      <AppShell.Navbar p="sm">
        <Splitter axis="x" size={left} onResize={setLeft} label="Resize the layer panel"
                  style={{ position: "absolute", top: 0, bottom: 0, right: -4 }} />
        <AppShell.Section>
          <Group justify="space-between" wrap="nowrap">
            <UnstyledButton onClick={() => setScreen("landing")} aria-label="About (back to the overview)">
              <Group gap={8} wrap="nowrap">
                <BrandMark className="workspace-mark" />
                <Title order={4}>{config?.branding.name ?? "VisiumVisualizer"}</Title>
              </Group>
            </UnstyledButton>
            <Group gap="sm" wrap="nowrap">
              {config && <Anchor size="xs" href={config.tutorial_url} target="_blank" rel="noopener">Guide ↗</Anchor>}
              {config?.branding.repository && (
                <Anchor size="xs" href={config.branding.repository} target="_blank" rel="noopener">GitHub ↗</Anchor>
              )}
              <Anchor size="xs" component="button" onClick={() => setScreen("landing")}>About</Anchor>
            </Group>
          </Group>
          <Button mt="xs" size="compact-xs" variant={analyses ? "light" : "default"} onClick={() => setAnalyses((a) => !a)}
                  aria-pressed={analyses}>Summary · embeddings · comparison</Button>
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
                    <Tabs.Panel value="display" pt="sm"><PanelHelp topics={["general", "color"]} /><LayerPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="clustering" pt="sm"><PanelHelp topics={["clustering"]} /><ClusteringPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="groups" pt="sm"><PanelHelp topics={["groups"]} /><GroupsPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="rois" pt="sm"><PanelHelp topics={["rois"]} /><ROIPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="de" pt="sm"><PanelHelp topics={["de"]} /><DEPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="deconvolution" pt="sm"><PanelHelp topics={["deconvolution"]} /><DeconvolutionPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="metadata" pt="sm"><PanelHelp topics={["metadata_comparison"]} /><MetadataPanel layer={layer} /></Tabs.Panel>
                    <Tabs.Panel value="parameters" pt="sm"><PanelHelp topics={["parameters"]} /><ParametersPanel layer={layer} /></Tabs.Panel>
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
        <Splitter axis="x" sign={-1} size={right} onResize={setRight} label="Resize the analyses panel"
                  style={{ position: "absolute", top: 0, bottom: 0, left: -4 }} />
        <Tabs defaultValue="summary" keepMounted={false} style={{ display: "flex", flexDirection: "column", minHeight: 0, flex: 1 }}>
          <Tabs.List grow>
            <Tabs.Tab value="summary">Summary</Tabs.Tab>
            <Tabs.Tab value="embeddings">Embeddings</Tabs.Tab>
            <Tabs.Tab value="comparison">Comparison</Tabs.Tab>
          </Tabs.List>
          <ScrollArea style={{ flex: 1 }} pt="sm">
            <Tabs.Panel value="summary"><PanelHelp topics={["metadata"]} /><SummaryPanel /></Tabs.Panel>
            <Tabs.Panel value="embeddings"><PanelHelp topics={["embedding"]} /><EmbeddingPanel /></Tabs.Panel>
            <Tabs.Panel value="comparison"><ComparisonPanel /></Tabs.Panel>
          </ScrollArea>
        </Tabs>
      </AppShell.Aside>
      <AppShell.Main style={{ height: "100vh", display: "flex", flexDirection: "column" }}>
        <div style={{ flex: 1, minHeight: 0, position: "relative" }}><MapView /></div>
        {/* DE and metadata comparison results below the map, as in the Bokeh UI */}
        {results && (
          <div style={{ height: below, position: "relative", flex: "none" }}>
            <Splitter axis="y" sign={-1} size={below} onResize={setBelow} label="Resize the results panel"
                      style={{ position: "absolute", left: 0, right: 0, top: -4 }} />
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

/** The help buttons of a tab (the old interface's "?" next to each section). */
function PanelHelp({ topics }: { topics: string[] }) {
  return (
    <Group gap={4} justify="flex-end" mb={4}>
      <Text size="xs" c="dimmed">Help</Text>
      {topics.map((topic) => <HelpButton key={topic} topic={topic} />)}
    </Group>
  );
}
