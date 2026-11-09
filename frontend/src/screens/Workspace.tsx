import { useState } from "react";
import type { ReactNode } from "react";
import { Anchor, AppShell, Button, Divider, Group, ScrollArea, Tabs, Text, Title, Tooltip, UnstyledButton } from "@mantine/core";
import { ClusteringPanel } from "../components/ClusteringPanel";
import { ComparisonPanel } from "../components/ComparisonPanel";
import { DEPanel } from "../components/DEPanel";
import { DEResults } from "../components/DEResults";
import { DeconvolutionPanel } from "../components/DeconvolutionPanel";
import { EmbeddingPanel } from "../components/EmbeddingPanel";
import { GroupsPanel } from "../components/GroupsPanel";
import { HelpButton } from "../components/HelpButton";
import {
  ClusteringIcon, DeconvIcon, DeIcon, DisplayIcon, GroupsIcon, MetadataIcon, PanelRightIcon, ParamsIcon, RoiIcon,
} from "../components/icons";
import { LayerPanel } from "../components/LayerPanel";
import { LayerSwitcher } from "../components/LayerSwitcher";
import { MapOverview } from "../components/MapOverview";
import { MapView } from "../components/MapView";
import { MetadataPanel } from "../components/MetadataPanel";
import { MetadataResults } from "../components/MetadataResults";
import { ParametersPanel } from "../components/ParametersPanel";
import { ROIPanel } from "../components/ROIPanel";
import { Splitter, useDragSize } from "../components/Splitter";
import { SummaryPanel } from "../components/SummaryPanel";
import { Layer, useApp } from "../store/app";
import { BrandMark } from "./Landing";

const HEADER = 52;
const RAIL = 76;

interface Tool {
  value: string;
  label: string;             // rail label (also the tab's accessible name)
  title: string;             // panel heading
  topics: string[];          // help topics of the panel
  icon: (props: { size?: number }) => ReactNode;
  render(layer: Layer): ReactNode;
}

// Tools that act on the active layer, in the order of a typical analysis
const TOOLS: Tool[] = [
  { value: "display", label: "Display", title: "Display", topics: ["general", "color"], icon: DisplayIcon,
    render: (layer) => <LayerPanel layer={layer} /> },
  { value: "clustering", label: "Clustering", title: "Clustering", topics: ["clustering"], icon: ClusteringIcon,
    render: (layer) => <ClusteringPanel layer={layer} /> },
  { value: "groups", label: "Groups", title: "Custom groups", topics: ["groups"], icon: GroupsIcon,
    render: (layer) => <GroupsPanel layer={layer} /> },
  { value: "rois", label: "ROIs", title: "Regions of interest", topics: ["rois"], icon: RoiIcon,
    render: (layer) => <ROIPanel layer={layer} /> },
  { value: "de", label: "DE", title: "Differential expression", topics: ["de"], icon: DeIcon,
    render: (layer) => <DEPanel layer={layer} /> },
  { value: "deconvolution", label: "Deconv.", title: "Deconvolution", topics: ["deconvolution"], icon: DeconvIcon,
    render: (layer) => <DeconvolutionPanel layer={layer} /> },
  { value: "metadata", label: "Metadata", title: "Metadata comparison", topics: ["metadata_comparison"],
    icon: MetadataIcon, render: (layer) => <MetadataPanel layer={layer} /> },
  { value: "parameters", label: "Params", title: "Drawing parameters", topics: ["parameters"], icon: ParamsIcon,
    render: (layer) => <ParametersPanel layer={layer} /> },
];

export function Workspace() {
  const { workspace, config, results, setScreen } = useApp();
  const activeId = useApp((s) => s.activeLayer);
  const [tool, setTool] = useState<string>("display");
  const [analyses, setAnalyses] = useState(true);
  const [left, setLeft] = useDragSize("tools", 360, 300, 720);
  const [right, setRight] = useDragSize("right", 440, 320, 900);
  const [below, setBelow] = useDragSize("results", Math.round(window.innerHeight * 0.42), 160, Math.round(window.innerHeight * 0.8));
  if (!workspace) return null;
  const layer = workspace.layers.find((l) => l.id === activeId) ?? workspace.layers[0];
  const current = TOOLS.find((t) => t.value === tool) ?? TOOLS[0];
  const branding = config?.branding ?? {};

  return (
    <AppShell header={{ height: HEADER }} navbar={{ width: RAIL + left, breakpoint: "sm" }} padding={0}
              aside={{ width: right, breakpoint: "md", collapsed: { desktop: !analyses, mobile: !analyses } }}>
      <AppShell.Header className="ws-header">
        <Group h="100%" px="md" justify="space-between" wrap="nowrap">
          <Group gap="md" wrap="nowrap">
            <UnstyledButton onClick={() => setScreen("landing")} aria-label="About (back to the overview)">
              <Group gap={8} wrap="nowrap">
                <BrandMark className="workspace-mark" />
                <Title order={4} className="ws-brand">{branding.name ?? "VisiumVisualizer"}</Title>
              </Group>
            </UnstyledButton>
            <Divider orientation="vertical" my={14} visibleFrom="sm" />
            <Text size="xs" c="dimmed" tt="uppercase" fw={600} visibleFrom="sm" style={{ letterSpacing: ".12em" }}>
              Analysis workspace
            </Text>
            {workspace.reference && (
              <Text size="xs" c="dimmed" visibleFrom="lg" truncate maw={320}>
                · Single-cell reference: {workspace.reference}
              </Text>
            )}
          </Group>
          <Group gap="md" wrap="nowrap">
            {config && <Anchor size="sm" href={config.tutorial_url} target="_blank" rel="noopener" c="dark.1">Guide ↗</Anchor>}
            {branding.repository && (
              <Anchor size="sm" href={branding.repository} target="_blank" rel="noopener" c="dark.1" visibleFrom="sm">GitHub ↗</Anchor>
            )}
            <Anchor size="sm" component="button" onClick={() => setScreen("landing")} c="dark.1">About</Anchor>
            <Button size="xs" variant={analyses ? "light" : "default"} leftSection={<PanelRightIcon size={16} />}
                    onClick={() => setAnalyses((a) => !a)} aria-pressed={analyses}>
              Summary · embeddings · comparison
            </Button>
          </Group>
        </Group>
      </AppShell.Header>

      {/* Tools: a rail of the analysis steps, and the selected tool for the active layer */}
      <AppShell.Navbar className="ws-nav">
        <Tabs value={current.value} onChange={(value) => value && setTool(value)} orientation="vertical"
              keepMounted={false} className="tool-tabs" variant="unstyled">
          <Tabs.List className="tool-rail" aria-label="Tools">
            {TOOLS.map(({ value, label, title, icon: Icon }) => (
              <Tooltip key={value} label={title} position="right" withArrow openDelay={500} disabled={title === label}>
                <Tabs.Tab value={value} className="tool-tab">
                  <span className="tool-tab-icon"><Icon size={20} /></span>
                  <span className="tool-tab-label">{label}</span>
                </Tabs.Tab>
              </Tooltip>
            ))}
          </Tabs.List>
          <div className="tool-panel">
            <div className="tool-panel-top">
              <Text size="xs" c="dimmed" tt="uppercase" fw={600} mb={6} style={{ letterSpacing: ".1em" }}>Layers</Text>
              <LayerSwitcher />
            </div>
            <Group justify="space-between" className="tool-panel-heading" wrap="nowrap">
              <div style={{ minWidth: 0 }}>
                <Title order={5} lineClamp={1}>{current.title}</Title>
                <Text size="xs" c="dimmed" truncate>{layer.title}</Text>
              </div>
              <HelpButton topics={current.topics} />
            </Group>
            <ScrollArea className="tool-panel-body" type="auto" offsetScrollbars>
              {TOOLS.map((t) => (
                <Tabs.Panel key={t.value} value={t.value} p="md" pt="sm">
                  {/* One instance per layer: switching layers starts from that layer's state */}
                  <div key={layer.id}>{t.render(layer)}</div>
                </Tabs.Panel>
              ))}
            </ScrollArea>
          </div>
        </Tabs>
        <Splitter axis="x" size={left} onResize={setLeft} label="Resize the layer panel"
                  style={{ position: "absolute", top: 0, bottom: 0, right: -4 }} />
      </AppShell.Navbar>

      {/* Workspace-wide analyses */}
      <AppShell.Aside className="ws-aside">
        <Splitter axis="x" sign={-1} size={right} onResize={setRight} label="Resize the analyses panel"
                  style={{ position: "absolute", top: 0, bottom: 0, left: -4 }} />
        <Tabs defaultValue="summary" keepMounted={false} className="aside-tabs">
          <Tabs.List grow px="sm">
            <Tabs.Tab value="summary">Summary</Tabs.Tab>
            <Tabs.Tab value="embeddings">Embeddings</Tabs.Tab>
            <Tabs.Tab value="comparison">Comparison</Tabs.Tab>
          </Tabs.List>
          <ScrollArea className="aside-body" type="auto" offsetScrollbars>
            <Tabs.Panel value="summary" p="md"><PanelHelp topics={["metadata"]} /><SummaryPanel /></Tabs.Panel>
            <Tabs.Panel value="embeddings" p="md"><PanelHelp topics={["embedding"]} /><EmbeddingPanel /></Tabs.Panel>
            <Tabs.Panel value="comparison" p="md"><ComparisonPanel /></Tabs.Panel>
          </ScrollArea>
        </Tabs>
      </AppShell.Aside>

      <AppShell.Main className="ws-main" style={{ height: `calc(100dvh - ${HEADER}px)` }}>
        <div className="map-area">
          <MapView />
          <MapOverview />
        </div>
        {/* DE and metadata comparison results open below the map */}
        {results && (
          <div className="results-area" style={{ height: below }}>
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

/** The help buttons of a panel, top right. */
function PanelHelp({ topics }: { topics: string[] }) {
  return (
    <Group gap={4} justify="flex-end" mb={-24} style={{ position: "relative", zIndex: 1 }}>
      {topics.map((topic) => <HelpButton key={topic} topic={topic} />)}
    </Group>
  );
}
