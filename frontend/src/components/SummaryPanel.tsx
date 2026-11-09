import { useEffect, useState } from "react";
import { Loader, ScrollArea, SegmentedControl, Stack, Table, Tabs, Text } from "@mantine/core";
import { getJson } from "../api/client";
import { useApp } from "../store/app";
import { BokehFigure } from "./BokehFigure";

interface Stats {
  id: string;
  title: string;
  kind: string;
  total: number;
  sequenced: number;
  not_sequenced: number;
  pct_sequenced: number;
  mean_umis: number;
  umis_per_um2: number;
}
interface Marker { cluster: string; gene: string; score: number; pval: number | null; logfc: number }

/** Summary panel: dataset table, then QC and marker genes per layer. */
export function SummaryPanel() {
  const layers = useApp((s) => s.workspace?.layers ?? []);
  const versions = useApp((s) => s.analysisVersion);
  const [stats, setStats] = useState<Stats[]>();
  const [layer, setLayer] = useState(layers[0]?.id ?? "");
  const [view, setView] = useState("qc");
  const [qc, setQc] = useState<any>(null);
  const [markers, setMarkers] = useState<{ clustered: boolean; markers: Marker[] }>();
  const version = versions[layer] ?? 0;
  const base = `layers/${encodeURIComponent(layer)}`;

  useEffect(() => { getJson<{ layers: Stats[] }>("summary").then((s) => setStats(s.layers)); }, []);
  useEffect(() => {
    setQc(null);
    if (layer) getJson(`${base}/qc`).then(setQc).catch((exc) => console.error(exc));
  }, [base]);
  useEffect(() => { if (layer) getJson<typeof markers>(`${base}/markers`).then(setMarkers); }, [base, version]);

  const unit = layers.find((l) => l.id === layer)?.kind.replace(/s$/, "") ?? "observation";
  return (
    <Stack gap="sm">
      <Text fw={600} size="sm">Dataset summary</Text>
      <Table.ScrollContainer minWidth={420}>
        <Table fz="xs" striped>
          <Table.Thead>
            <Table.Tr>
              {["Sample", "Total", "Sequenced", "Not seq", "% Seq", "Mean UMIs", "UMIs/µm²"].map((h) => <Table.Th key={h}>{h}</Table.Th>)}
            </Table.Tr>
          </Table.Thead>
          <Table.Tbody>
            {(stats ?? []).map((s) => (
              <Table.Tr key={s.id}>
                <Table.Td fw={600}>{s.title}</Table.Td>
                <Table.Td>{s.total.toLocaleString()}</Table.Td>
                <Table.Td>{s.sequenced.toLocaleString()}</Table.Td>
                <Table.Td>{s.not_sequenced.toLocaleString()}</Table.Td>
                <Table.Td>{s.pct_sequenced}%</Table.Td>
                <Table.Td>{s.mean_umis}</Table.Td>
                <Table.Td>{s.umis_per_um2}</Table.Td>
              </Table.Tr>
            ))}
          </Table.Tbody>
        </Table>
      </Table.ScrollContainer>

      {layers.length > 1 && (
        <SegmentedControl size="xs" value={layer} onChange={setLayer} data={layers.map((l) => ({ value: l.id, label: l.title }))} />
      )}
      <Tabs value={view} onChange={(v) => setView(v ?? "qc")}>
        <Tabs.List>
          <Tabs.Tab value="qc">QC</Tabs.Tab>
          <Tabs.Tab value="markers">Markers</Tabs.Tab>
        </Tabs.List>
        <Tabs.Panel value="qc" pt="xs">
          {!qc && <Loader size="sm" />}
          <BokehFigure item={qc} docKey="qc" />
          {qc && qc.shown < qc.total && (
            <Text size="xs" c="dark.2">Scatters show a random {qc.shown.toLocaleString()} of {qc.total.toLocaleString()} {unit}s.</Text>
          )}
        </Tabs.Panel>
        <Tabs.Panel value="markers" pt="xs">
          {markers && !markers.markers.length && (
            <Text size="xs" c="dark.2">
              {markers.clustered ? "Marker genes appear here once they are ranked after clustering." : "Run clustering to populate."}
            </Text>
          )}
          {!!markers?.markers.length && (
            <ScrollArea.Autosize mah={350}>
              <Text fw={600} size="sm" mb={4}>Top marker genes</Text>
              <Table fz="xs" stickyHeader>
                <Table.Thead>
                  <Table.Tr>{["Cluster", "Gene", "Z-Score", "Adj P-val", "Log2FC"].map((h) => <Table.Th key={h}>{h}</Table.Th>)}</Table.Tr>
                </Table.Thead>
                <Table.Tbody>
                  {markers.markers.map((m) => (
                    <Table.Tr key={`${m.cluster}:${m.gene}`}>
                      <Table.Td c="dark.1">{m.cluster}</Table.Td>
                      <Table.Td>{m.gene}</Table.Td>
                      <Table.Td c="blue.4">{m.score}</Table.Td>
                      <Table.Td c="dark.1">{m.pval && m.pval > 0 ? m.pval.toExponential(2) : "0.00e+00"}</Table.Td>
                      <Table.Td style={{ color: m.logfc >= 0 ? "#4ef1a1" : "#f14e4e" }}>{m.logfc}</Table.Td>
                    </Table.Tr>
                  ))}
                </Table.Tbody>
              </Table>
            </ScrollArea.Autosize>
          )}
        </Tabs.Panel>
      </Tabs>
    </Stack>
  );
}
