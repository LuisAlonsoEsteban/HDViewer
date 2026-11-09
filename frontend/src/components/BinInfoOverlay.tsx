import { useState } from "react";
import { ActionIcon, Group, Paper, Select, Stack, Text } from "@mantine/core";
import { requestViewport, setStyle } from "../map/layers";
import { Layer, useApp } from "../store/app";
import { HelpButton } from "./HelpButton";

const fmt = (value: number | null, digits = 2) =>
  value === null ? "–" : value.toLocaleString(undefined, { minimumFractionDigits: digits, maximumFractionDigits: digits });

/** Readout of the bins on the map (top right): the aggregation in use, its
 * geometry, and the resolution mode (follow zoom or a fixed level). */
export function BinInfoOverlay({ layer }: { layer: Layer }) {
  const info = useApp((s) => s.binInfo[layer.id]);
  const updateLayer = useApp((s) => s.updateLayer);
  const [open, setOpen] = useState(false);       // folded: the title line shows the level in use
  const [error, setError] = useState<string | null>(null);
  if (!info || !layer.visible) return null;

  const options = [{ value: "0", label: "Dynamic · follow zoom" },
    ...info.factors.filter((f) => info.minimum === null || f >= info.minimum)
      .map((f) => ({ value: String(f), label: f === 1 ? "Fixed · raw/native" : `Fixed · ${f} × ${f}` }))];
  const choose = async (value: string | null) => {
    if (value === null) return;
    setError(null);
    try {
      await setStyle(layer.id, { aggregation: Number(value) });
      updateLayer(layer.id, { fixed_aggregation: Number(value) || null });
      requestViewport(layer.id);
    } catch (exc: any) {
      setError(exc.message);
    }
  };
  const title = `${layer.title.toUpperCase()} · ${info.factor} × ${info.factor}${info.factor === 1 ? " · native" : ""}`;

  return (
    <Paper className="bin-info" p="xs" withBorder data-testid={`bin-info-${layer.id}`}>
      <Group justify="space-between" gap={6} wrap="nowrap">
        <Text ff="monospace" size="xs" fw={600}>{title}</Text>
        <Group gap={4} wrap="nowrap">
          <HelpButton topic="bin_resolution" />
          <ActionIcon size="sm" variant="default" onClick={() => setOpen((o) => !o)}
                      aria-label={open ? "Collapse bin resolution" : "Expand bin resolution"}>{open ? "−" : "+"}</ActionIcon>
        </Group>
      </Group>
      {open && (
        <Stack gap={6} mt={6}>
          <Select size="xs" label="Bin resolution" data={options} allowDeselect={false}
                  value={info.fixed === null ? "0" : String(info.factor)} onChange={choose} />
          <Text ff="monospace" size="xs" style={{ whiteSpace: "pre" }}>
            {[`Width      ${fmt(info.width)} plot units`,
              `Height     ${fmt(info.height)} plot units`,
              `Rotation   ${fmt(info.rotation)}°`,
              ...(info.bin_size_um !== null ? [`Native bin ${info.bin_size_um} µm`] : []),
              info.fixed === null ? "Resolution Zoom-driven" : "Resolution Fixed",
              ...(info.minimum !== null ? [`Cluster minimum ${info.minimum} × ${info.minimum}`] : [])].join("\n")}
          </Text>
          {info.minimum !== null && (
            <Text size="xs" c="dimmed">
              Cluster columns require at least {info.minimum} × {info.minimum}. Dynamic zoom and coarser fixed
              levels remain available.
            </Text>
          )}
          {error && <Text size="xs" c="red.4">{error}</Text>}
        </Stack>
      )}
    </Paper>
  );
}
