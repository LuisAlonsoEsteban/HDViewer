import { useState } from "react";
import { ActionIcon, Group, Stack, Text, Tooltip, UnstyledButton } from "@mantine/core";
import { requestViewport, setStyle } from "../map/layers";
import { Layer, useApp } from "../store/app";
import { EyeIcon, EyeOffIcon } from "./icons";

/** The workspace's layers: pick the one the tools work on, and show or hide each. */
export function LayerSwitcher() {
  const layers = useApp((s) => s.workspace?.layers ?? []);
  const active = useApp((s) => s.activeLayer);
  const setActive = useApp((s) => s.setActiveLayer);
  return (
    <Stack gap={4} role="radiogroup" aria-label="Active layer">
      {layers.map((layer) => (
        <LayerRow key={layer.id} layer={layer} active={layer.id === active} onSelect={() => setActive(layer.id)} />
      ))}
    </Stack>
  );
}

function LayerRow({ layer, active, onSelect }: { layer: Layer; active: boolean; onSelect(): void }) {
  const updateLayer = useApp((s) => s.updateLayer);
  const [busy, setBusy] = useState(false);
  const toggle = async () => {
    const visible = !layer.visible;
    setBusy(true);
    try {
      await setStyle(layer.id, { visible });
      updateLayer(layer.id, { visible });
      if (visible) requestViewport(layer.id);
    } finally {
      setBusy(false);
    }
  };
  const colourBy = layer.colour_by.name
    ? `${layer.colour_by.kind === "clusters" ? "Clusters" : layer.colour_by.kind === "groups" ? "Groups" : layer.colour_by.name}`
    : "Fixed colour";
  return (
    <Group gap={0} wrap="nowrap" className="layer-row" data-active={active || undefined} data-hidden={!layer.visible || undefined}>
      <UnstyledButton role="radio" aria-checked={active} onClick={onSelect} className="layer-row-select"
                      aria-label={`${layer.title} layer`}>
        <Group gap={10} wrap="nowrap">
          <span className="layer-swatch" style={{ background: layer.colour_by.kind === "gene" || layer.colour_by.kind === "metadata" || !layer.colour_by.name ? layer.colour : undefined }}
                data-categorical={layer.colour_by.kind === "category" || layer.colour_by.kind === "clusters" || layer.colour_by.kind === "groups" || undefined} />
          <div style={{ minWidth: 0 }}>
            <Text size="sm" fw={600} truncate>{layer.title}</Text>
            <Text size="xs" c="dimmed" truncate>{layer.visible ? colourBy : "Hidden"}</Text>
          </div>
        </Group>
      </UnstyledButton>
      <Tooltip label={layer.visible ? "Hide layer" : "Show layer"} withArrow openDelay={400}>
        <ActionIcon variant="subtle" color="gray" onClick={toggle} loading={busy}
                    aria-label={`${layer.visible ? "Hide" : "Show"} ${layer.title.toLowerCase()}`} aria-pressed={layer.visible}>
          {layer.visible ? <EyeIcon size={18} /> : <EyeOffIcon size={18} />}
        </ActionIcon>
      </Tooltip>
    </Group>
  );
}
