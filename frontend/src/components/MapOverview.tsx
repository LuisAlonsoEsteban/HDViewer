import { useState } from "react";
import { ActionIcon, Group, Paper, Slider, Stack, Text } from "@mantine/core";
import { imageOpacity, setImageAlpha } from "../map/layers";
import { MapIcon } from "./icons";
import { Thumbnail } from "./Thumbnail";

/** Map overlay (bottom left): the slide overview with the visible region, the
 * viewport coordinates and the tissue image's opacity. Folds to a button. */
export function MapOverview() {
  const [open, setOpen] = useState(true);
  if (!open) {
    return (
      <ActionIcon className="map-overview-button" size="lg" variant="default" onClick={() => setOpen(true)}
                  aria-label="Show the overview">
        <MapIcon size={18} />
      </ActionIcon>
    );
  }
  return (
    <Paper className="map-overview" p="xs" withBorder shadow="md">
      <Stack gap={6}>
        <Group justify="space-between" gap={4}>
          <Text size="xs" fw={600} tt="uppercase" c="dimmed" style={{ letterSpacing: ".08em" }}>Overview</Text>
          <ActionIcon size="sm" variant="subtle" color="gray" onClick={() => setOpen(false)} aria-label="Hide the overview">−</ActionIcon>
        </Group>
        <Thumbnail />
        <div>
          <Text size="xs" c="dimmed">Image opacity</Text>
          <Slider size="sm" min={0} max={1} step={0.05} defaultValue={imageOpacity()}
                  label={(v) => `${Math.round(v * 100)}%`} onChange={setImageAlpha} thumbLabel="Image opacity" />
        </div>
      </Stack>
    </Paper>
  );
}
