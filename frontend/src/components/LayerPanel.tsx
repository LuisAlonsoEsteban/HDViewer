import { useEffect, useState } from "react";
import { Autocomplete, Button, ColorInput, Group, Select, Slider, Stack, Text } from "@mantine/core";
import { getJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { requestViewport, sendSelection, setStyle } from "../map/layers";
import { Layer, useApp } from "../store/app";

/** Display tool: colouring, opacity and data of one layer. */
export function LayerPanel({ layer }: { layer: Layer }) {
  const { config, updateLayer } = useApp();
  const [query, setQuery] = useState(layer.colour_by.name ?? "");
  const [options, setOptions] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const timer = setTimeout(() => {
      getJson<string[]>(`layers/${encodeURIComponent(layer.id)}/columns?q=${encodeURIComponent(query)}`)
        .then(setOptions).catch(() => setOptions([]));
    }, 150);
    return () => clearTimeout(timer);
  }, [query, layer.id]);

  const change = async (changes: Record<string, unknown>, local: Partial<Layer>, reload = false) => {
    setError(null);
    try {
      await setStyle(layer.id, changes);
      updateLayer(layer.id, local);
      if (reload) requestViewport(layer.id);
    } catch (exc: any) {
      setError(exc.message);
    }
  };

  return (
    <Stack gap="sm">
      {!layer.visible && <Text size="xs" c="yellow.4">This layer is hidden: show it with the eye in the layer list.</Text>}
      <Autocomplete label="Gene / metadata" placeholder="Start typing" data={options} value={query} onChange={setQuery}
                    disabled={!layer.visible}
                    onOptionSubmit={(name) => change({ colour_by: name }, { colour_by: { kind: "column", name } })} />
      <Select label="Expression values" disabled={!layer.visible} value={layer.scale}
              data={(config?.expression_scales ?? []).map(({ value, label }) => ({ value, label }))}
              onChange={(scale) => scale && change({ scale }, { scale })} />
      <ColorInput label="Colour" value={layer.colour} disabled={!layer.visible || layer.colour_by.kind === "category"}
                  onChangeEnd={(colour) => change({ colour }, { colour })} />
      <div>
        <Text size="sm">Opacity</Text>
        <Slider min={0} max={1} step={0.05} defaultValue={layer.alpha} disabled={!layer.visible}
                onChangeEnd={(alpha) => change({ alpha }, { alpha })} />
      </div>
      {layer.kind === "bins" && (
        <Text size="xs" c="dimmed">Bin resolution is set in the box at the top right of the map.</Text>
      )}
      {error && <Text size="xs" c="red.4">{error}</Text>}
      <Group grow gap="xs" mt={4}>
        <Button size="xs" variant="default" disabled={!layer.visible}
                onClick={async () => {
                  bridge.selectQuietly(bridge.model(layer.id, "source"), []);
                  await sendSelection(layer.id, []);
                }}>Clear selection</Button>
        {/* The layer's AnnData with its clusters, groups and deconvolution results */}
        <Button component="a" size="xs" variant="default" href={`api/layers/${encodeURIComponent(layer.id)}/export`}
                download={`${layer.id}.h5ad`} title="The layer's data as AnnData, with its clusters, groups and deconvolution results">
          Download .h5ad
        </Button>
      </Group>
    </Stack>
  );
}
