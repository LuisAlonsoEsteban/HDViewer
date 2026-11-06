import { useEffect, useState } from "react";
import { Autocomplete, Button, Checkbox, ColorInput, Select, Slider, Stack, Switch, Text } from "@mantine/core";
import { getJson } from "../api/client";
import { requestViewport, setStyle } from "../map/layers";
import { Layer, useApp } from "../store/app";

/** Display controls of one layer (the Bokeh UI's General and Color tabs). */
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
      <Switch label={`Show ${layer.title.toLowerCase()}`} checked={layer.visible}
              onChange={(e) => change({ visible: e.currentTarget.checked }, { visible: e.currentTarget.checked }, true)} />
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
      <Checkbox label="Show not sequenced" checked={layer.show_not_sequenced} disabled={!layer.visible}
                onChange={(e) => change({ show_not_sequenced: e.currentTarget.checked },
                                       { show_not_sequenced: e.currentTarget.checked })} />
      {layer.kind === "bins" && (
        <Select label="Bin resolution" disabled={!layer.visible}
                value={String(layer.fixed_aggregation ?? 0)}
                data={[{ value: "0", label: "Follow zoom" },
                       ...layer.aggregation_factors.map((f) => ({ value: String(f), label: f === 1 ? "Native" : `${f} × ${f} bins` }))]}
                onChange={(value) => change({ aggregation: Number(value) },
                                            { fixed_aggregation: Number(value) || null }, true)} />
      )}
      {error && <Text size="xs" c="red.4">{error}</Text>}
      {/* The layer's AnnData with its clusters, groups and deconvolution results */}
      <Button component="a" size="xs" variant="default" href={`api/layers/${encodeURIComponent(layer.id)}/export`}
              download={`${layer.id}.h5ad`}>Download data (.h5ad)</Button>
    </Stack>
  );
}
