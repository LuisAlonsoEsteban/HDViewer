import { useEffect, useState } from "react";
import {
  Autocomplete, Button, Checkbox, ColorInput, Group, NumberInput, Paper, SegmentedControl, Stack, Text,
} from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { rebuildMap } from "../map/layers";
import { Layer } from "../store/app";

interface Spec {
  value: string | number;
  not_sequenced: string | number;
  low?: number;
  high?: number;
  minv?: number;
  maxv?: number;
}

interface Current {
  attr: string;
  mode: "static" | "data";
  value?: string | number;
  column?: string;
  min?: number;
  max?: number;
  not_sequenced: string | number;
  hovered?: boolean;
}

interface Entry {
  attr: string;
  color: boolean;
  spec: Spec;
  current: Current | null;
}

const LABELS: Record<string, string> = {
  fill_alpha: "Fill opacity",
  line_width: "Line width",
  line_color: "Line colour",
  line_alpha: "Line opacity",
};

/** Extra glyph attributes of a layer (the Bokeh UI's Parameters tab): each one
 * follows a fixed value or a column (colours: categories; others: genes or
 * numeric metadata scaled to a range). */
export function ParametersPanel({ layer }: { layer: Layer }) {
  const [entries, setEntries] = useState<Entry[] | null>(null);
  const [error, setError] = useState<string | null>(null);

  const reload = () =>
    getJson<Entry[]>(`layers/${encodeURIComponent(layer.id)}/parameters`).then(setEntries)
      .catch((exc) => setError(exc.message));
  useEffect(() => {
    reload();
  }, [layer.id]);

  if (!entries) return <Text size="sm" c="dimmed">Loading…</Text>;
  return (
    <Stack gap="sm">
      <Text size="xs" c="dimmed">
        Colour and opacity are set in the Display tab; these attributes refine the glyphs' fill and outline.
      </Text>
      {entries.map((entry) => (
        <AttributeEditor key={entry.attr} layer={layer} entry={entry} onChanged={reload} onError={setError} />
      ))}
      {error && <Text size="xs" c="red.4">{error}</Text>}
    </Stack>
  );
}

function AttributeEditor({ layer, entry, onChanged, onError }: {
  layer: Layer; entry: Entry; onChanged(): void; onError(message: string | null): void;
}) {
  const { attr, color, spec, current } = entry;
  const [mode, setMode] = useState<"static" | "data">(current?.mode ?? "static");
  const [value, setValue] = useState<string | number>(current?.value ?? spec.value);
  const [column, setColumn] = useState(current?.column ?? "");
  const [min, setMin] = useState<number | string>(current?.min ?? spec.minv ?? 0);
  const [max, setMax] = useState<number | string>(current?.max ?? spec.maxv ?? 1);
  const [notSequenced, setNotSequenced] = useState<string | number>(current?.not_sequenced ?? spec.not_sequenced);
  const [hovered, setHovered] = useState(current?.hovered ?? false);
  const [options, setOptions] = useState<string[]>([]);
  const [busy, setBusy] = useState(false);
  const path = `layers/${encodeURIComponent(layer.id)}/parameters/${attr}`;

  useEffect(() => {
    if (mode !== "data") return;
    const timer = setTimeout(() => {
      getJson<string[]>(`${path}/columns?q=${encodeURIComponent(column)}`).then(setOptions).catch(() => setOptions([]));
    }, 150);
    return () => clearTimeout(timer);
  }, [mode, column, path]);

  const run = async (action: () => Promise<unknown>) => {
    setBusy(true);
    onError(null);
    try {
      await action();
      await rebuildMap();
      onChanged();
    } catch (exc: any) {
      onError(`${LABELS[attr]}: ${exc.message}`);
    } finally {
      setBusy(false);
    }
  };

  const apply = () => run(() => sendJson("PUT", path, mode === "static"
    ? { mode, value, not_sequenced: notSequenced }
    : { mode, column, min: color ? undefined : Number(min), max: color ? undefined : Number(max),
        not_sequenced: notSequenced, hovered }));
  const reset = () => run(async () => {
    const response = await fetch(`api/${path}`, { method: "DELETE", credentials: "same-origin" });
    if (!response.ok) throw new Error(response.statusText);
  });

  return (
    <Paper withBorder p="xs" data-testid={`parameter-${attr}`}>
      <Stack gap={6}>
        <Group justify="space-between">
          <Text size="sm" fw={600}>{LABELS[attr]}</Text>
          <Text size="xs" c="dimmed">{current ? (current.mode === "static" ? "fixed" : current.column) : "default"}</Text>
        </Group>
        <SegmentedControl size="xs" value={mode} onChange={(m) => setMode(m as "static" | "data")}
                          data={[{ value: "static", label: "Fixed" }, { value: "data", label: "From data" }]} />
        {mode === "static" && (color
          ? <ColorInput size="xs" label="Value" value={String(value)} onChange={setValue} />
          : <NumberInput size="xs" label="Value" value={value} onChange={setValue} min={spec.low} max={spec.high}
                         step={attr.includes("alpha") ? 0.05 : 0.5} decimalScale={2} />)}
        {mode === "data" && (
          <>
            <Autocomplete size="xs" label={color ? "Categorical column" : "Gene or numeric column"}
                          placeholder="Start typing" data={options} value={column} onChange={setColumn} />
            {!color && (
              <Group grow>
                <NumberInput size="xs" label="Min" value={min} onChange={setMin} min={spec.low} max={spec.high}
                             step={0.1} decimalScale={2} />
                <NumberInput size="xs" label="Max" value={max} onChange={setMax} min={spec.low} max={spec.high}
                             step={0.1} decimalScale={2} />
              </Group>
            )}
            <Checkbox size="xs" label="Show the column on hover" checked={hovered}
                      onChange={(e) => setHovered(e.currentTarget.checked)} />
          </>
        )}
        {color
          ? <ColorInput size="xs" label="Not sequenced" value={String(notSequenced)} onChange={setNotSequenced} />
          : <NumberInput size="xs" label="Not sequenced" value={notSequenced} onChange={setNotSequenced}
                         min={spec.low} max={spec.high} step={0.1} decimalScale={2} />}
        <Group gap="xs">
          <Button size="compact-xs" onClick={apply} loading={busy} disabled={mode === "data" && !column}>Apply</Button>
          {current && <Button size="compact-xs" variant="default" onClick={reset} disabled={busy}>Reset</Button>}
        </Group>
      </Stack>
    </Paper>
  );
}
