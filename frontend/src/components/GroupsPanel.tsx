import { useEffect, useState } from "react";
import {
  ActionIcon, Autocomplete, Button, Checkbox, ColorInput, Group, Select, Stack, Text, TextInput,
} from "@mantine/core";
import { getJson, sendJson } from "../api/client";
import { patchFrame, refreshLayer } from "../map/layers";
import { Layer, useApp } from "../store/app";

interface GroupInfo { name: string; colour: string; size: number; visible: boolean; conditions: string[][] }
interface Groups { groups: GroupInfo[]; operators: string[]; shown: boolean }
interface Condition { column: string; op: string; value: string; logic: string }

/** The Bokeh UI's Groups tab: build conditions, create groups, show/hide them. */
export function GroupsPanel({ layer }: { layer: Layer }) {
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const analysesChanged = useApp((s) => s.analysesChanged);   // DE and ROI tabs list the groups
  const [state, setState] = useState<Groups>();
  const [column, setColumn] = useState("");
  const [options, setOptions] = useState<string[]>([]);
  const [op, setOp] = useState<string | null>(">");
  const [value, setValue] = useState("");
  const [logic, setLogic] = useState<string | null>("AND");
  const [conditions, setConditions] = useState<Condition[]>([]);
  const [name, setName] = useState("Group A");
  const [message, setMessage] = useState<{ text: string; error: boolean } | null>(null);
  const [busy, setBusy] = useState(false);
  const base = `layers/${encodeURIComponent(layer.id)}`;

  useEffect(() => { getJson<Groups>(`${base}/groups`).then(setState); }, [version, base]);
  useEffect(() => {
    const timer = setTimeout(() => getJson<string[]>(`${base}/columns?q=${encodeURIComponent(column)}`)
      .then(setOptions).catch(() => setOptions([])), 150);
    return () => clearTimeout(timer);
  }, [column, base]);

  const addCondition = async () => {
    try {
      const { hint } = await sendJson<{ hint: string }>("POST", `${base}/groups/check`, { column, op, value });
      setConditions((c) => [...c, { column, op: op!, value, logic: c.length ? logic! : "AND" }]);
      setMessage(hint ? { text: hint, error: false } : null);
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    }
  };

  const create = async () => {
    setBusy(true);
    try {
      setState(await sendJson<Groups>("POST", `${base}/groups`, { name, conditions }));
      setConditions([]);
      setMessage({ text: `Group '${name}' created.`, error: false });
      await refreshLayer(layer.id);
      analysesChanged(layer.id);
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    } finally {
      setBusy(false);
    }
  };

  const display = async (group: string, changes: Record<string, unknown>) => {
    setState((s) => s && { ...s, groups: s.groups.map((g) => (g.name === group ? { ...g, ...changes } : g)) });
    await patchFrame(`${base}/groups/${encodeURIComponent(group)}`, changes);
    setState(await getJson<Groups>(`${base}/groups`));
  };
  const remove = async (group: string) => {
    setState(await sendJson<Groups>("DELETE", `${base}/groups/${encodeURIComponent(group)}`));
    await refreshLayer(layer.id);
    analysesChanged(layer.id);
  };

  return (
    <Stack gap="sm">
      <Text fw={600} size="sm">Build conditions</Text>
      <Autocomplete label="Gene / metadata" placeholder="Start typing" data={options} value={column} onChange={setColumn} />
      <Group grow>
        <Select label="Operator" data={state?.operators ?? [">", "<", ">=", "<=", "=="]} value={op} onChange={setOp} />
        <TextInput label="Value" value={value} onChange={(e) => setValue(e.currentTarget.value)} />
      </Group>
      {conditions.length > 0 && <Select label="Combine with previous" data={["AND", "OR"]} value={logic} onChange={setLogic} />}
      <Group>
        <Button size="xs" variant="default" onClick={addCondition} disabled={!column || !value}>Add condition</Button>
        <Button size="xs" variant="subtle" onClick={() => setConditions([])} disabled={!conditions.length}>Clear</Button>
      </Group>
      {conditions.map((c, i) => (
        <Text key={i} size="xs" ff="monospace">{i ? `${c.logic} ` : ""}{c.column} {c.op} {c.value}</Text>
      ))}
      <TextInput label="Group name" value={name} onChange={(e) => setName(e.currentTarget.value)} />
      <Button onClick={create} loading={busy} disabled={!conditions.length}>Create group</Button>
      {message && <Text size="xs" c={message.error ? "red.4" : "teal.3"}>{message.text}</Text>}

      {state && state.groups.length > 0 && (
        <Stack gap={6}>
          <Text fw={600} size="sm">Created groups</Text>
          {!state.shown && <Button size="xs" variant="light" onClick={async () => {
            await sendJson("POST", `${base}/groups/show`); await refreshLayer(layer.id);
            setState(await getJson<Groups>(`${base}/groups`));
          }}>Colour the map by groups</Button>}
          {state.groups.map((group) => (
            <Group key={group.name} gap="xs" wrap="nowrap">
              <Checkbox checked={group.visible} onChange={(e) => display(group.name, { visible: e.currentTarget.checked })}
                        label={<Text size="sm" c={group.colour} fw={600}>{group.name} <Text span size="xs" c="dark.2">({group.size})</Text></Text>} />
              <ColorInput size="xs" w={110} value={group.colour} onChangeEnd={(colour) => display(group.name, { colour })} />
              <ActionIcon variant="subtle" color="red" aria-label={`Delete ${group.name}`} onClick={() => remove(group.name)}>✕</ActionIcon>
            </Group>
          ))}
        </Stack>
      )}
    </Stack>
  );
}
