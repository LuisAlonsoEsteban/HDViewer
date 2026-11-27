import { useEffect, useState } from "react";
import {
  ActionIcon, Autocomplete, Button, Checkbox, ColorInput, Group, SegmentedControl, Select, Stack, Text, TextInput,
} from "@mantine/core";
import { frame, getJson, sendJson } from "../api/client";
import { applyLayerFrame, patchFrame, refreshLayer, requestViewport } from "../map/layers";
import { Layer, useApp } from "../store/app";
import { notifyError, reportError } from "../notify";
import { create as createStore } from "zustand";

interface GroupInfo {
  name: string; colour: string; size: number; visible: boolean; conditions: string[][];
  description: string | null;     // how a combined group was made, e.g. "A ∩ B"
}
interface Groups { groups: GroupInfo[]; operators: string[]; shown: boolean }
interface Condition { column: string; op: string; value: string; logic: string }

/** The group being built on each layer and whether its creation is running. Kept
 * outside the panel: switching tool tabs unmounts it while the request goes on. */
interface Draft { conditions: Condition[]; name: string; message: { text: string; error: boolean } | null; busy: boolean }
const EMPTY: Draft = { conditions: [], name: "Group A", message: null, busy: false };
const useDrafts = createStore<{ drafts: Record<string, Draft>; update(layer: string, changes: Partial<Draft>): void }>(
  (set) => ({
    drafts: {},
    update: (layer, changes) => set((s) => ({ drafts: { ...s.drafts, [layer]: { ...(s.drafts[layer] ?? EMPTY), ...changes } } })),
  }));
type Update<T> = T | ((current: T) => T);

/** Groups tool: build conditions, create groups, show/hide them. */
export function GroupsPanel({ layer }: { layer: Layer }) {
  const version = useApp((s) => s.analysisVersion[layer.id] ?? 0);
  const analysesChanged = useApp((s) => s.analysesChanged);   // DE and ROI tabs list the groups
  const [state, setState] = useState<Groups>();
  const [column, setColumn] = useState("");
  const [options, setOptions] = useState<string[]>([]);
  const [op, setOp] = useState<string | null>(">");
  const [value, setValue] = useState("");
  const [logic, setLogic] = useState<string | null>("AND");
  const { conditions, name, message, busy } = useDrafts((s) => s.drafts[layer.id] ?? EMPTY);
  const draft = (changes: Partial<Draft>) => useDrafts.getState().update(layer.id, changes);
  const setConditions = (next: Update<Condition[]>) => draft({
    conditions: typeof next === "function" ? next(useDrafts.getState().drafts[layer.id]?.conditions ?? []) : next });
  const setName = (value: string) => draft({ name: value });
  const setMessage = (value: Draft["message"]) => draft({ message: value });
  const base = `layers/${encodeURIComponent(layer.id)}`;

  const reload = () => getJson<Groups>(`${base}/groups`).then(setState).catch(reportError("Groups could not be loaded"));
  useEffect(() => { reload(); }, [version, base]);
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
    const id = layer.id;              // the panel may unmount (another tab) before the reply
    const update = (changes: Partial<Draft>) => useDrafts.getState().update(id, changes);
    update({ busy: true, message: { text: `Creating group '${name}'…`, error: false } });
    try {
      setState(await sendJson<Groups>("POST", `${base}/groups`, { name, conditions }));
      update({ conditions: [], message: { text: `Group '${name}' created.`, error: false } });
      await refreshLayer(id);
      analysesChanged(id);            // the remounted panel (and DE/ROI tabs) reload the groups
    } catch (exc: any) {
      update({ message: { text: exc.message, error: true } });
    } finally {
      update({ busy: false });
    }
  };

  const display = async (group: string, changes: Record<string, unknown>) => {
    setState((s) => s && { ...s, groups: s.groups.map((g) => (g.name === group ? { ...g, ...changes } : g)) });
    try {
      await patchFrame(`${base}/groups/${encodeURIComponent(group)}`, changes);
    } catch (exc) {
      notifyError(exc, "The group could not be updated");
    }
    await reload();                // the server's state, also undoing a failed change
  };
  const remove = async (group: string) => {
    try {
      setState(await sendJson<Groups>("DELETE", `${base}/groups/${encodeURIComponent(group)}`));
      await refreshLayer(layer.id);
      analysesChanged(layer.id);
    } catch (exc) {
      notifyError(exc, "The group could not be deleted");
    }
  };

  return (
    <Stack gap="sm">
      <Text fw={600} size="sm">Build conditions</Text>
      <Autocomplete label="Gene / metadata" placeholder="Start typing" data={options} value={column} onChange={setColumn}
                    wrapperProps={{ "data-tour": "groups-form" }} />
      <Group grow data-tour="groups-form">
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
      <TextInput label="Group name" value={name} onChange={(e) => setName(e.currentTarget.value)}
                 wrapperProps={{ "data-tour": "groups-form" }} />
      <Button onClick={create} loading={busy} disabled={!conditions.length} data-tour="groups-form">Create group</Button>
      {message && <Text size="xs" c={message.error ? "red.4" : "teal.3"}>{message.text}</Text>}

      {state && state.groups.length > 0 && (
        <Stack gap={6} data-tour="groups-list">
          <Text fw={600} size="sm">Created groups</Text>
          {!state.shown && <Button size="xs" variant="light" onClick={async () => {
            try {
              // The reply is the layer's recoloured map frame (binary), not JSON
              await applyLayerFrame(await frame("POST", `${base}/groups/show`));
              useApp.getState().updateLayer(layer.id, { colour_by: { kind: "groups", name: "custom_groups" } });
              await requestViewport(layer.id);
            } catch (exc) {
              notifyError(exc, "The map could not be coloured by groups");
            }
            await reload();
          }}>Colour the map by groups</Button>}
          {state.groups.map((group) => (
            <Group key={group.name} gap="xs" wrap="nowrap">
              <Checkbox checked={group.visible} onChange={(e) => display(group.name, { visible: e.currentTarget.checked })}
                        label={<Text size="sm" c={group.colour} fw={600}>{group.name} <Text span size="xs" c="dark.2">({group.size})</Text></Text>}
                        description={group.description && group.description !== group.name ? group.description : undefined}
                        style={{ flex: 1, minWidth: 0 }} />
              <ColorInput size="xs" w={110} value={group.colour} onChangeEnd={(colour) => display(group.name, { colour })} />
              <ActionIcon variant="subtle" color="red" aria-label={`Delete ${group.name}`} onClick={() => remove(group.name)}>✕</ActionIcon>
            </Group>
          ))}
        </Stack>
      )}

      {state && state.groups.length > 0 && (
        <CombineGroups base={base} groups={state.groups.map((g) => g.name)}
                       onCreated={async (next) => {
                         setState(next);
                         await refreshLayer(layer.id);
                         analysesChanged(layer.id);
                       }} />
      )}
    </Stack>
  );
}

const OPERATIONS = [
  { value: "intersection", label: "A ∩ B", help: "In both groups" },
  { value: "union", label: "A ∪ B", help: "In either group" },
  { value: "difference", label: "A − B", help: "In A but not in B" },
  { value: "complement", label: "not A", help: "Every sequenced observation not in A" },
];

/** A new group from a set operation on existing groups, or the opposite set. */
function CombineGroups({ base, groups, onCreated }: {
  base: string; groups: string[]; onCreated(state: Groups): Promise<void>;
}) {
  const [operation, setOperation] = useState("intersection");
  const [first, setFirst] = useState<string | null>(groups[0] ?? null);
  const [second, setSecond] = useState<string | null>(groups[1] ?? null);
  const [opposite, setOpposite] = useState(false);
  const [name, setName] = useState("");
  const [preview, setPreview] = useState<{ size: number; description: string; total: number } | null>(null);
  const [message, setMessage] = useState<{ text: string; error: boolean } | null>(null);
  const [busy, setBusy] = useState(false);
  const unary = operation === "complement";
  const body = { operation, first, second: unary ? null : second, opposite: !unary && opposite };
  const ready = Boolean(first && (unary || (second && second !== first)));

  // Keep the choices valid as groups are created and deleted
  useEffect(() => {
    if (!first || !groups.includes(first)) setFirst(groups[0] ?? null);
    if (!second || !groups.includes(second)) setSecond(groups.find((g) => g !== first) ?? null);
  }, [groups.join("|")]);

  useEffect(() => {
    setPreview(null);
    if (!ready) return;
    const timer = setTimeout(() => {
      sendJson<{ size: number; description: string; total: number }>("POST", `${base}/groups/combine`, { ...body, preview: true })
        .then(setPreview).catch(() => setPreview(null));
    }, 150);
    return () => clearTimeout(timer);
  }, [base, operation, first, second, opposite, ready]);

  const create = async () => {
    setBusy(true);
    setMessage(null);
    try {
      const reply = await sendJson<Groups & { group: string; size: number }>("POST", `${base}/groups/combine`, { ...body, name });
      setMessage({ text: `Group '${reply.group}' created (${reply.size.toLocaleString()} observations).`, error: false });
      setName("");
      await onCreated(reply);
    } catch (exc: any) {
      setMessage({ text: exc.message, error: true });
    } finally {
      setBusy(false);
    }
  };

  return (
    <Stack gap={8} data-testid="combine-groups">
      <Text fw={600} size="sm">Combine groups</Text>
      <SegmentedControl size="xs" value={operation} onChange={setOperation} fullWidth
                        data={OPERATIONS.map(({ value, label }) => ({ value, label }))} aria-label="Set operation" />
      <Text size="xs" c="dimmed">{OPERATIONS.find((o) => o.value === operation)?.help}</Text>
      <Group grow align="end">
        <Select size="xs" label="A" data={groups} value={first} onChange={setFirst} allowDeselect={false} />
        {!unary && (
          <Select size="xs" label="B" data={groups.filter((g) => g !== first)} value={second} onChange={setSecond}
                  allowDeselect={false} />
        )}
      </Group>
      {!unary && (
        <Checkbox size="xs" label="Take the opposite set (everything outside the result)" checked={opposite}
                  onChange={(e) => setOpposite(e.currentTarget.checked)} />
      )}
      <TextInput size="xs" label="Name (optional)" placeholder={preview?.description ?? "e.g. A ∩ B"} value={name}
                 onChange={(e) => setName(e.currentTarget.value)} />
      {preview && (
        <Text size="xs" c={preview.size ? "dark.1" : "yellow.4"} role="status">
          <b>{preview.description}</b>: {preview.size.toLocaleString()} of {preview.total.toLocaleString()} observations
          {preview.size ? "" : " (empty)"}
        </Text>
      )}
      <Button size="xs" onClick={create} loading={busy} disabled={!ready || !preview?.size}>Create combined group</Button>
      {message && <Text size="xs" c={message.error ? "red.4" : "teal.3"}>{message.text}</Text>}
      <Text size="xs" c="dimmed">The new group keeps its members even if A or B change later.</Text>
    </Stack>
  );
}
