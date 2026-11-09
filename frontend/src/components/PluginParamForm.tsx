import { Checkbox, NumberInput, Select, Stack } from "@mantine/core";

export interface Param {
  name: string;
  type: "float" | "int" | "bool" | "choice" | "obs_column";
  default: unknown;
  label: string;
  help?: string | null;
  min?: number | null;
  max?: number | null;
  step?: number | null;
  choices?: unknown[] | null;
  required?: boolean;
}

/** Form generated from a plugin manifest's parameters;
 * `columns` are the choices of obs_column parameters (e.g. reference annotations). */
export function PluginParamForm({ params, values, onChange, columns = [], disabled = false }: {
  params: Param[];
  values: Record<string, unknown>;
  onChange(values: Record<string, unknown>): void;
  columns?: string[];
  disabled?: boolean;
}) {
  const set = (name: string, value: unknown) => onChange({ ...values, [name]: value });
  return (
    <Stack gap="xs">
      {params.map((param) => {
        const value = values[param.name] ?? param.default;
        if (param.type === "obs_column")
          return <Select key={param.name} label={param.label} description={param.help ?? undefined} disabled={disabled}
                         data={[{ value: "", label: param.required ? "Choose column" : "None" },
                                ...columns.map((column) => ({ value: column, label: column }))]}
                         value={columns.includes(String(value ?? "")) ? String(value) : ""} allowDeselect={false}
                         onChange={(column) => set(param.name, column ?? "")} />;
        if (param.type === "bool")
          return <Checkbox key={param.name} label={param.label} description={param.help ?? undefined} disabled={disabled}
                           checked={Boolean(value)} onChange={(e) => set(param.name, e.currentTarget.checked)} />;
        if (param.type === "choice")
          return <Select key={param.name} label={param.label} description={param.help ?? undefined} disabled={disabled}
                         data={(param.choices ?? []).map(String)} value={String(value)}
                         onChange={(choice) => set(param.name, (param.choices ?? []).find((c) => String(c) === choice))} />;
        return <NumberInput key={param.name} label={param.label} description={param.help ?? undefined} disabled={disabled}
                            value={value as number | undefined} min={param.min ?? undefined} max={param.max ?? undefined}
                            step={param.step ?? (param.type === "int" ? 1 : 0.1)} allowDecimal={param.type === "float"}
                            placeholder="default" onChange={(v) => set(param.name, v === "" ? null : v)} />;
      })}
    </Stack>
  );
}
