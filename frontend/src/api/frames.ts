// Binary frames from the API (see hd_viewer/api/frames.py): uint32 header
// length, JSON header, then 8-byte aligned numeric buffers viewed in place.

export type Column = ArrayLike<unknown>;
export type Op =
  | { op: "append" | "replace"; columns: Record<string, Column> }
  | { op: "set"; name: string; column: Column };

export interface FrameHeader {
  layer: string;
  base: number;
  version: number;
  row_version: number;
  rows: number;
  level: string | null;
  constants: Record<string, number>;
  palette: string[];
  palettes?: Record<string, string[]>;   // Parameters tab colour attributes
  tooltips: [string, string][];
  visible: boolean;
  filtered: boolean;
}

export interface Frame {
  header: FrameHeader;
  ops: Op[];
  extra: Record<string, Column>;
}

type Encoded = { dtype: string; offset: number; length: number } | { values: unknown[] };

const ARRAYS: Record<string, (b: ArrayBuffer, o: number, n: number) => ArrayLike<number>> = {
  float32: (b, o, n) => new Float32Array(b, o, n),
  float64: (b, o, n) => new Float64Array(b, o, n),
  uint8: (b, o, n) => new Uint8Array(b, o, n),
  bool: (b, o, n) => new Uint8Array(b, o, n),
  uint16: (b, o, n) => new Uint16Array(b, o, n),
  uint32: (b, o, n) => new Uint32Array(b, o, n),
  int8: (b, o, n) => new Int8Array(b, o, n),
  int16: (b, o, n) => new Int16Array(b, o, n),
  int32: (b, o, n) => new Int32Array(b, o, n),
};

export function decodeFrame(buffer: ArrayBuffer): Frame {
  const length = new DataView(buffer).getUint32(0, true);
  const raw = JSON.parse(new TextDecoder().decode(new Uint8Array(buffer, 4, length)));
  const start = 4 + length;
  const body = buffer.slice(start); // own the bytes so views start 8-byte aligned
  const resolve = (column: Encoded): Column => {
    if ("values" in column) return column.values;
    const make = ARRAYS[column.dtype];
    if (!make) throw new Error(`Unsupported dtype ${column.dtype}`);
    return make(body, column.offset, column.length);
  };
  const ops: Op[] = (raw.ops ?? []).map((op: any) =>
    op.op === "set"
      ? { op: "set", name: op.name, column: resolve(op.column) }
      : { op: op.op, columns: Object.fromEntries(Object.entries(op.columns).map(([k, v]) => [k, resolve(v as Encoded)])) },
  );
  const extra = Object.fromEntries(Object.entries(raw.extra ?? {}).map(([k, v]) => [k, resolve(v as Encoded)]));
  delete raw.ops;
  delete raw.extra;
  return { header: raw as FrameHeader, ops, extra };
}
