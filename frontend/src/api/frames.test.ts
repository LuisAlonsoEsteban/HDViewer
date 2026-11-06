import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import { decodeFrame } from "./frames";

// Fixture written by hd_viewer.api.frames.encode (see tests/test_api.py)
const bytes = readFileSync(new URL("./__fixtures__/frame.bin", import.meta.url));

describe("decodeFrame", () => {
  it("reads the Python encoder's header, typed columns, JSON columns and extras", () => {
    const buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
    const { header, ops, extra } = decodeFrame(buffer);
    expect(header.layer).toBe("bins");
    expect(header.constants.width).toBe(2);
    expect(header.palette).toEqual(["#ffffff", "#123456"]);
    const append = ops[0];
    if (append.op !== "append") throw new Error("expected append");
    expect(append.columns.x).toBeInstanceOf(Float32Array);
    expect(Array.from(append.columns.x as Float32Array)).toEqual([1.5, 2.5, 3.5]);
    expect(Array.from(append.columns.color as Uint8Array)).toEqual([1, 0, 1]);
    expect(Array.from(append.columns.codes as Uint16Array)).toEqual([0, 300, 65535]);
    expect(append.columns.x_rings).toEqual([[0, 1], [2], [3, 4, 5]]);
    const set = ops[1];
    if (set.op !== "set") throw new Error("expected set");
    expect(set.name).toBe("alpha");
    expect(Array.from(set.column as Float32Array)).toEqual([0.25, 0.5, 1]);
    expect(Array.from(extra.visible as Uint8Array)).toEqual([1, 0, 1]);
  });
});
