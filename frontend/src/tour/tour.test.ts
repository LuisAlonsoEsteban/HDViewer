import { beforeEach, describe, expect, it } from "vitest";
import { P1CRC, peripheryRoi, rctdRun, tumorGroup } from "./p1crc";
import { useTour } from "./store";
import type { Facts } from "./types";

const facts = (changes: Partial<Facts> = {}): Facts => ({
  screen: "workspace", loadState: "ready", dataset: "Human CRC", datasets: ["Human CRC"], hasWorkspace: true,
  binsLayer: "bins", activeKind: "bins", colourBy: null, deconvJob: undefined, methods: [], runs: [],
  groups: [], rois: [], results: null, now: 0, ...changes,
});
const step = (id: string) => P1CRC.steps.find((s) => s.id === id)!;
const level2 = { id: "2", method: "spacexr", cell_types: ["CAF", "Macrophage", "Tumor II"], fraction_columns: ["deconv_2:CAF"] };

describe("tour store", () => {
  beforeEach(() => useTour.getState().stop());

  it("starts, moves and stops", () => {
    const tour = useTour.getState();
    tour.start("p1crc");
    useTour.getState().next();
    useTour.getState().next();
    useTour.getState().back();
    expect(useTour.getState()).toMatchObject({ tour: "p1crc", step: 1 });
    useTour.getState().back();
    useTour.getState().back();
    expect(useTour.getState().step).toBe(0);
    useTour.getState().stop();
    expect(useTour.getState().tour).toBeNull();
  });
});

describe("P1CRC steps", () => {
  it("has unique step ids", () => {
    const ids = P1CRC.steps.map((s) => s.id);
    expect(new Set(ids).size).toBe(ids.length);
  });

  it("needs the Human CRC example to be loading", () => {
    expect(step("load").done!(facts({ screen: "loading", loadState: "running", dataset: "Mouse brain (8um)" }))).toBe(false);
    expect(step("load").done!(facts({ screen: "loading", loadState: "running" }))).toBe(true);
  });

  it("finds the Level2 RCTD run, not an NNLS or Level1 run", () => {
    const nnls = { ...level2, id: "3", method: "nnls" };
    const level1 = { ...level2, id: "1", cell_types: ["Fibroblast", "Myeloid", "Tumor"] };
    expect(rctdRun(facts({ runs: [level1, level2, nnls] }))?.id).toBe("2");
    expect(step("deconvolution-wait").done!(facts({ runs: [level1] }))).toBe(false);
    expect(step("deconvolution-wait").done!(facts({ runs: [level2] }))).toBe(true);
  });

  it("follows the tumor group into the periphery ROI and its group", () => {
    const roi = { name: "P1 periphery 50um", group: null, seed_group: "P1 tumor", group_excluded: true };
    expect(tumorGroup(facts({ groups: [" p1 Tumor "] }))).toBe(" p1 Tumor ");
    expect(peripheryRoi(facts({ groups: ["P1 tumor"], rois: [{ ...roi, group_excluded: false }] }))).toBeUndefined();
    expect(step("periphery-roi").done!(facts({ groups: ["P1 tumor"], rois: [roi] }))).toBe(true);
    expect(step("periphery-group").done!(facts({ groups: ["P1 tumor"], rois: [roi] }))).toBe(false);
    expect(step("periphery-group").done!(facts({ groups: ["P1 tumor"], rois: [{ ...roi, group: roi.name }] }))).toBe(true);
  });

  it("checks the displayed genes and the metadata comparison", () => {
    expect(step("c1qc").done!(facts({ colourBy: "C1QC" }))).toBe(true);
    expect(step("col1a1").done!(facts({ colourBy: "C1QC" }))).toBe(false);
    expect(step("composition").done!(facts({ results: "metadata" }))).toBe(true);
  });
});
