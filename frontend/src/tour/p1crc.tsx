import { Code, List, Text } from "@mantine/core";
import type { Facts, RoiInfo, RunInfo, Tour } from "./types";

/** Guided version of the biological use case (Oliveira et al. 2025, Fig. 4a–b,
 * patient P1CRC); the static page biological-use-case.html has the full text. */

export const DATASET = "Human CRC";
const TUMOR = "P1 tumor";
const PERIPHERY = "P1 periphery 50um";

const same = (a: string | null | undefined, b: string) => (a ?? "").trim().toLowerCase() === b.toLowerCase();

/** The newest RCTD run with Level2 labels (Level2 separates CAF and Macrophage). */
export function rctdRun(facts: Facts): RunInfo | undefined {
  return [...facts.runs].reverse().find((r) => r.method === "spacexr" && r.cell_types.includes("CAF"));
}

export function tumorGroup(facts: Facts): string | undefined {
  return facts.groups.find((g) => same(g, TUMOR));
}

/** The periphery ROI: seeded by the tumor group, with the tumor bins removed. */
export function peripheryRoi(facts: Facts): RoiInfo | undefined {
  const tumor = tumorGroup(facts);
  return tumor ? facts.rois.find((r) => r.seed_group === tumor && r.group_excluded) : undefined;
}

const dominantColumn = (run: RunInfo | undefined) => `deconv_${run?.id ?? "N"}:most_abundant`;
const ui = (label: string) => <Text span fw={600} c="green.3">{label}</Text>;
const RCTD_PATIENCE_MS = 2 * 60 * 1000;

export const P1CRC: Tour = {
  id: "p1crc",
  title: "Colorectal tumor periphery",
  dataset: DATASET,
  steps: [
    {
      id: "intro",
      title: "Which cells surround a colorectal tumor?",
      body: () => (
        <>
          <Text size="sm">
            This tour reproduces, in the app, Figure 4a–b of Oliveira et al. (<i>Nature Genetics</i>, 2025) for
            patient P1CRC: find the tumor, define a 50 µm band around it, look at marker genes along the boundary,
            and compare the cell types in that band with the rest of the tissue.
          </Text>
          <Text size="sm" mt="xs">
            Each step highlights what to use. The tour moves on by itself once you have done it; you can also
            use Next and Back at any time.
          </Text>
        </>
      ),
      readMore: "#overview",
    },
    {
      id: "open-loading",
      title: "Set up a dataset",
      screen: "landing",
      anchor: ["explore"],
      body: () => <Text size="sm">Click {ui("Explore your data")} to open the loading screen.</Text>,
      done: (f) => f.screen !== "landing",
    },
    {
      id: "load",
      title: "Load the P1CRC sample",
      screen: "loading",
      anchor: ["dataset", "pseudospots", "load"],
      body: (f) => (
        <>
          <Text size="sm">
            Choose the {ui(DATASET)} example: 8 µm bins of P1CRC with the study's single-cell reference. Pseudospots
            and cell segmentation are not needed for this figure, so leave them off, then click {ui("Load data")}.
          </Text>
          {f.dataset && f.dataset !== DATASET && (
            <Text size="sm" c="yellow.4" mt="xs">You loaded {f.dataset}; this tour needs {DATASET}.</Text>
          )}
        </>
      ),
      done: (f) => f.dataset === DATASET && (f.loadState === "running" || f.loadState === "ready"),
      readMore: "#load",
    },
    {
      id: "loading",
      title: "Loading",
      anchor: ["progress"],
      body: (f) => f.loadState === "failed"
        ? <Text size="sm" c="red.4">Loading failed (see the message on the page). Go Back to try again.</Text>
        : (
          <Text size="sm">
            The server reads the expression matrix, aligns the bins to the H&amp;E image and builds the image tiles.
            The workspace opens by itself when everything is ready.
          </Text>
        ),
      done: (f) => f.screen === "workspace" && f.loadState === "ready",
    },
    {
      id: "bins",
      title: "Check the bins over the tissue",
      screen: "workspace",
      tool: "display",
      anchor: ["layers"],
      body: () => (
        <Text size="sm">
          Work on the {ui("Bins")} layer (selected in the layer list). Zoom into the map and check that the 8 µm bins
          sit on the H&amp;E tissue. Click Next when you are done looking.
        </Text>
      ),
      done: (f) => f.activeKind === "bins",
      readMore: "#load",
    },
    {
      id: "deconvolution",
      title: "Infer cell types with RCTD",
      screen: "workspace",
      tool: "deconvolution",
      anchor: ["deconv-form", "deconv-run"],
      body: (f) => {
        const spacexr = f.methods.find((m) => m.name === "spacexr");
        return (
          <>
            <List size="sm" spacing={2}>
              <List.Item>Set {ui("Reference cell-type column (.obs)")} to <b>Level2</b> (it separates CAF and Macrophage).</List.Item>
              <List.Item>Choose {ui("spacexr (RCTD)")}; keep {ui("RCTD mode")} = doublet and {ui("Minimum target UMIs")} = 100, as in the paper.</List.Item>
              <List.Item>Click {ui("Run deconvolution")}.</List.Item>
            </List>
            {f.methods.length > 0 && (!spacexr || spacexr.unavailable) && (
              <Text size="sm" c="yellow.4" mt="xs">
                RCTD is not available on this server{spacexr?.unavailable ? `: ${spacexr.unavailable}` : ""}. Run the
                tour on a server with spacexr installed, or locally.
              </Text>
            )}
          </>
        );
      },
      done: (f) => f.deconvJob?.state === "running" || rctdRun(f) !== undefined,
      readMore: "#analysis",
    },
    {
      id: "deconvolution-wait",
      title: "Wait for RCTD",
      screen: "workspace",
      tool: "deconvolution",
      anchor: ["deconv-run", "deconv-status"],
      body: (f) => {
        const job = f.deconvJob;
        const slow = job?.state === "running" && job.started && f.now - job.started > RCTD_PATIENCE_MS;
        return (
          <>
            <Text size="sm">
              With the same data, reference and parameters as the server's saved run, RCTD reuses that result and
              finishes in a few minutes.
            </Text>
            {slow && (
              <Text size="sm" c="yellow.4" mt="xs">
                Still running: there is no saved result for these settings, so RCTD is computing from scratch, which
                takes several hours. Keep this tab open, or check that the column is Level2 and the defaults are unchanged.
              </Text>
            )}
            {job?.state === "failed" && <Text size="sm" c="red.4" mt="xs">The run failed: {job.message}</Text>}
          </>
        );
      },
      done: (f) => rctdRun(f) !== undefined,
    },
    {
      id: "tumor-group",
      title: "Define the tumor",
      screen: "workspace",
      tool: "groups",
      anchor: ["groups-form"],
      body: (f) => (
        <>
          <Text size="sm">The paper uses the <b>Tumor II</b> subtype for P1CRC. Build that group from each bin's dominant cell type:</Text>
          <List size="sm" spacing={2} mt={4}>
            <List.Item>{ui("Gene / metadata")}: <Code>{dominantColumn(rctdRun(f))}</Code></List.Item>
            <List.Item>{ui("Operator")} <Code>==</Code>, {ui("Value")} <Code>Tumor II</Code>, then {ui("Add condition")}</List.Item>
            <List.Item>{ui("Group name")}: <b>{TUMOR}</b>, then {ui("Create group")}</List.Item>
          </List>
        </>
      ),
      done: (f) => tumorGroup(f) !== undefined,
      readMore: "#analysis",
    },
    {
      id: "periphery-roi",
      title: "Draw the 50 µm periphery",
      screen: "workspace",
      tool: "rois",
      anchor: ["roi-form"],
      body: () => (
        <>
          <List size="sm" spacing={2}>
            <List.Item>{ui("ROI name")}: <b>{PERIPHERY}</b></List.Item>
            <List.Item>{ui("Geometry")}: Seed + radius; {ui("Radius (µm)")}: <b>50</b></List.Item>
            <List.Item>{ui("Seed source")}: All bins in a group; {ui("Seed group")}: {TUMOR}</List.Item>
            <List.Item>{ui("ROI region")}: ROI minus group (drops the tumor bins themselves)</List.Item>
            <List.Item>Check the preview, then {ui("Save ROI")}</List.Item>
          </List>
          <Text size="xs" c="dimmed" mt={4}>
            If the radius is shown in px instead of µm, the image calibration is missing.
          </Text>
        </>
      ),
      done: (f) => peripheryRoi(f) !== undefined,
      readMore: "#periphery",
    },
    {
      id: "periphery-group",
      title: "Turn the periphery into a group",
      screen: "workspace",
      tool: "rois",
      anchor: ["roi-list"],
      body: () => (
        <Text size="sm">
          In the saved ROI, click {ui("Create group")}. The band should follow every tumor area across the section,
          including disconnected ones, not a single circle.
        </Text>
      ),
      done: (f) => Boolean(peripheryRoi(f)?.group),
      readMore: "#periphery",
    },
    {
      id: "show-groups",
      title: "Show tumor and periphery",
      screen: "workspace",
      tool: "groups",
      anchor: ["groups-list"],
      body: () => (
        <Text size="sm">
          Colour the map by groups and, to echo the paper's map, make {TUMOR} red and the periphery blue. Look at the
          whole section, then zoom into a boundary. Click Next to continue.
        </Text>
      ),
      readMore: "#markers",
    },
    {
      id: "c1qc",
      title: "Macrophage marker: C1QC",
      screen: "workspace",
      tool: "display",
      anchor: ["gene"],
      body: () => (
        <Text size="sm">
          In {ui("Gene / metadata")}, type and pick <b>C1QC</b>. Is the signal concentrated near the tumor boundary?
        </Text>
      ),
      done: (f) => same(f.colourBy, "C1QC"),
      readMore: "#markers",
    },
    {
      id: "col1a1",
      title: "Fibroblast marker: COL1A1",
      screen: "workspace",
      tool: "display",
      anchor: ["gene"],
      body: () => (
        <Text size="sm">
          Without moving the map, pick <b>COL1A1</b> and compare it with C1QC and the H&amp;E. Neither marker alone
          establishes a cell's identity; together with the deconvolution they support it.
        </Text>
      ),
      done: (f) => same(f.colourBy, "COL1A1"),
      readMore: "#markers",
    },
    {
      id: "composition",
      title: "Compare the periphery with the rest",
      screen: "workspace",
      tool: "metadata",
      anchor: ["metadata-form"],
      body: (f) => (
        <>
          <Text size="sm">
            The comparison uses only the custom groups that exist: delete any other groups first.
          </Text>
          <List size="sm" spacing={2} mt={4}>
            <List.Item>Under {ui("Groups")}, select {peripheryRoi(f)?.group ?? "the periphery group"}, and tick {ui("The rest")}</List.Item>
            <List.Item>{ui("Metadata set")}: the fractions of run {rctdRun(f)?.id ?? "N"}; include <b>CAF</b> and <b>Macrophage</b></List.Item>
            <List.Item>{ui("Summary")}: Mean, then {ui("Compare metadata")}</List.Item>
          </List>
        </>
      ),
      done: (f) => f.results === "metadata",
      readMore: "#composition",
    },
    {
      id: "result",
      title: "Read the result",
      screen: "workspace",
      anchor: ["results"],
      body: () => (
        <>
          <Text size="sm">
            The paper reports <b>CAFs</b> as the most prominent cell type in the periphery and <b>macrophages</b> as its
            most abundant immune population. Hover over the bars for values and counts. Do you see the same here?
          </Text>
          <Text size="xs" c="dimmed" mt="xs">
            These are mean inferred fractions, not the paper's percentages of singlet bins; the use-case page explains
            how to match the published numbers.
          </Text>
        </>
      ),
      readMore: "#results",
    },
  ],
};

export const TOURS: Record<string, Tour> = { [P1CRC.id]: P1CRC };
