import { expect, Page, test } from "@playwright/test";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

async function mapBox(page: Page) {
  const box = await page.locator(".map-host canvas").first().boundingBox();
  if (!box) throw new Error("no map canvas");
  return { cx: box.x + box.width / 2, cy: box.y + box.height / 2 };
}

const roiData = (page: Page, name: string) => page.evaluate((name) => {
  const doc = (window as any).hdv.docs.get("map");
  return doc.get_model_by_name(name).data.xs.length;
}, name);

test("draw ROIs, group one and run differential expression on it", async ({ page }) => {
  test.setTimeout(480_000);
  const errors: string[] = [];
  page.on("console", (message) => { if (message.type() === "error") errors.push(message.text()); });
  page.on("pageerror", (error) => errors.push(String(error)));

  await page.goto("/");
  await page.getByRole("button", { name: "Launch" }).click();
  await page.getByRole("combobox", { name: "Example dataset" }).click();
  await page.getByRole("option", { name: DATASET }).click();
  await page.getByLabel("Generate Visium v2 pseudospots (55 µm)").uncheck();
  const firstFrame = page.waitForResponse((r) => r.url().includes("/viewport") && r.ok(), { timeout: 200_000 });
  await page.getByRole("button", { name: "Load data" }).click();
  await firstFrame;

  // Bounding box dragged on the map
  await page.getByRole("tab", { name: "ROIs" }).click();
  await page.getByRole("textbox", { name: "ROI name (optional)" }).fill("Centre");
  await page.getByRole("button", { name: "Draw" }).click();
  await expect(page.getByText(/drag the region on the map/)).toBeVisible();
  const { cx, cy } = await mapBox(page);
  await page.mouse.move(cx - 80, cy - 60);
  await page.mouse.down();
  await page.mouse.move(cx, cy, { steps: 5 });
  await page.mouse.move(cx + 80, cy + 60, { steps: 5 });
  await page.mouse.up();
  await expect(page.getByText("Draft ready. Save ROI or cancel.")).toBeVisible({ timeout: 30_000 });
  expect(await roiData(page, "draft")).toBe(1);
  // Drawing an ROI must not select bins
  expect(await page.evaluate(() => (window as any).hdv.model("bins", "source").selected.indices.length)).toBe(0);
  await page.getByRole("button", { name: "Save ROI" }).click();
  await expect(page.getByText(/^Centre: [\d,]+ observations$/)).toBeVisible({ timeout: 60_000 });
  expect(await roiData(page, "draft")).toBe(0);
  expect(await roiData(page, "bins:rois")).toBe(1);

  // Seed + radius clicked on the map
  await page.getByRole("combobox", { name: "Geometry" }).click();
  await page.getByRole("option", { name: "Seed + radius" }).click();
  await page.getByRole("button", { name: "Draw" }).click();
  await page.mouse.click(cx + 150, cy + 100);
  await expect(page.getByText("Draft ready. Save ROI or cancel.")).toBeVisible({ timeout: 30_000 });
  await page.getByRole("button", { name: "Save ROI" }).click();
  await expect(page.getByText(/^ROI \d+: [\d,]+ observations$/)).toBeVisible({ timeout: 60_000 });
  expect(await roiData(page, "bins:rois")).toBe(2);

  // Group, select and outline toggling of the first ROI
  const centre = page.getByTestId("roi").filter({ hasText: "Centre" });
  await centre.getByRole("button", { name: "Create group" }).click();
  await expect(page.getByText("Created group: ROI: Centre")).toBeVisible();
  await expect(centre.getByRole("button", { name: "Group created" })).toBeDisabled();
  await centre.getByRole("button", { name: "Select", exact: true }).click();
  await expect(page.getByText(/[\d,]+ observations selected\./)).toBeVisible();
  await expect.poll(() => page.evaluate(() => (window as any).hdv.model("bins", "source").selected.indices.length))
    .toBeGreaterThan(0);
  await centre.getByRole("checkbox", { name: "Outline" }).uncheck();
  await expect.poll(() => roiData(page, "bins:rois")).toBe(1);

  // DE: the ROI's group against the rest
  await page.getByRole("tab", { name: "DE", exact: true }).click();
  await page.getByRole("checkbox", { name: "ROI: Centre" }).check();
  await page.getByRole("button", { name: "Run DE" }).click();
  const results = page.getByTestId("de-results");
  await expect(results.getByText("DE vs Rest")).toBeVisible({ timeout: 300_000 });
  await expect(results.locator("tbody tr").first()).toBeVisible();
  await expect.poll(() => results.locator(".bk-Figure, .bk-Column").count(), { timeout: 30_000 }).toBeGreaterThan(0);
  const sorted = page.waitForResponse((r) => r.url().includes("/de/table") && r.url().includes("desc=true") && r.ok());
  await results.getByRole("button", { name: /^L2FC/ }).click();          // descending first
  await sorted;
  await expect(results.getByRole("button", { name: /^L2FC ▼/ })).toBeVisible();
  const lfcs = await results.locator("tbody tr td:nth-child(2)").allTextContents();
  const values = lfcs.map(Number).filter((v) => !Number.isNaN(v));
  expect(values).toEqual([...values].sort((a, b) => b - a));
  await results.getByPlaceholder("Search gene...").fill("Ttr");
  await expect(results.getByText(/genes match/)).toBeVisible();

  // Deleting the compared ROI deletes its group and clears the results
  await page.getByRole("tab", { name: "ROIs" }).click();
  await centre.getByRole("button", { name: "Delete Centre" }).click();
  await expect(page.getByText("ROI and group 'ROI: Centre' deleted.")).toBeVisible();
  await expect(results.getByText(/DE results were cleared/)).toBeVisible();

  expect(errors.filter((e) => !/favicon/.test(e))).toEqual([]);
});
