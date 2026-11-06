import { expect, Page, test } from "@playwright/test";
import { load } from "./helpers";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

const selected = (page: Page, doc: string, model: string) => page.evaluate(([doc, model]) => {
  const source = (window as any).hdv.docs.get(doc)?.get_model_by_name(model);
  return source ? source.selected.indices.length : -1;
}, [doc, model]);

const select = (page: Page, doc: string, model: string, indices: number[]) => page.evaluate(([doc, model, indices]) => {
  (window as any).hdv.docs.get(doc).get_model_by_name(model).selected.indices = indices;
}, [doc, model, indices] as const);

test("summary, embeddings linked to the map, and comparison of spots with bins", async ({ page }) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on("console", (message) => {
    // errors, and BokehJS warnings of a source that does not match its glyph
    if (message.type() === "error" || /inconsistent lengths|nonexistent field/.test(message.text())) errors.push(message.text());
  });
  page.on("pageerror", (error) => errors.push(String(error)));
  await load(page, DATASET, { pseudospots: true });
  const aside = page.locator(".mantine-AppShell-aside");

  // Summary: dataset table, QC figures, markers before clustering
  await expect(aside.getByRole("cell", { name: "Bins", exact: true })).toBeVisible({ timeout: 60_000 });
  await expect(aside.getByRole("cell", { name: "Spots", exact: true })).toBeVisible();
  await expect.poll(() => aside.locator(".bk-Figure").count(), { timeout: 120_000 }).toBe(4);
  await aside.getByRole("tab", { name: "Markers" }).click();
  await expect(aside.getByText("Run clustering to populate.")).toBeVisible();

  // Embeddings of the bins
  await aside.getByRole("tab", { name: "Embeddings" }).click();
  await aside.getByRole("radio", { name: "Bins" }).check({ force: true });
  await aside.getByRole("button", { name: "Run embedding" }).click();
  await expect(aside.getByText(/^Loaded [\d,]+ bin embeddings\.$/)).toBeVisible({ timeout: 300_000 });
  await expect.poll(() => selected(page, "embedding", "embedding:source"), { timeout: 60_000 }).toBe(0);
  const recoloured = page.waitForResponse((r) => r.url().includes("/embedding/colours") && r.ok());
  await aside.getByRole("combobox", { name: "Color by" }).fill("Ttr");
  await page.getByRole("option", { name: "Ttr", exact: true }).click();
  await recoloured;

  // Embedding -> map, then map (an ROI's members) -> embedding
  const mapSelection = page.waitForResponse((r) => r.url().includes("/embedding/selection") && r.request().method() === "POST");
  await select(page, "embedding", "embedding:source", Array.from({ length: 50 }, (_, i) => i));
  await mapSelection;
  await expect.poll(() => selected(page, "map", "bins:source")).toBeGreaterThan(0);
  await page.getByRole("tab", { name: "ROIs" }).first().click();
  await page.getByRole("button", { name: "Draw" }).click();
  const box = (await page.locator(".map-host canvas").first().boundingBox())!;
  const cx = box.x + box.width / 2, cy = box.y + box.height / 2;
  await page.mouse.move(cx - 60, cy - 60);
  await page.mouse.down();
  await page.mouse.move(cx + 60, cy + 60, { steps: 6 });
  await page.mouse.up();
  await page.getByRole("button", { name: "Save ROI" }).click();
  await expect(page.getByText(/^ROI 1: [\d,]+ observations$/)).toBeVisible({ timeout: 60_000 });
  const embeddingSelection = page.waitForResponse((r) => r.url().includes("/embedding/selection") && r.request().method() === "GET");
  await page.getByTestId("roi").getByRole("button", { name: "Select", exact: true }).click();
  await embeddingSelection;
  await expect.poll(() => selected(page, "embedding", "embedding:source")).toBeGreaterThan(50);

  // Comparison: spots (low-res reference) against bins
  await aside.getByRole("tab", { name: "Comparison" }).click();
  await expect(aside.getByRole("combobox", { name: "Sample 1", exact: true })).toHaveValue("Spots");
  await expect(aside.getByRole("combobox", { name: "Sample 2", exact: true })).toHaveValue("Bins");
  await aside.getByRole("button", { name: "Compare samples" }).click();
  await expect(aside.getByText(/^Compared [\d,]+ spots with [\d,]+ bins\./)).toBeVisible({ timeout: 600_000 });
  await expect(aside.getByRole("cell", { name: /pseudobulk_pearson/ })).toBeVisible();
  await aside.getByRole("cell", { name: /pseudobulk_pearson/ }).click();
  await expect(aside.getByText(/^pseudobulk_pearson: Correlation of each gene/)).toBeVisible();

  await aside.getByRole("tab", { name: "Local" }).click();
  await expect.poll(() => aside.locator(".bk-Figure").count(), { timeout: 60_000 }).toBe(2);
  const shown = page.waitForResponse((r) => r.url().includes("/layers/spots/style") && r.ok());
  await aside.getByRole("button", { name: "Show on map" }).click();
  await shown;
  const spotsSelected = page.waitForResponse((r) => r.url().includes("/comparison/selection") && r.ok());
  await select(page, "comparison", "comparison:source", [0, 1, 2, 3, 4]);
  await spotsSelected;
  await expect.poll(() => selected(page, "map", "spots:source")).toBeGreaterThan(0);

  await aside.getByRole("tab", { name: "Genes" }).click();
  const bothShown = Promise.all(["spots", "bins"].map((layer) =>
    page.waitForResponse((r) => r.url().includes(`/layers/${layer}/style`) && r.ok())));
  await aside.getByRole("tabpanel", { name: "Genes" }).locator("tbody tr").first().click();
  await bothShown;

  expect(errors.filter((e) => !/favicon/.test(e))).toEqual([]);
});
