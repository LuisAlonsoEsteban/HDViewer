import { expect, Page, test } from "@playwright/test";
import { enableHover, load } from "./helpers";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

const map = (page: Page) => page.evaluate(() => {
  const doc = (window as any).hdv.docs.get("map");
  const plot = doc.get_model_by_name("map");
  const renderer = doc.get_model_by_name("bins:renderer");
  const source = doc.get_model_by_name("bins:source");
  const tiles = doc.get_model_by_name("image:tiles");
  const single = doc.get_model_by_name("image:url");
  return {
    x0: plot.x_range.start, x1: plot.x_range.end,
    lineWidth: JSON.stringify(renderer.glyph.line_width),
    rows: (source.data.line_width ?? []).length,
    imageAlpha: tiles ? tiles.alpha : single?.glyph.global_alpha,
    hover: doc.get_model_by_name("bins:hover").active,
  };
});

test("parameters tab, image opacity and data download", async ({ page }) => {
  test.setTimeout(600_000);
  const errors: string[] = [];
  page.on("console", (message) => {
    // errors, and BokehJS warnings of a source that does not match its glyph
    if (message.type() === "error" || /inconsistent lengths|nonexistent field/.test(message.text())) errors.push(message.text());
  });
  page.on("pageerror", (error) => errors.push(String(error)));
  await load(page, DATASET);

  // Zoom in a little, so the rebuilt figure must keep a non-default view
  const canvas = page.locator(".map-host canvas").first();
  const box = (await canvas.boundingBox())!;
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.wheel(0, -300);
  await page.waitForTimeout(500);
  await enableHover(page);                       // the user's hover choice survives the rebuild
  const before = await map(page);

  // Line width fixed at 3: the figure is rebuilt with a per-row column
  await page.getByRole("tab", { name: "Params" }).first().click();
  const width = page.getByTestId("parameter-line_width");
  await width.getByRole("textbox", { name: "Value", exact: true }).fill("3");
  const rebuilt = page.waitForResponse((r) => r.url().includes("/figures/map?") && r.ok());
  await width.getByRole("button", { name: "Apply" }).click();
  await rebuilt;
  await expect(width.getByText("fixed")).toBeVisible();
  await expect.poll(async () => (await map(page)).rows, { timeout: 60_000 }).toBeGreaterThan(0);
  const after = await map(page);
  expect(after.lineWidth).toContain("line_width");
  expect(after.x0).toBeCloseTo(before.x0, 3);
  expect(after.x1).toBeCloseTo(before.x1, 3);
  expect(after.hover).toBe(true);

  // Reset brings back the default glyph
  await width.getByRole("button", { name: "Reset" }).click();
  await expect(width.getByText("default")).toBeVisible({ timeout: 60_000 });
  expect((await map(page)).lineWidth).not.toContain("field");

  // Image opacity is applied in the browser
  const slider = page.getByRole("slider", { name: "Image opacity" });
  await slider.focus();
  for (let i = 0; i < 10; i++) await page.keyboard.press("ArrowLeft");
  await expect.poll(async () => (await map(page)).imageAlpha).toBeCloseTo(0.5, 2);

  // Download the layer's data
  await page.getByRole("tab", { name: "Display" }).first().click();
  const download = page.waitForEvent("download");
  await page.getByRole("link", { name: "Download .h5ad" }).first().click();
  expect((await download).suggestedFilename()).toBe("bins.h5ad");

  expect(errors).toEqual([]);
});
