import { expect, test } from "@playwright/test";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

test("load an example, browse the map, hover and recolour", async ({ page }) => {
  const errors: string[] = [];
  page.on("console", (message) => { if (message.type() === "error") errors.push(message.text()); });
  page.on("pageerror", (error) => errors.push(String(error)));

  await page.goto("/");
  await page.getByRole("button", { name: /Explore your data/ }).click();
  await page.getByRole("combobox", { name: "Example dataset" }).click();
  await page.getByRole("option", { name: DATASET }).click();
  await page.getByLabel("Generate Visium v2 pseudospots (55 µm)").uncheck();

  const firstFrame = page.waitForResponse((r) => r.url().includes("/viewport") && r.ok(), { timeout: 200_000 });
  await page.getByRole("button", { name: "Load data" }).click();
  await firstFrame;
  await expect(page.getByRole("button", { name: "Copy viewport coordinates" })).toBeVisible();

  // The map has drawn bins: the layer's BokehJS source holds rows
  await expect.poll(async () => page.evaluate(() => {
    const doc = (window as any).hdv.docs.get("map");
    return doc?.get_model_by_name("bins:source")?.data?.x?.length ?? 0;
  })).toBeGreaterThan(100);

  // Hover a loaded bin: the tooltip asks the server for its barcode
  const target = await page.evaluate(() => {
    const B = (window as any).Bokeh;
    const doc = (window as any).hdv.docs.get("map");
    const plot = doc.get_model_by_name("map");
    const source = doc.get_model_by_name("bins:source");
    const hover = doc.get_model_by_name("bins:hover");
    const view = B.index.find_one(plot);
    const alpha = source.data.alpha;
    let row = 0;
    for (let i = 0; i < alpha.length; i++) if (alpha[i] > 0.3) { row = i; break; }
    const sx = view.frame.x_scale.compute(source.data.x[row]);
    const sy = view.frame.y_scale.compute(source.data.y[row]);
    const rect = view.canvas_view.el.getBoundingClientRect();
    return { x: rect.left + sx, y: rect.top + sy, active: hover.active };
  });
  expect(target.active).toBe(true);
  const names = page.waitForResponse((r) => r.url().includes("/names") && r.ok());
  await page.mouse.move(target.x - 2, target.y - 2);
  await page.mouse.move(target.x, target.y);
  await names;
  await expect(page.locator(".bk-Tooltip, .bk-tooltip").first()).toContainText(/Barcode/);

  // Colour by a gene
  const style = page.waitForResponse((r) => r.url().includes("/style") && r.ok());
  await page.getByRole("combobox", { name: "Gene / metadata" }).fill("Ttr");
  await page.getByRole("option", { name: "Ttr", exact: true }).click();
  await style;

  expect(errors.filter((e) => !/favicon/.test(e))).toEqual([]);
});
