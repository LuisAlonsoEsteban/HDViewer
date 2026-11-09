import { expect, Page, test } from "@playwright/test";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

const binsRows = (page: Page) => page.evaluate(() =>
  (window as any).hdv.docs.get("map")?.get_model_by_name("bins:source")?.data?.x?.length ?? 0);

test("landing page, loading details, help, bin overlay, resume and resizing", async ({ page, context }) => {
  test.setTimeout(600_000);
  await context.grantPermissions(["clipboard-read", "clipboard-write"]);
  const errors: string[] = [];
  page.on("console", (message) => {
    if (message.type() === "error" || /inconsistent lengths|nonexistent field/.test(message.text())) errors.push(message.text());
  });
  page.on("pageerror", (error) => errors.push(String(error)));

  // Landing page: content ported from the old interface
  await page.goto("/");
  await expect(page.getByRole("heading", { name: /See the tissue/ })).toBeVisible();
  for (const logo of ["University of Navarra", "Cima Universidad de Navarra"]) {
    await expect(page.getByRole("img", { name: logo })).toBeVisible();
  }
  await page.getByRole("button", { name: "Spatial bins" }).click();
  await expect(page.getByText("Illustrative spatial grid")).toBeVisible();
  await page.getByRole("button", { name: "Copy citation" }).click();
  await expect(page.getByRole("status")).toHaveText(/Citation (copied|selected)/);
  await expect(page.getByRole("link", { name: /Download citation/ })).toHaveAttribute("href", /citation\.cff$/);
  await expect(page.getByRole("heading", { name: "Bring your own method." })).toBeVisible();
  await expect(page.getByRole("link", { name: /Read the plugin guide/ })).toHaveAttribute("href", /docs\/plugins\.md$/);
  await expect(page.getByRole("link", { name: /Propose a method/ })).toHaveAttribute("href", /method-plugin\.yml$/);
  const icon = await page.locator("link[rel=icon]").getAttribute("href");
  expect((await page.request.get(icon!)).ok()).toBeTruthy();

  // Loading screen: tutorial card, example note, segmentation runner, readable task names
  await page.getByRole("button", { name: /Explore your data/ }).click();
  await expect(page.getByRole("heading", { name: "New here?" })).toBeVisible();
  await page.getByRole("combobox", { name: "Example dataset" }).click();
  await page.getByRole("option", { name: DATASET }).click();
  await expect(page.getByText(/single-cell reference/).first()).toBeVisible();
  await expect(page.getByText(/Runner:/)).toBeVisible({ timeout: 60_000 });
  await page.getByLabel("Generate Visium v2 pseudospots (55 µm)").uncheck();
  const firstFrame = page.waitForResponse((r) => r.url().includes("/viewport") && r.ok(), { timeout: 400_000 });
  await page.getByRole("button", { name: "Load data" }).click();
  await expect(page.getByText("Image tiles")).toBeVisible();
  await firstFrame;
  await expect.poll(() => binsRows(page), { timeout: 60_000 }).toBeGreaterThan(100);

  // Help popover of a tab
  await page.getByRole("button", { name: "Help: Showing a layer" }).first().click();
  await expect(page.getByRole("link", { name: /Read more in the tutorial/ }).first()).toHaveAttribute("href", /tutorial\.html#app-general$/);
  await expect(page.getByRole("link", { name: /Read more in the tutorial/ }).nth(1)).toHaveAttribute("href", /tutorial\.html#app-color$/);
  await page.keyboard.press("Escape");

  // Bin resolution overlay: the level in use, and a fixed native level
  const overlay = page.getByTestId("bin-info-bins");
  await expect(overlay).toContainText(/BINS · \d+ × \d+/);
  await expect(overlay).not.toContainText("Resolution Zoom-driven");        // starts folded to its title
  await overlay.getByRole("button", { name: "Expand bin resolution" }).click();
  await expect(overlay).toContainText("Resolution Zoom-driven");
  await overlay.getByRole("combobox", { name: "Bin resolution" }).click();
  await page.getByRole("option", { name: "Fixed · raw/native" }).click();
  await expect(overlay).toContainText("BINS · 1 × 1 · native", { timeout: 60_000 });
  await expect(overlay).toContainText("Resolution Fixed");
  await overlay.getByRole("button", { name: "Collapse bin resolution" }).click();
  await expect(overlay).not.toContainText("Resolution Fixed");

  // Clear selection (once the native rows have arrived)
  await expect.poll(() => binsRows(page), { timeout: 180_000 }).toBeGreaterThan(100);
  await page.evaluate(() => {
    // every loaded row (some are sequenced; not-sequenced bins have no barcode to select)
    const source = (window as any).hdv.docs.get("map").get_model_by_name("bins:source");
    source.selected.indices = Array.from({ length: source.data.x.length }, (_, i) => i);
  });
  const selectedOnServer = () => page.evaluate(() =>
    fetch("api/layers/bins/selection", { credentials: "same-origin" }).then((r) => r.json()).then((s) => s.selected));
  await expect.poll(selectedOnServer, { timeout: 60_000 }).toBeGreaterThan(0);
  await page.getByRole("button", { name: "Clear selection" }).first().click();
  await expect.poll(selectedOnServer, { timeout: 60_000 }).toBe(0);

  // About -> overview -> Resume workspace: the map is rebuilt and refilled
  await page.getByRole("button", { name: "About", exact: true }).click();
  await page.getByRole("button", { name: /Resume workspace/ }).first().click();
  await expect.poll(() => binsRows(page), { timeout: 60_000 }).toBeGreaterThan(100);

  // Resizable panels (keyboard on the drag handle)
  const handle = page.getByRole("separator", { name: "Resize the layer panel" });
  const before = Number(await handle.getAttribute("aria-valuenow"));
  await handle.focus();
  await page.keyboard.press("Shift+ArrowRight");
  await expect(handle).toHaveAttribute("aria-valuenow", String(before + 50));

  expect(errors).toEqual([]);
});
