import { expect, Page } from "@playwright/test";

export async function load(page: Page, dataset: string, { pseudospots = false } = {}) {
  await page.goto("/");
  await page.getByRole("button", { name: /Explore your data/ }).click();
  await page.getByRole("combobox", { name: "Example dataset" }).click();
  await page.getByRole("option", { name: dataset }).click();
  await page.getByLabel("Generate Visium v2 pseudospots (55 µm)").setChecked(pseudospots);
  const firstFrame = page.waitForResponse((r) => r.url().includes("/viewport") && r.ok(), { timeout: 600_000 });
  await page.getByRole("button", { name: "Load data" }).click();
  await firstFrame;
}

/** Draw a box ROI (offsets from the map centre, in screen pixels) and make it a group. */
export async function roiGroup(page: Page, name: string, dx: number) {
  await page.getByRole("tab", { name: "ROIs" }).click();
  await page.getByRole("textbox", { name: "ROI name (optional)" }).fill(name);
  await page.getByRole("button", { name: "Draw" }).click();
  const box = (await page.locator(".map-host canvas").first().boundingBox())!;
  const cx = box.x + box.width / 2 + dx, cy = box.y + box.height / 2;
  await page.mouse.move(cx - 50, cy - 50);
  await page.mouse.down();
  await page.mouse.move(cx, cy, { steps: 4 });
  await page.mouse.move(cx + 50, cy + 50, { steps: 4 });
  await page.mouse.up();
  await expect(page.getByText("Draft ready. Save ROI or cancel.")).toBeVisible({ timeout: 30_000 });
  await page.getByRole("button", { name: "Save ROI" }).click();
  await expect(page.getByText(new RegExp(`^${name}: [\\d,]+ observations$`))).toBeVisible({ timeout: 60_000 });
  await page.getByTestId("roi").filter({ hasText: name }).getByRole("button", { name: "Create group" }).click();
  await expect(page.getByText(`Created group: ROI: ${name}`)).toBeVisible({ timeout: 60_000 });
}

/** Turn on the map's hover tooltips (off by default) with the toolbar button. */
export async function enableHover(page: Page): Promise<void> {
  await page.locator(".map-host .bk-tool-icon-hover").first().click();
  await expect.poll(() => page.evaluate(() =>
    (window as any).hdv.docs.get("map").get_model_by_name("bins:hover").active)).toBe(true);
}
