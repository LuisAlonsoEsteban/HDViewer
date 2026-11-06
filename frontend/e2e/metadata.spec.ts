import { expect, test } from "@playwright/test";
import { load, roiGroup } from "./helpers";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

test("deconvolution needs a reference; compare numeric metadata between ROI groups", async ({ page }) => {
  test.setTimeout(480_000);
  const errors: string[] = [];
  page.on("console", (message) => { if (message.type() === "error") errors.push(message.text()); });
  page.on("pageerror", (error) => errors.push(String(error)));
  await load(page, DATASET);

  await page.getByRole("tab", { name: "Deconv." }).click();
  await expect(page.getByText(/No single-cell reference was provided/)).toBeVisible();
  await expect(page.getByRole("button", { name: "Run deconvolution" })).toBeDisabled();

  await roiGroup(page, "Left", -150);
  await roiGroup(page, "Right", 150);

  await page.getByRole("tab", { name: "Metadata" }).click();
  await page.getByRole("checkbox", { name: "ROI: Left" }).check();
  await page.getByRole("button", { name: "Compare metadata" }).click();
  await expect(page.getByText(/Select at least two groups or clusters/)).toBeVisible();
  await page.getByRole("checkbox", { name: "ROI: Right" }).check();
  await page.getByRole("combobox", { name: "Metadata columns" }).click();
  const first = page.getByRole("option").first();
  const column = (await first.textContent())!;
  await first.click();
  await page.keyboard.press("Escape");
  await page.getByRole("button", { name: "Compare metadata" }).click();
  await expect(page.getByText("Metadata comparison ready below the map.")).toBeVisible({ timeout: 60_000 });
  const results = page.getByTestId("metadata-results");
  await expect(results.getByText(/Metadata comparison · Bins/)).toBeVisible({ timeout: 30_000 });
  await expect.poll(() => results.locator(".bk-Figure").count()).toBe(1);
  expect(column.length).toBeGreaterThan(0);

  // Against The rest: one chart per selected group
  await page.getByRole("checkbox", { name: /^The rest/ }).check();
  await page.getByRole("button", { name: "Compare metadata" }).click();
  await expect.poll(() => results.locator(".bk-Figure").count(), { timeout: 30_000 }).toBe(2);

  // A new group changes The rest: the comparison is cleared
  await roiGroup(page, "Centre", 0);
  await expect(results.getByText(/metadata comparison was cleared/)).toBeVisible();

  expect(errors.filter((e) => !/favicon/.test(e))).toEqual([]);
});
