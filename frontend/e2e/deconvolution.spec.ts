import { expect, test } from "@playwright/test";
import { load, roiGroup } from "./helpers";

// A dataset with a single-cell reference; NNLS runs on every native bin (slow
// unless the result is cached under the dataset's deconvolution folder).
const DATASET = process.env.HDV_DECONV_DATASET ?? "Human CRC";

test("deconvolve with NNLS, show a cell type and compare fractions between groups", async ({ page }) => {
  test.setTimeout(7_200_000);                 // NNLS on a whole 8 µm slide: ~80 min
  const errors: string[] = [];
  page.on("console", (message) => { if (message.type() === "error") errors.push(message.text()); });
  page.on("pageerror", (error) => errors.push(String(error)));
  await load(page, DATASET);

  await page.getByRole("tab", { name: "Deconv." }).click();
  await expect(page.getByText(/^Reference ready: /)).toBeVisible({ timeout: 60_000 });
  await expect(page.getByRole("combobox", { name: "Deconvolution method" })).toHaveValue("NNLS");
  const run = page.getByRole("button", { name: "Run deconvolution" });
  await expect(run).toBeEnabled({ timeout: 30_000 });
  const recoloured = page.waitForResponse((r) => r.url().includes("/style") && r.request().method() === "PATCH" && r.ok(),
                                          { timeout: 6_900_000 });
  await run.click();
  await expect(page.getByText(/^Deconvolved [\d,]+ observations into \d+ cell types$/)).toBeVisible({ timeout: 6_900_000 });
  await recoloured;                                         // the run's first cell type is shown
  await expect(page.getByRole("combobox", { name: "Saved deconvolution run" })).toHaveValue(/^Run \d+: nnls$/);
  const cellType = page.getByRole("combobox", { name: "Display cell type" });
  await expect(cellType).not.toHaveValue("");
  await cellType.click();
  const second = page.getByRole("option").nth(1);
  const name = (await second.textContent())!;
  const shown = page.waitForResponse((r) => r.url().includes("/style") && r.ok());
  await second.click();
  await shown;
  await expect(cellType).toHaveValue(name);

  await roiGroup(page, "Left", -150);
  await roiGroup(page, "Right", 150);
  await page.getByRole("tab", { name: "Metadata" }).click();
  await expect(page.getByRole("combobox", { name: "Metadata set" })).toHaveValue(/nnls fractions$/);
  await page.getByRole("checkbox", { name: "ROI: Left" }).check();
  await page.getByRole("checkbox", { name: "ROI: Right" }).check();
  await page.getByRole("button", { name: "Compare metadata" }).click();
  const results = page.getByTestId("metadata-results");
  await expect.poll(() => results.locator(".bk-Figure").count(), { timeout: 60_000 }).toBe(1);

  expect(errors.filter((e) => !/favicon/.test(e))).toEqual([]);
});

test("cancel a running deconvolution", async ({ page }) => {
  test.setTimeout(900_000);
  await load(page, DATASET);
  await page.getByRole("tab", { name: "Deconv." }).click();
  await expect(page.getByText(/^Reference ready: /)).toBeVisible({ timeout: 60_000 });
  const run = page.getByRole("button", { name: "Run deconvolution" });
  await expect(run).toBeEnabled({ timeout: 30_000 });
  // An unusual batch size is not cached in the dataset folder, so the run is still going when cancelled
  await page.getByRole("textbox", { name: "Observations per batch" }).fill("7");
  await expect(run).toBeEnabled({ timeout: 30_000 });
  await run.click();
  const cancel = page.getByRole("button", { name: "Cancel", exact: true });
  await expect(cancel).toBeVisible({ timeout: 60_000 });
  await page.waitForTimeout(3_000);                  // past the input export, into the plugin process
  await cancel.click();
  await expect(page.getByText("Deconvolution cancelled. No result was saved.")).toBeVisible({ timeout: 60_000 });
  await expect(cancel).toBeHidden();
  await expect(run).toBeEnabled();
});
