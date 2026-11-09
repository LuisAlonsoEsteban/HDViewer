import { test } from "@playwright/test";
import { load } from "../helpers";

// Not a test: captures the workspace for design review (run on demand)
test("capture workspace screens", async ({ page }) => {
  test.setTimeout(600_000);
  await page.setViewportSize({ width: 1600, height: 950 });
  await load(page, process.env.HDV_DATASET ?? "Mouse brain (16um)", { pseudospots: true });
  await page.waitForTimeout(4000);
  const out = process.env.SHOTS ?? "test-results/shots";
  await page.screenshot({ path: `${out}/workspace.png` });
  for (const tab of ["Clustering", "ROIs", "Deconv."]) {
    await page.getByRole("tab", { name: tab, exact: true }).first().click();
    await page.waitForTimeout(800);
    await page.screenshot({ path: `${out}/tab-${tab.replace(".", "")}.png` });
  }
});
