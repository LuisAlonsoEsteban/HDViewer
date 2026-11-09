import { expect, test } from "@playwright/test";

const DATASET = process.env.HDV_DATASET ?? "Mouse brain (16um)";

test("cluster the layer, hide a cluster and create a group", async ({ page }) => {
  test.setTimeout(480_000);
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

  // Clustering
  await page.getByRole("tab", { name: "Clustering" }).click();
  await page.getByRole("button", { name: "Run clustering" }).click();
  await expect(page.getByText(/clusters found, marker genes ready/)).toBeVisible({ timeout: 300_000 });
  const first = page.getByRole("checkbox", { name: /^Cluster 0/ });
  await expect(first).toBeChecked();
  const hidden = page.waitForResponse((r) => r.url().includes("/clusters") && r.request().method() === "PATCH" && r.ok());
  await first.uncheck();
  await hidden;
  await expect.poll(() => page.evaluate(() => {
    const doc = (window as any).hdv.docs.get("map");
    const booleans = doc.get_model_by_name("bins:filter").booleans;
    return booleans ? booleans.filter((b: boolean) => !b).length : 0;
  })).toBeGreaterThan(0);

  // Groups
  await page.getByRole("tab", { name: "Groups" }).click();
  await page.getByRole("combobox", { name: "Gene / metadata" }).fill("Ttr");
  await page.getByRole("option", { name: "Ttr", exact: true }).click();
  await page.getByRole("textbox", { name: "Value" }).fill("1");
  await page.getByRole("button", { name: "Add condition" }).click();
  await expect(page.getByText("Ttr > 1")).toBeVisible();
  await page.getByRole("textbox", { name: "Group name" }).fill("Choroid plexus");
  await page.getByRole("button", { name: "Create group" }).click();
  await expect(page.getByText("Group 'Choroid plexus' created.")).toBeVisible({ timeout: 120_000 });
  await expect(page.getByRole("checkbox", { name: /Choroid plexus/ })).toBeChecked();

  // A second group, then combined groups: the difference and a complement
  await page.getByRole("textbox", { name: "Value" }).fill("5");
  await page.getByRole("button", { name: "Add condition" }).click();
  await page.getByRole("textbox", { name: "Group name" }).fill("High Ttr");
  await page.getByRole("button", { name: "Create group" }).click();
  await expect(page.getByText("Group 'High Ttr' created.")).toBeVisible({ timeout: 120_000 });

  const combine = page.getByTestId("combine-groups");
  await expect(combine.getByRole("combobox", { name: "A", exact: true })).toHaveValue("Choroid plexus");
  await expect(combine.getByRole("combobox", { name: "B", exact: true })).toHaveValue("High Ttr");
  await combine.getByText("A − B", { exact: true }).click();
  await expect(combine.getByRole("status")).toContainText(/Choroid plexus − High Ttr: [\d,]+ of [\d,]+ observations/);
  await combine.getByRole("button", { name: "Create combined group" }).click();
  await expect(combine.getByText(/Group 'Choroid plexus − High Ttr' created/)).toBeVisible({ timeout: 60_000 });
  await expect(page.getByRole("checkbox", { name: /Choroid plexus − High Ttr/ })).toBeChecked();

  await combine.getByText("not A", { exact: true }).click();
  await expect(combine.getByRole("combobox", { name: "B", exact: true })).toHaveCount(0);
  await combine.getByRole("textbox", { name: "Name (optional)" }).fill("Outside the plexus");
  await combine.getByRole("button", { name: "Create combined group" }).click();
  await expect(page.getByRole("checkbox", { name: /Outside the plexus/ })).toBeChecked({ timeout: 60_000 });
  await expect(page.locator(".mantine-Checkbox-description", { hasText: "not Choroid plexus" })).toBeVisible();   // how it was made

  expect(errors.filter((e) => !/favicon/.test(e))).toEqual([]);
});
