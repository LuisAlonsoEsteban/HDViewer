import { expect, test } from "@playwright/test";
import { load } from "./helpers";

// The P1CRC guided tour, up to the RCTD step (RCTD itself is not run here)
test("guided tour: load P1CRC and reach the deconvolution step", async ({ page, request }) => {
  test.setTimeout(900_000);
  const datasets: { id: string }[] = await (await request.get("api/datasets")).json();
  test.skip(!datasets.some((d) => d.id === "Human CRC"), "Human CRC example not on this server");

  await page.goto("/?tour=p1crc");
  const card = page.getByRole("dialog", { name: /Guided tour/ });
  await expect(card).toContainText("Which cells surround a colorectal tumor?");
  expect(page.url()).not.toContain("tour=");
  await card.getByRole("button", { name: "Next" }).click();

  await expect(card).toContainText("Set up a dataset");
  await expect(page.locator(".tour-spotlight")).toBeVisible();
  await page.getByRole("button", { name: /Explore your data/ }).click();

  // Done by the user, so the tour moves on by itself, with the example preselected
  await expect(card).toContainText("Load the P1CRC sample");
  await expect(page.getByRole("combobox", { name: "Example dataset" })).toHaveValue("Human CRC");
  await expect(page.getByLabel("Generate Visium v2 pseudospots (55 µm)")).not.toBeChecked();
  await page.getByRole("button", { name: "Load data" }).click();

  await expect(card).toContainText("Check the bins over the tissue", { timeout: 600_000 });
  await expect(card).toContainText("✓ Done");          // the bins layer is already active
  await card.getByRole("button", { name: "Next" }).click();

  // The step opens its tool and outlines its controls
  await expect(card).toContainText("Infer cell types with RCTD");
  await expect(page.getByRole("tab", { name: "Deconv." })).toHaveAttribute("aria-selected", "true");
  await expect(page.locator(".tour-spotlight")).toBeVisible();

  // A reload resumes the tour where it was
  await page.reload();
  await expect(card).toContainText("Infer cell types with RCTD");
  await card.getByRole("button", { name: "Exit the guided tour" }).click();
  await expect(card).toHaveCount(0);
});

// Tour mechanics in the workspace, on any example (it warns that the data is not Human CRC)
test("guided tour: workspace steps open their tool, outline it and survive a reload", async ({ page }) => {
  test.setTimeout(900_000);
  await load(page, process.env.HDV_DATASET ?? "Mouse brain (16um)");
  await page.getByRole("button", { name: "Guided tour" }).click();
  const card = page.getByRole("dialog", { name: /Guided tour/ });
  for (let i = 0; i < 10 && !(await card.textContent())?.includes("Infer cell types with RCTD"); i++) {
    await card.getByRole("button", { name: /Next|Skip/ }).click();
  }
  await expect(card).toContainText("Infer cell types with RCTD");
  await expect(card).toContainText("The loaded data is not Human CRC");
  await expect(page.getByRole("tab", { name: "Deconv." })).toHaveAttribute("aria-selected", "true");
  await expect(page.locator(".tour-spotlight")).toBeVisible();

  await card.getByRole("button", { name: /Next|Skip/ }).click();
  await card.getByRole("button", { name: /Next|Skip/ }).click();
  await expect(card).toContainText("Define the tumor");
  await expect(page.getByRole("tab", { name: "Groups" })).toHaveAttribute("aria-selected", "true");

  await page.reload();
  await expect(card).toContainText("Define the tumor");
  // The app reopens the loaded workspace, and the step its tool
  await expect(page.getByRole("tab", { name: "Groups" })).toHaveAttribute("aria-selected", "true");
  await card.getByRole("button", { name: "Exit the guided tour" }).click();
  await expect(card).toHaveCount(0);
});
