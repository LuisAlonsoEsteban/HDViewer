import { expect, test } from "@playwright/test";

// Needs a server whose session limit is 1 (VHD_MAX_SESSIONS=1)
test("a full server shows a friendly notice that retries and offers running locally", async ({ browser }) => {
  test.setTimeout(300_000);
  const first = await browser.newContext();
  const firstPage = await first.newPage();
  await firstPage.goto("/");
  await expect(firstPage.getByRole("heading", { name: /See the tissue/ })).toBeVisible();
  await expect.poll(() => first.cookies().then((c) => c.some((x) => x.name === "hdv_sid"))).toBe(true);

  const second = await browser.newContext();
  const page = await second.newPage();
  const errors: string[] = [];
  page.on("pageerror", (error) => errors.push(String(error)));
  await page.goto("/");
  const notice = page.getByRole("dialog", { name: "The server is full right now" });
  await expect(notice).toBeVisible();
  await expect(notice).toContainText("Nothing went wrong");
  await expect(notice.getByRole("status")).toContainText(/Next check in \d+ s/);
  await expect(notice.getByRole("link", { name: "run it on your own computer" })).toHaveAttribute("href", /#run-it-locally$/);
  await page.screenshot({ path: "test-results/capacity-notice.png" });

  expect(errors).toEqual([]);
});
