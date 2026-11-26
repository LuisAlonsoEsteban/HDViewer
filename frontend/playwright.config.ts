import { defineConfig, devices } from "@playwright/test";

// Runs against a live server: uvicorn hd_viewer.api.main:app (HDV_URL, default :8765).
// Chromium by default; HDV_BROWSERS=all also runs Firefox and WebKit (Safari's engine),
// after `npx playwright install firefox webkit`.
const viewport = { width: 1400, height: 900 };
const browsers = process.env.HDV_BROWSERS === "all"
  ? [["chromium", devices["Desktop Chrome"]], ["firefox", devices["Desktop Firefox"]], ["webkit", devices["Desktop Safari"]]] as const
  : [["chromium", devices["Desktop Chrome"]]] as const;

export default defineConfig({
  testDir: "e2e",
  // e2e/tools: screenshot scripts, not tests; run them with HDV_TOOLS=1
  testIgnore: process.env.HDV_TOOLS ? [] : ["tools/**"],
  timeout: 240_000,
  use: { baseURL: process.env.HDV_URL ?? "http://127.0.0.1:8765", viewport },
  projects: browsers.map(([name, device]) => ({ name, use: { ...device, viewport } })),
});
