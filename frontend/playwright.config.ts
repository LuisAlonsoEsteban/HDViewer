import { defineConfig } from "@playwright/test";

// Runs against a live server: uvicorn hd_viewer.api.main:app (HDV_URL, default :8765)
export default defineConfig({
  testDir: "e2e",
  timeout: 240_000,
  use: { baseURL: process.env.HDV_URL ?? "http://127.0.0.1:8765", viewport: { width: 1400, height: 900 } },
});
