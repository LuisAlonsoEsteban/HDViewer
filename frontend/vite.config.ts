import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

// `npm run dev` proxies the API (uvicorn hd_viewer.api.main:app --port 8000);
// `npm run build` writes the bundle FastAPI serves (hd_viewer/api/dist).
export default defineConfig({
  plugins: [react()],
  base: "./",
  build: { outDir: "../hd_viewer/api/dist", emptyOutDir: true },
  server: {
    proxy: {
      "/api": { target: "http://127.0.0.1:8000", ws: true },
      "/static": "http://127.0.0.1:8000",
    },
  },
  test: { environment: "node", include: ["src/**/*.test.ts"] },
});
