import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Relative asset paths: the built page is served from any static folder.
export default defineConfig({ plugins: [react()], base: "./", build: { outDir: "dist", chunkSizeWarningLimit: 100000 } });
