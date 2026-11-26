import "@mantine/core/styles.css";
import "./app.css";
import { MantineProvider } from "@mantine/core";
import { createRoot } from "react-dom/client";
import { App } from "./App";
import { ErrorBoundary } from "./ErrorBoundary";
import { theme } from "./theme";

createRoot(document.getElementById("root")!).render(
  <MantineProvider theme={theme} defaultColorScheme="dark">
    <ErrorBoundary>
      <App />
    </ErrorBoundary>
  </MantineProvider>,
);
