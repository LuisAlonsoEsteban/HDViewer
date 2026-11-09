import { createTheme } from "@mantine/core";

// Workspace colours (dark greens, shared with the landing page's palette)
export const theme = createTheme({
  primaryColor: "teal",
  fontFamily: '"Aptos", "Segoe UI", system-ui, sans-serif',
  defaultRadius: "md",
  colors: {
    dark: ["#e8f0e9", "#adc5b7", "#779887", "#3c584c", "#2a4738", "#253e32", "#192d24", "#14251f", "#11231d", "#0b1712"],
  },
});
