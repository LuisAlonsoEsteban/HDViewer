import { useState } from "react";
import { Anchor, Button, Group, Paper, Text } from "@mantine/core";
import { useApp } from "../store/app";

const KEY = "hdv:cookie-notice";

/** Informs visitors about cookies. Only a strictly necessary session cookie is set,
 * so no consent is collected; the dismissal itself lives in sessionStorage, which is
 * cleared with the tab, so nothing on the visitor's device outlasts the visit. */
export function CookieNotice() {
  const [open, setOpen] = useState(() => {
    try { return sessionStorage.getItem(KEY) !== "seen"; } catch { return true; }
  });
  const privacy = useApp((s) => s.config?.privacy_url) ?? "static/assets/privacy.html";
  if (!open) return null;

  const dismiss = () => {
    try { sessionStorage.setItem(KEY, "seen"); } catch { /* storage blocked: just hide */ }
    setOpen(false);
  };

  return (
    <Paper role="region" aria-label="Cookie notice" shadow="md" withBorder radius="md" p="sm"
           style={{ position: "fixed", left: 16, right: 16, bottom: 16, zIndex: 1000,
                    maxWidth: 720, marginInline: "auto" }}>
      <Group justify="space-between" wrap="nowrap" gap="md" align="center">
        <Text size="sm">
          This site uses only a session cookie, needed to keep your analysis workspace while you use it;
          it is deleted when you close the browser. No permanent or third-party cookies are set by HDViewer.{" "}
          <Anchor href={privacy} target="_blank" rel="noopener" size="sm">Privacy and data retention</Anchor>
        </Text>
        <Button size="xs" variant="light" onClick={dismiss} style={{ flexShrink: 0 }}>Got it</Button>
      </Group>
    </Paper>
  );
}
