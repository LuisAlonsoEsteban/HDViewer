import { useEffect, useState } from "react";
import { Anchor, Button, Group, Modal, Progress, Stack, Text, Title } from "@mantine/core";
import { useApp } from "../store/app";

/** Shown when every analysis session on the server is in use: a normal, temporary
 * state of a shared public server, not an error. Retries by itself. */
export function CapacityNotice({ retryAfter, onRetry }: { retryAfter: number; onRetry(): Promise<boolean> }) {
  const repository = useApp((s) => s.config?.branding.repository);
  const [left, setLeft] = useState(retryAfter);
  const [checking, setChecking] = useState(false);

  const retry = async () => {
    setChecking(true);
    const free = await onRetry();
    setChecking(false);
    if (!free) setLeft(retryAfter);
  };

  useEffect(() => {
    if (checking) return;
    if (left <= 0) {
      retry();
      return;
    }
    const timer = setTimeout(() => setLeft((s) => s - 1), 1000);
    return () => clearTimeout(timer);
  }, [left, checking]);

  return (
    <Modal opened onClose={() => {}} withCloseButton={false} centered size="lg"
           title={<Title order={3}>The server is full right now</Title>}>
      <Stack gap="sm">
        <Text>
          Nothing went wrong: this public server keeps a limited number of analysis sessions open at the same
          time, so that every visitor gets a responsive workspace, and right now they are all in use.
        </Text>
        <Text>
          A place frees up as soon as another visitor finishes or leaves. This page checks again automatically;
          you can keep it open.
        </Text>
        <div>
          <Group justify="space-between" mb={4}>
            <Text size="sm" c="dimmed" role="status" aria-live="polite">
              {checking ? "Checking for a free session…" : `Next check in ${left} s`}
            </Text>
            <Button size="xs" variant="light" onClick={retry} loading={checking}>Try now</Button>
          </Group>
          <Progress value={checking ? 100 : ((retryAfter - left) / retryAfter) * 100} animated={checking} />
        </div>
        <Text size="sm" c="dimmed">
          Have large datasets, or prefer not to wait? VisiumVisualizer is free and open source: you can{" "}
          {repository
            ? <Anchor href={`${repository}#run-it-locally`} target="_blank" rel="noopener" size="sm">run it on your own computer</Anchor>
            : "run it on your own computer"}
          , with the same features and no session limit.
        </Text>
      </Stack>
    </Modal>
  );
}
