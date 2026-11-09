import { useEffect, useState } from "react";
import { Group, Loader, Stack, Text } from "@mantine/core";
import { Job } from "../store/app";

/** While a job runs: spinner, elapsed time and its recent steps (the old
 * interface's loading log). Nothing once it has finished. */
export function JobProgress({ job, message = true }: { job?: Job; message?: boolean }) {
  const [now, setNow] = useState(Date.now());
  const running = job?.state === "running";
  useEffect(() => {
    if (!running) return;
    const timer = setInterval(() => setNow(Date.now()), 1000);
    return () => clearInterval(timer);
  }, [running]);
  if (!running || !job) return null;
  const seconds = Math.max(0, Math.round((now - (job.started ?? now)) / 1000));
  const earlier = (job.log ?? []).slice(0, -1);
  return (
    <Stack gap={2} data-testid="job-progress">
      <Group gap={6}>
        <Loader size="xs" />
        <Text size="xs" c="dark.1">{message ? `${job.message || "Working…"} · ` : "Running · "}{seconds} s</Text>
      </Group>
      {earlier.map((line, i) => <Text key={i} size="xs" c="dimmed" pl={22}>{line}</Text>)}
    </Stack>
  );
}
