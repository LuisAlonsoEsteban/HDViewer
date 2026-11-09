import { useEffect, useState } from "react";
import { ActionIcon, Anchor, Group, Popover, Stack, Text, Title } from "@mantine/core";
import { getJson } from "../api/client";

interface Topic {
  title: string;
  html: string;           // server-owned help text (hd_viewer/help/help.py)
  tutorial: string | null;
}

let topics: Promise<Record<string, Topic>> | null = null;

function loadTopics() {
  topics ??= getJson<Record<string, Topic>>("help").catch(() => {
    topics = null;
    return {};
  });
  return topics;
}

/** "?" button explaining `topic` (or several `topics` in one popover). */
export function HelpButton({ topic, topics }: { topic?: string; topics?: string[] }) {
  const keys = topics ?? (topic ? [topic] : []);
  const [entries, setEntries] = useState<Topic[]>([]);
  useEffect(() => {
    loadTopics().then((all) => setEntries(keys.map((key) => all[key]).filter(Boolean)));
  }, [keys.join("|")]);
  if (!entries.length) return null;
  return (
    <Popover width={360} position="right-start" withArrow shadow="md">
      <Popover.Target>
        <ActionIcon size="sm" radius="xl" variant="default" aria-label={`Help: ${entries.map((e) => e.title).join(" · ")}`}>?</ActionIcon>
      </Popover.Target>
      <Popover.Dropdown>
        <Stack gap="sm">
          {entries.map((entry) => (
            <div key={entry.title}>
              <Title order={6} mb={4}>{entry.title}</Title>
              <Text size="xs" component="div" className="help-body" dangerouslySetInnerHTML={{ __html: entry.html }} />
              {entry.tutorial && (
                <Anchor size="xs" href={entry.tutorial} target="_blank" rel="noopener" mt={6} display="block">
                  Read more in the tutorial →
                </Anchor>
              )}
            </div>
          ))}
        </Stack>
      </Popover.Dropdown>
    </Popover>
  );
}

/** A section title with its help button. */
export function HelpTitle({ title, topic }: { title: string; topic: string }) {
  return (
    <Group gap={6} wrap="nowrap">
      <Text size="sm" fw={600}>{title}</Text>
      <HelpButton topic={topic} />
    </Group>
  );
}
