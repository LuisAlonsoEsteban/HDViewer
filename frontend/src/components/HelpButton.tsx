import { useEffect, useState } from "react";
import { ActionIcon, Anchor, Group, Popover, Text, Title } from "@mantine/core";
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

/** "?" button explaining `topic` (the old interface's help popovers). */
export function HelpButton({ topic }: { topic: string }) {
  const [entry, setEntry] = useState<Topic | null>(null);
  useEffect(() => {
    loadTopics().then((all) => setEntry(all[topic] ?? null));
  }, [topic]);
  if (!entry) return null;
  return (
    <Popover width={340} position="right-start" withArrow shadow="md">
      <Popover.Target>
        <ActionIcon size="xs" radius="xl" variant="default" aria-label={`Help: ${entry.title}`}>?</ActionIcon>
      </Popover.Target>
      <Popover.Dropdown>
        <Title order={6} mb={4}>{entry.title}</Title>
        <Text size="xs" component="div" className="help-body" dangerouslySetInnerHTML={{ __html: entry.html }} />
        {entry.tutorial && (
          <Anchor size="xs" href={entry.tutorial} target="_blank" rel="noopener" mt={6} display="block">
            Read more in the tutorial →
          </Anchor>
        )}
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
