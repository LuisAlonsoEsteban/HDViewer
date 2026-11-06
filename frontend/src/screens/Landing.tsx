import { Anchor, Button, Container, Group, Stack, Text, Title } from "@mantine/core";
import { useApp } from "../store/app";

export function Landing() {
  const { config, setScreen } = useApp();
  const branding = config?.branding ?? {};
  return (
    <Container size="md" py={80}>
      <Stack gap="lg">
        <Text c="teal.3" fw={600} tt="uppercase" size="sm">Visium HD · spatial transcriptomics</Text>
        <Title order={1}>{branding.name ?? "VisiumVisualizer"}</Title>
        <Text size="lg" c="dark.1">
          Explore Visium HD expression, cell populations and tissue neighborhoods on top of the original
          tissue image: 2–16 µm bins, Visium-like pseudospots and segmented cells, in your browser.
          Free and open to everyone, no login required. Uploaded data stays private to your session and
          is deleted when it ends.
        </Text>
        <Group>
          <Button size="md" onClick={() => setScreen("loading")}>Launch</Button>
          {config && <Button size="md" variant="default" component="a" href={config.tutorial_url} target="_blank">Tutorial</Button>}
          {config && <Button size="md" variant="subtle" component="a" href={config.use_case_url} target="_blank">Biological use case</Button>}
        </Group>
        <Text size="sm" c="dark.2">
          {branding.citation}{" "}
          {branding.repository && <Anchor href={branding.repository} target="_blank">Source code</Anchor>}
          {" · "}
          {branding.issues && <Anchor href={branding.issues} target="_blank">Report a problem</Anchor>}
          {branding.email && <> · <Anchor href={`mailto:${branding.email}`}>{branding.email}</Anchor></>}
        </Text>
      </Stack>
    </Container>
  );
}
