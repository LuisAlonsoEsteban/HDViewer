import { Component, ReactNode } from "react";
import { Alert, Button, Container, Stack, Text } from "@mantine/core";

/** Last resort for a rendering error: a message and a reload button instead of a
 * blank page. The session (and a loaded dataset) survives the reload. */
export class ErrorBoundary extends Component<{ children: ReactNode }, { error: Error | null }> {
  state = { error: null as Error | null };

  static getDerivedStateFromError(error: Error) {
    return { error };
  }

  componentDidCatch(error: Error) {
    console.error(error);
  }

  render() {
    if (!this.state.error) return this.props.children;
    return (
      <Container size="sm" py="xl">
        <Alert color="red" title="The page stopped working">
          <Stack gap="sm">
            <Text size="sm">
              An unexpected error occurred while drawing the page. Reloading keeps your session and loaded data.
            </Text>
            <Text size="xs" c="dimmed" style={{ fontFamily: "monospace" }}>{this.state.error.message}</Text>
            <Button w="fit-content" onClick={() => window.location.reload()}>Reload</Button>
          </Stack>
        </Alert>
      </Container>
    );
  }
}
