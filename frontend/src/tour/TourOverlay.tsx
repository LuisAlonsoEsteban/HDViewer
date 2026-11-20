import { useEffect, useLayoutEffect, useRef, useState } from "react";
import { ActionIcon, Anchor, Badge, Button, Group, Paper, Portal, Progress, Text } from "@mantine/core";
import { useApp } from "../store/app";
import { useTourFacts } from "./facts";
import { TOURS } from "./p1crc";
import { useTour } from "./store";
import type { Facts, Tour } from "./types";

const CARD_WIDTH = 360;
const GAP = 14;
const MARGIN = 16;
const PAD = 6;                       // spotlight padding around the anchors
const ADVANCE_MS = 1200;

interface Box { left: number; top: number; right: number; bottom: number }

/** The union of the visible elements tagged with any of the keys. */
function anchorBox(keys: string[] | undefined): Box | null {
  if (!keys?.length) return null;
  let box: Box | null = null;
  for (const key of keys) {
    document.querySelectorAll(`[data-tour="${key}"]`).forEach((element) => {
      const r = element.getBoundingClientRect();
      if (r.width === 0 && r.height === 0) return;
      box = box
        ? { left: Math.min(box.left, r.left), top: Math.min(box.top, r.top),
            right: Math.max(box.right, r.right), bottom: Math.max(box.bottom, r.bottom) }
        : { left: r.left, top: r.top, right: r.right, bottom: r.bottom };
    });
  }
  return box;
}

const sameBox = (a: Box | null, b: Box | null) =>
  a === b || (!!a && !!b && a.left === b.left && a.top === b.top && a.right === b.right && a.bottom === b.bottom);

/** Follows the anchors as the layout changes (tabs mounting, panels resizing, scrolling). */
function useAnchorBox(keys: string[] | undefined): Box | null {
  const [box, setBox] = useState<Box | null>(null);
  const key = keys?.join(" ");
  useEffect(() => {
    let frame = 0;
    let last: Box | null = null;
    const tick = () => {
      const next = anchorBox(keys);
      if (!sameBox(next, last)) setBox((last = next));
      frame = requestAnimationFrame(tick);
    };
    tick();
    return () => cancelAnimationFrame(frame);
  }, [key]);
  return box;
}

/** Beside the anchors if there is room (right, left, below, above), else docked bottom right. */
function place(box: Box | null, width: number, height: number) {
  const vw = window.innerWidth;
  const vh = window.innerHeight;
  const clampTop = (top: number) => Math.min(Math.max(MARGIN, top), vh - height - MARGIN);
  const clampLeft = (left: number) => Math.min(Math.max(MARGIN, left), vw - width - MARGIN);
  if (!box) return { left: (vw - width) / 2, top: Math.max(MARGIN, (vh - height) / 2) };
  if (vw >= 700) {
    if (box.right + GAP + width + MARGIN <= vw) return { left: box.right + GAP, top: clampTop(box.top) };
    if (box.left - GAP - width >= MARGIN) return { left: box.left - GAP - width, top: clampTop(box.top) };
    if (box.bottom + GAP + height + MARGIN <= vh) return { left: clampLeft(box.left), top: box.bottom + GAP };
    if (box.top - GAP - height >= MARGIN) return { left: clampLeft(box.left), top: box.top - GAP - height };
  }
  return { left: vw - width - MARGIN, top: vh - height - MARGIN };
}

export function TourOverlay() {
  const { tour: tourId, step: index, next, back, stop, go } = useTour();
  const tour: Tour | undefined = tourId ? TOURS[tourId] : undefined;
  const facts = useTourFacts(Boolean(tour), index);
  if (!tour) return null;
  const step = tour.steps[Math.min(index, tour.steps.length - 1)];
  return <TourCard key={`${tour.id}:${step.id}`} tour={tour} index={Math.min(index, tour.steps.length - 1)}
                   facts={facts} next={next} back={back} stop={stop} go={go} />;
}

function TourCard({ tour, index, facts, next, back, stop, go }: {
  tour: Tour; index: number; facts: Facts; next(): void; back(): void; stop(): void; go(step: number): void;
}) {
  const step = tour.steps[index];
  const config = useApp((s) => s.config);
  const setScreen = useApp((s) => s.setScreen);
  const setTool = useApp((s) => s.setTool);
  const last = index === tour.steps.length - 1;
  const done = step.done?.(facts) ?? false;
  // Auto-advance only when the user completes the step here, not when returning to a finished one
  const doneOnEntry = useRef(done);
  useEffect(() => {
    if (!done || doneOnEntry.current || last) return;
    const timer = setTimeout(next, ADVANCE_MS);
    return () => clearTimeout(timer);
  }, [done]);

  // Open the step's tool once the workspace is shown
  useEffect(() => {
    if (step.tool && facts.screen === "workspace") setTool(step.tool);
  }, [step.tool, facts.screen === "workspace"]);

  const elsewhere = step.screen && step.screen !== facts.screen;
  const wrongData = facts.screen === "workspace" && facts.hasWorkspace && facts.dataset !== undefined
    && facts.dataset !== tour.dataset && index > loadStep(tour);
  const box = useAnchorBox(elsewhere ? undefined : step.anchor);

  const card = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState({ width: CARD_WIDTH, height: 220 });
  useLayoutEffect(() => {
    const element = card.current;
    if (!element) return;
    const observer = new ResizeObserver(() => setSize({ width: element.offsetWidth, height: element.offsetHeight }));
    observer.observe(element);
    return () => observer.disconnect();
  }, []);
  const position = place(box, size.width, size.height);

  return (
    <Portal>
      {box && (
        <div className="tour-spotlight" aria-hidden
             style={{ left: box.left - PAD, top: box.top - PAD,
                      width: box.right - box.left + 2 * PAD, height: box.bottom - box.top + 2 * PAD }} />
      )}
      <Paper ref={card} className="tour-card" shadow="xl" withBorder radius="md" p="md" role="dialog"
             aria-label={`Guided tour: ${tour.title}`} aria-live="polite"
             style={{ left: position.left, top: position.top, width: `min(${CARD_WIDTH}px, calc(100vw - ${2 * MARGIN}px))` }}>
        <Group justify="space-between" wrap="nowrap" mb={4}>
          <Text size="xs" c="dimmed" tt="uppercase" fw={600} style={{ letterSpacing: ".08em" }}>
            Guided tour · {index + 1} of {tour.steps.length}
          </Text>
          <ActionIcon variant="subtle" color="gray" size="sm" onClick={stop} aria-label="Exit the guided tour">✕</ActionIcon>
        </Group>
        <Progress value={((index + 1) / tour.steps.length) * 100} size="xs" mb="sm" />
        <Group gap="xs" wrap="nowrap" mb={6} align="center">
          <Text fw={700}>{step.title}</Text>
          {done && <Badge color="teal" variant="light" size="sm">✓ Done</Badge>}
        </Group>
        {step.body(facts)}

        {wrongData && (
          <Text size="sm" c="yellow.4" mt="xs">
            The loaded data is not {tour.dataset}.{" "}
            <Anchor size="sm" component="button" onClick={() => { setScreen("loading"); go(loadStep(tour)); }}>
              Load {tour.dataset}
            </Anchor>
          </Text>
        )}
        {elsewhere && step.screen === "workspace" && facts.hasWorkspace && (
          <Button size="xs" variant="light" mt="xs" onClick={() => setScreen("workspace")}>Back to the workspace</Button>
        )}
        {elsewhere && step.screen === "loading" && (
          <Button size="xs" variant="light" mt="xs" onClick={() => setScreen("loading")}>Open the loading screen</Button>
        )}

        <Group justify="space-between" mt="md" wrap="nowrap">
          {step.readMore && config
            ? <Anchor size="xs" href={`${config.use_case_url}${step.readMore}`} target="_blank" rel="noopener">Details ↗</Anchor>
            : <span />}
          <Group gap="xs" wrap="nowrap">
            {index > 0 && <Button size="xs" variant="default" onClick={back}>Back</Button>}
            {last
              ? <Button size="xs" onClick={stop}>Finish</Button>
              : <Button size="xs" variant={done || !step.done ? "filled" : "light"} onClick={next}>
                  {done || !step.done ? "Next" : "Skip"}
                </Button>}
          </Group>
        </Group>
      </Paper>
    </Portal>
  );
}

const loadStep = (tour: Tour) => Math.max(0, tour.steps.findIndex((s) => s.id === "load"));
