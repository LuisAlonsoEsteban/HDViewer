import { create } from "zustand";

/** The running guided tour (if any) and its step. Mirrored to sessionStorage, so a
 * reload resumes the tour; nothing outlives the browser session. */
const KEY = "hdv:tour";

interface Saved { tour: string | null; step: number }

function read(): Saved {
  try {
    const saved = JSON.parse(sessionStorage.getItem(KEY) ?? "null");
    if (saved && typeof saved.tour === "string" && Number.isInteger(saved.step)) return saved;
  } catch { /* storage blocked or malformed: start without a tour */ }
  return { tour: null, step: 0 };
}

function write(saved: Saved) {
  try {
    if (saved.tour) sessionStorage.setItem(KEY, JSON.stringify(saved));
    else sessionStorage.removeItem(KEY);
  } catch { /* storage blocked: the tour just won't survive a reload */ }
}

interface TourState extends Saved {
  start(tour: string): void;
  go(step: number): void;
  next(): void;
  back(): void;
  stop(): void;
}

export const useTour = create<TourState>((set, get) => {
  const update = (saved: Saved) => { write(saved); set(saved); };
  return {
    ...(typeof window === "undefined" ? { tour: null, step: 0 } : read()),
    start: (tour) => update({ tour, step: 0 }),
    go: (step) => update({ tour: get().tour, step: Math.max(0, step) }),
    next: () => update({ tour: get().tour, step: get().step + 1 }),
    back: () => update({ tour: get().tour, step: Math.max(0, get().step - 1) }),
    stop: () => update({ tour: null, step: 0 }),
  };
});
