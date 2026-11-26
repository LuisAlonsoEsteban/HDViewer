// Error notices for failed background requests (a toast in the bottom-right
// corner), so failures are never silent and never leave a spinner forever.
import { Notification, Stack } from "@mantine/core";
import { create } from "zustand";
import { ApiError } from "./api/client";

interface Notice { id: number; title: string; message: string }

const useNotices = create<{ notices: Notice[]; close: (id: number) => void }>((set) => ({
  notices: [],
  close: (id) => set((s) => ({ notices: s.notices.filter((n) => n.id !== id) })),
}));

let next = 0;

/** Show `exc` (an ApiError, Error or text) as an error notice. Aborted requests
 * and capacity responses (which have their own notice) are ignored. */
export function notifyError(exc: unknown, title = "Something went wrong"): void {
  if (exc instanceof DOMException && exc.name === "AbortError") return;
  if (exc instanceof ApiError && exc.reason === "capacity") return;
  const message = exc instanceof Error ? exc.message : String(exc);
  console.error(exc);
  const id = ++next;
  useNotices.setState((s) => ({
    // The same failure repeated (e.g. while the server restarts) is shown once
    notices: [...s.notices.filter((n) => n.message !== message), { id, title, message }].slice(-3),
  }));
  setTimeout(() => useNotices.getState().close(id), 10000);
}

/** A rejection handler for promises whose failure only needs reporting. */
export const reportError = (title?: string) => (exc: unknown) => notifyError(exc, title);

export function ErrorNotices() {
  const { notices, close } = useNotices();
  if (!notices.length) return null;
  return (
    <Stack gap="xs" style={{ position: "fixed", right: 16, bottom: 16, zIndex: 1000, maxWidth: 380 }}>
      {notices.map((n) => (
        <Notification key={n.id} color="red" title={n.title} onClose={() => close(n.id)}>{n.message}</Notification>
      ))}
    </Stack>
  );
}
