// Thin client for the VisiumVisualizer API. Paths are relative, so the app
// works under any prefix (e.g. /v2/ behind nginx).
import { decodeFrame, Frame } from "./frames";

export class ApiError extends Error {
  constructor(public status: number, message: string, public reason?: string) {
    super(message);
  }
}

/** Called when the server has no free analysis session (an expected, temporary
 * state): the app shows a notice instead of an error. Seconds until a retry. */
let capacityHandler: ((retryAfter: number) => void) | null = null;
export function onCapacity(handler: (retryAfter: number) => void): void {
  capacityHandler = handler;
}

async function check(response: Response): Promise<Response> {
  if (!response.ok) {
    let message = response.statusText;
    let reason: string | undefined;
    let retryAfter = Number(response.headers.get("Retry-After")) || 60;
    try {
      const body = await response.json();
      message = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail);
      reason = body.reason;
      if (typeof body.retry_after === "number") retryAfter = body.retry_after;
    } catch {
      /* not JSON */
    }
    if (response.status === 503 && reason === "capacity") capacityHandler?.(retryAfter);
    throw new ApiError(response.status, message, reason);
  }
  return response;
}

export async function getJson<T>(path: string): Promise<T> {
  return (await check(await fetch(`api/${path}`, { credentials: "same-origin" }))).json();
}

export async function sendJson<T>(method: string, path: string, body?: unknown, signal?: AbortSignal): Promise<T> {
  const response = await fetch(`api/${path}`, {
    method,
    signal,
    credentials: "same-origin",
    headers: { "Content-Type": "application/json" },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  return (await check(response)).json();
}

export async function frame(method: string, path: string, body?: unknown): Promise<Frame> {
  const response = await fetch(`api/${path}`, {
    method,
    credentials: "same-origin",
    headers: body === undefined ? undefined : { "Content-Type": "application/json" },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  return decodeFrame(await (await check(response)).arrayBuffer());
}

/** Upload `file` into `field` in raw-byte chunks; `onProgress(sentBytes)`. */
export async function uploadFile(field: string, file: File, onProgress: (sent: number) => void,
                                 signal?: AbortSignal): Promise<void> {
  const { upload_id, chunk_bytes } = await sendJson<{ upload_id: string; chunk_bytes: number }>(
    "POST", "uploads", { field, filename: file.name, size: file.size });
  for (let index = 0, offset = 0; offset < file.size; index++, offset += chunk_bytes) {
    if (signal?.aborted) throw new DOMException("Upload cancelled", "AbortError");
    const chunk = file.slice(offset, offset + chunk_bytes);
    await check(await fetch(`api/uploads/${field}/${upload_id}/${index}`, {
      method: "PUT", credentials: "same-origin", body: chunk, signal,
    }));
    onProgress(Math.min(offset + chunk_bytes, file.size));
  }
}

export type ServerEvent =
  | { type: "status"; load: LoadStatus }
  | { type: "progress"; task: string; pct: number; status: string }
  | { type: "failed"; task: string | null; message: string; fatal: boolean }
  | { type: "ready" }
  | { type: "job"; layer: string; name: string; state: "running" | "done" | "failed"; message: string }
  | { type: "layer_stale"; layer: string }
  | { type: "ping" };

export interface LoadStatus {
  state?: "running" | "ready" | "failed";
  tasks?: Record<string, { pct: number; status: string }>;
  error?: string | null;
}

/** Events socket with automatic reconnection (the session survives reloads). */
export function connectEvents(onEvent: (event: ServerEvent) => void): () => void {
  let socket: WebSocket | null = null;
  let closed = false;
  let retry = 0;
  const open = () => {
    const url = new URL("api/events", window.location.href);
    url.protocol = url.protocol === "https:" ? "wss:" : "ws:";
    socket = new WebSocket(url);
    socket.onmessage = (message) => onEvent(JSON.parse(message.data));
    socket.onopen = () => (retry = 0);
    socket.onclose = () => {
      if (!closed) setTimeout(open, Math.min(10000, 500 * 2 ** retry++));
    };
  };
  open();
  return () => {
    closed = true;
    socket?.close();
  };
}
