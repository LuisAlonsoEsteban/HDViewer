// Body of a Bokeh CustomJS attached to an upload field's "Choose file…" button
// (see hd_viewer/startup/uploads.py). Args: channel, accept, chunk_bytes,
// load_button, source.
//
// Opens a native file picker, then streams the chosen file to the server in
// slices of chunk_bytes: each slice is base64-encoded and written to
// channel.request, and the next slice is only read once the server has
// acknowledged the previous one on channel.response. Choosing another file
// for the same field abandons the upload in progress.
const state = channel.__hd_upload ?? (channel.__hd_upload = {
  generation: 0,
  waiter: null,
  input: null,
  on_response(text) {
    let message;
    try { message = JSON.parse(text); } catch (error) { return; }
    const waiter = this.waiter;
    if (waiter == null || message.id !== waiter.id) return;
    if (message.index === waiter.index || !message.ok) {
      this.waiter = null;
      waiter.resolve(message);
    }
  },
});

function wait_for(id, index) {
  return new Promise((resolve) => {
    state.waiter?.resolve({ok: false, error: "superseded"});
    state.waiter = {id, index, resolve};
  });
}

function read_base64(blob) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => {
      const url = reader.result;
      resolve(url.slice(url.indexOf(",") + 1));
    };
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(blob);
  });
}

function send(message) {
  channel.setv({request: JSON.stringify(message)});
}

async function upload(file) {
  const generation = ++state.generation;
  const id = Date.now().toString(36) + "-" + Math.random().toString(36).slice(2);
  if (source.active === 0) load_button.disabled = true;
  // The server validates the name, extension and size limit before any data is sent.
  let reply = wait_for(id, -1);
  send({action: "start", id, name: file.name, size: file.size});
  if (!(await reply).ok || generation !== state.generation) return;
  for (let index = 0, offset = 0; offset < file.size; index++, offset += chunk_bytes) {
    let data;
    try {
      data = await read_base64(file.slice(offset, offset + chunk_bytes));
    } catch (error) {
      // Tell the server so it reports the failure instead of waiting forever.
      data = "!";
    }
    if (generation !== state.generation) return;
    reply = wait_for(id, index);
    send({action: "chunk", id, index, data});
    data = null;
    if (!(await reply).ok || generation !== state.generation) return;
  }
}

state.input?.remove();
const input = document.createElement("input");
input.type = "file";
input.accept = accept;
input.style.display = "none";
state.input = input;
input.addEventListener("change", () => {
  const file = input.files?.[0];
  input.remove();
  if (state.input === input) state.input = null;
  if (file != null) upload(file);
});
input.addEventListener("cancel", () => input.remove());
document.body.appendChild(input);
input.click();
