
from ..settings.config import grid_fit_kwargs, image_kwargs, bins_to_spots_kwargs, segmentation_kwargs

from .. import settings
from bokeh.models import (
    FileInput, Button, Div, Select,
    Checkbox, CustomJS, RadioButtonGroup, TextInput, NumericInput,
)
from bokeh.layouts import row, column
from bokeh.events import DocumentReady

import base64
import hashlib
import os
import sys
import shutil
import uuid
from pathlib import Path
from html import escape

import threading
import traceback

from .data.VisiumHD import VisiumHD
from .data.image.segmentation import SegmentationService
from .datasets import AVAILABLE_DATASETS

from ..board.app import App
from ..help import TUTORIAL_CARD_HTML, help_button


class StartupController:

    PROGRESS_WIDTH = settings.PROGRESS_WIDTH

    def __init__(self, doc):
        self.doc = doc
        session = doc.session_context
        session_id = session.id if session is not None else uuid.uuid4().hex
        self.SESSION_UPLOAD_DIR = settings.UPLOAD_DIR / session_id
        self.doc.template_variables["layout_settings"] = settings.PANEL_SPLITS

        self.loading = False
        self.preparing = False
        self.failed = False
        self.finished = False
        self._progress_cache = {}
        self._progress_displayed = {}
        self._progress_lock = threading.Lock()
        self._progress_timer = None
        self._ready_to_finish = False
        self._running_workers = 0
        self.dataset = None

        self.image_processor = None
        self.b2s = None

        self._build_widgets()
        self._build_layout()
        self._add_callbacks()
        self._update_source_mode()

        self.tasks = {}

        self.app = App(self.doc)
        self.app.deconvolution_work_dir = self.SESSION_UPLOAD_DIR / "deconvolution"

    def restart(self):
        self.app.reset_geometries()
        self.failed = False
        self.finished = False
        self._ready_to_finish = False
        self.dataset = None
        self.image_processor = None
        self.b2s = None

    def _build_widgets(self):
        self.data_source = RadioButtonGroup(
            labels=["Upload files", "Load example", "ZIP URL"],
            active=0,
            sizing_mode="stretch_width",
        )
        self.zip_url_input = TextInput(
            title="Dataset ZIP URL",
            placeholder="https://example.org/dataset.zip",
        )

        self.h5_input = FileInput(accept=".h5", value="", filename="")
        self.h5_input_widget = column(
            Div(text="<b>1. Expression data (.h5)</b>"),
            self.h5_input
        )

        self.positions_input = FileInput(accept=".parquet", value="", filename="")
        self.positions_input_widget = column(
            Div(text="<b>2. Positions data (.parquet)</b>"),
            self.positions_input
        )

        self.resolution_select = Select(
            title="Resolution (µm)",
            options=[str(size) for size in settings.UPLOAD_BIN_SIZES],
            value=str(settings.UPLOAD_DEFAULT_BIN_SIZE)
        )
        self.resolution_select_widget = column(
            Div(text="<b>3. Resolution</b>"),
            self.resolution_select
        )

        self.image_input = FileInput(accept=".btf, .png, .tiff", value="", filename="")
        self.tiles_progress_bar = Div(text="", width=1000)
        self.image_input_widget = column(
            Div(text="<b>4. Tissue image (.btf, .png, .tiff)</b>"),
            self.image_input
        )

        self.reference_input = FileInput(accept=".h5ad, .h5", value="", filename="")
        self.reference_annotations_input = FileInput(accept=".csv,.parquet", value="", filename="")
        self.reference_upload_status = Div(text=
            "No single-cell reference received yet. After selecting a file, wait for the upload confirmation before loading data.")
        self.reference_input_widget = column(
            Div(text="<b>Single-cell reference (.h5ad or .h5, optional)</b>"
                "<p>Upload once to reuse it for deconvolution in every layer. "
                "AnnData annotations in .obs are loaded automatically. "
                "For a 10x .h5 counts file, supply cell-type annotations in a CSV or Parquet file below.</p>"),
            self.reference_input,
            self.reference_upload_status,
            Div(text="<b>Reference annotations (.csv or .parquet, optional)</b>"
                "<p>First column: cell barcode matching the reference. Parquet may also store "
                "barcodes in its index or a barcode column. "
                "Remaining columns: cell-type labels or other annotations.</p>"),
            self.reference_annotations_input,
        )

        self.pseudospots_checkbox = Checkbox(
            label="Generate Visium v2 pseudospots",
            active=True
        )
        self._build_segmentation_widgets()

        self.example_select = Select(
            title="Choose example dataset",
            options=[(key, key) for key in AVAILABLE_DATASETS.keys()],
            value=list(AVAILABLE_DATASETS.keys())[0]
        )

        self.bins_progress = Div(text="", width=200)
        self.pseudospots_progress = Div(text="", width=200)
        self.status_text = Div(text="Upload your data")

        self.load_button = Button(label="Load Data", button_type="success", disabled=True)
        self.load_button.on_click(self.load_data_callback)

        # Track selection separately: FileInput.value arrives only after the browser
        # has read and transferred the file, including when replacing an old upload.
        self._uploads = [
            (self.h5_input, "Expression data", True, (".h5",)),
            (self.positions_input, "Positions data", True, (".parquet",)),
            (self.image_input, "Tissue image", True, (".btf", ".png", ".tiff")),
            (self.reference_input, "Single-cell reference", False, (".h5ad", ".h5")),
            (self.reference_annotations_input, "Reference annotations", False, (".csv", ".parquet")),
        ]
        self._upload_pending = {widget.id: TextInput(value="", visible=False)
                                for widget, *_ in self._uploads}
        self.upload_status = Div(text="", styles={"overflow_wrap": "anywhere"})

    def _build_segmentation_widgets(self):
        self.segmentation_controls = {}
        try:
            self.segmentation = SegmentationService(**segmentation_kwargs())
        except Exception as exc:
            traceback.print_exc()
            self.segmentation = None
            self.segmentation_info = Div(
                text=f'<span style="color:#fb9292">Segmentation is unavailable: {escape(str(exc))}</span>')
            return

        runner = self.segmentation.runner.name
        hint = ("" if runner != "local" else
                " No container runtime found: plugins use this Python environment "
                "(see segmentation/config.yaml).")
        self.segmentation_info = Div(text=(
            '<p style="color:#b7cbd1;margin:0">Select one or more methods to compare; '
            f'each becomes its own cells layer. Runner: <b>{runner}</b>.{escape(hint)}</p>'))

        for plugin in self.segmentation.plugins.values():
            checkbox = Checkbox(label=plugin.title, active=plugin.name == "stardist")
            params = {param.name: self._param_widget(param) for param in plugin.params}
            about = escape(plugin.description)
            if plugin.url:
                about += f' <a href="{escape(plugin.url)}" target="_blank" style="color:#88f5dc">More</a>'
            box = column(
                Div(text=f'<p style="color:#b7cbd1;margin:0">{about}</p>'),
                *params.values(),
                visible=checkbox.active, sizing_mode="stretch_width", spacing=6,
                styles=dict(padding_left="28px"),
            )
            self.segmentation_controls[plugin.name] = dict(
                plugin=plugin, checkbox=checkbox, params=params, box=box,
                progress=Div(text="", width=200),
            )

    @staticmethod
    def _param_widget(param):
        if param.type in ("float", "int"):
            widget = NumericInput(
                title=param.label, value=param.default, mode=param.type,
                low=param.min, high=param.max,
                placeholder="default" if param.default is None else "",
            )
        elif param.type == "choice":
            widget = Select(title=param.label, value=str(param.default),
                            options=[str(choice) for choice in param.choices])
        else:
            widget = Checkbox(label=param.label, active=bool(param.default))
        if param.help and hasattr(widget, "description"):
            widget.description = param.help
        return widget

    def _read_segmentation_params(self):
        """Selected methods' parameters, validated before loading starts."""
        values = {}
        for name, controls in self.segmentation_controls.items():
            if not controls["checkbox"].active:
                continue
            plugin = controls["plugin"]
            raw = {}
            for param in plugin.params:
                widget = controls["params"][param.name]
                if param.type == "choice":
                    raw[param.name] = next(c for c in param.choices if str(c) == widget.value)
                elif param.type == "bool":
                    raw[param.name] = widget.active
                else:
                    raw[param.name] = widget.value
            try:
                values[name] = plugin.resolve_params(raw)
            except ValueError as exc:
                raise ValueError(f"{plugin.title}: {exc}") from None
        return values

    @property
    def source_mode(self):
        return ("upload", "example", "url")[self.data_source.active]

    def _build_layout(self):
        widget_css = """
            :host { color: #e5edf0; font-family: Arial, sans-serif; }
            .bk-input, .bk-btn { border-radius: 10px; }
            .bk-input { background: #182327; color: #e5edf0;
                border: 1px solid #53666c; padding: 10px; }
            label { color: #c4d3d7; }
            .bk-btn { background: #243338; color: #e5edf0;
                border: 1px solid #53666c; padding: 13px 8px; }
            .bk-btn:hover { background: #30464b; }
            .bk-btn.bk-active, .bk-btn-success { background: #5ee0c2;
                color: #102a24; border-color: #5ee0c2; font-weight: 600; }
            .bk-btn:disabled { opacity: .45; cursor: not-allowed; }
            .bk-btn:focus-visible, .bk-input:focus-visible,
            input:focus-visible { outline: 3px solid #88f5dc; outline-offset: 3px; }
            input[type=checkbox] { accent-color: #5ee0c2; }
        """
        card_style = dict(
            background="#202c31", border="1px solid #3b5057",
            border_radius="18px", padding="20px", box_sizing="border-box",
        )

        def card(*children):
            return column(*children, sizing_mode="stretch_width",
                          spacing=14, styles=card_style)

        self.upload_panel = column(
            self.h5_input_widget, self.positions_input_widget,
            self.resolution_select_widget, self.image_input_widget,
            sizing_mode="stretch_width", spacing=12,
        )
        self.example_panel = column(
            Div(text="Explore a ready-to-load dataset with your preprocessing choices."),
            self.example_select, sizing_mode="stretch_width",
        )
        self.url_panel = column(
            self.zip_url_input,
            Div(text='<span style="color:#88f5dc;font-weight:600">Coming soon</span>'
                '<p>Paste a ZIP URL here. Loading from a URL is not available yet.</p>'),
            sizing_mode="stretch_width",
        )
        def heading(html, topic):
            return row(Div(text=html), help_button(topic), sizing_mode="stretch_width")

        self.preprocessing_panel = card(
            heading('<h2 style="margin:0;font-size:19px">Preprocessing</h2>'
                    '<p style="color:#b7cbd1;margin-bottom:0">Shared across all data sources. '
                    'Choose the representations to generate.</p>', "preprocessing"),
            self.pseudospots_checkbox,
            Div(text='<h3 style="margin:8px 0 0;font-size:16px">Cell segmentation</h3>'),
            self.segmentation_info,
            *[widget for controls in self.segmentation_controls.values()
              for widget in (controls["checkbox"], controls["box"])],
        )
        self.progress_panel = card(
            heading('<h2 style="margin:0;font-size:19px">Loading progress</h2>'
                    '<p style="color:#b7cbd1;margin-bottom:0">Tiles report completed tiles. '
                    'Preprocessing percentages track stages completed, not time remaining.</p>',
                    "progress"),
            *self._progress_divs,
        )
        self.layout = column(
            Div(text='<div style="color:#88f5dc;letter-spacing:2px;font-size:12px;'
                'font-weight:700">SPATIAL TRANSCRIPTOMICS</div>'
                '<h1 style="font-size:32px;margin:10px 0">Visium HD Viewer</h1>'
                '<p style="color:#b7cbd1;font-size:15px">Start with your data. '
                'Bring spatial expression into focus.</p>'),
            card(Div(text=TUTORIAL_CARD_HTML)),
            card(
                heading('<h2 style="margin:0;font-size:19px">Choose your data source</h2>', "data_source"),
                self.data_source, self.upload_panel, self.example_panel, self.url_panel,
                self.reference_input_widget,
            ),
            self.preprocessing_panel,
            self.upload_status, self.load_button, self.status_text, self.progress_panel,
            *self._upload_pending.values(),
            sizing_mode="stretch_width", spacing=20,
            styles=dict(color="#e5edf0"), name="startup",
        )

        for model in self.layout.references():
            if isinstance(model, (FileInput, Select, TextInput, NumericInput, Checkbox,
                                  Button, RadioButtonGroup, Div)):
                model.stylesheets = [widget_css]
                model.sizing_mode = "stretch_width"
                model.width = None

        for widget in (self.h5_input_widget, self.positions_input_widget,
                       self.resolution_select_widget, self.image_input_widget):
            widget.sizing_mode = "stretch_width"

        self.load_button.label = "Load data"

        for div in self._progress_divs:
            div.height = 76
            div.stylesheets = [widget_css]

        self._reset_progress()


    def _update_source_mode(self):
        mode = self.source_mode

        self.upload_panel.visible = mode == "upload"
        self.example_panel.visible = mode == "example"
        self.url_panel.visible = mode == "url"
        self.upload_status.visible = mode == "upload"
        self.load_button.disabled = self.loading or mode == "url"

        bundled = self._example_config().get("reference_path")
        self.reference_input_widget.visible = mode == "upload"
        self.status_text.text = {
            "upload": "Choose your files to get started.",
            "example": "This example will load without a single-cell reference.",
            "url": "URL loading is coming soon.",
        }[mode]
        if bundled:
            self.status_text.text = ("The example's single-cell reference and spatial metadata will load from disk. "
                                     "No upload is required.")
            if self._example_config().get("reference_annotations_path"):
                self.status_text.text += " Single-cell annotations are included for deconvolution."
        self._update_upload_status()

    def _upload_state(self, widget, required, suffixes):
        pending = self._upload_pending[widget.id].value
        name = pending or widget.filename or ""
        if name and not name.lower().endswith(suffixes):
            return False, f"Invalid file: {name}. Expected {' or '.join(suffixes)}."
        if pending or (name and not widget.value):
            return False, f"Uploading {name}… Please wait for server confirmation."
        if widget.value and name:
            return True, f"Uploaded: {name}"
        return not required, "Required — waiting for a file." if required else "Optional — no file selected."

    def _uploads_ready(self):
        return (all(self._upload_state(widget, required, suffixes)[0]
                    for widget, _, required, suffixes in self._uploads)
                and not self._reference_upload_error())

    def _update_upload_status(self):
        rows = []
        for widget, label, required, suffixes in self._uploads:
            _, message = self._upload_state(widget, required, suffixes)
            rows.append(f"<li><b>{escape(label)}:</b> {escape(message)}</li>")
        reference_error = self._reference_upload_error()
        if reference_error:
            rows.append(f"<li>{escape(reference_error)}</li>")
        self.upload_status.text = '<ul style="margin:0;padding-left:20px">' + "".join(rows) + "</ul>"
        if self.loading or self.source_mode != "upload":
            return
        ready = self._uploads_ready()
        self.load_button.disabled = not ready
        self.status_text.text = ("All selected files received. Ready to load." if ready else
                                 "Complete the required uploads and wait for all selected files to arrive.")

    def _example_config(self):
        return AVAILABLE_DATASETS.get(self.example_select.value, {}) if self.source_mode == "example" else {}

    def _add_callbacks(self):
        for index, (widget, label, _, _) in enumerate(self._uploads):
            pending = self._upload_pending[widget.id]
            css_class = f"startup-upload-{index}"
            widget.css_classes = [*widget.css_classes, css_class]

            def received(attr, old, new, pending=pending):
                if new:
                    pending.value = ""
                self._update_upload_status()

            widget.on_change("value", received)
            widget.on_change("filename", lambda attr, old, new: self._update_upload_status())
            pending.on_change("value", lambda attr, old, new: self._update_upload_status())
            # Listen before Bokeh's asynchronous FileReader starts. Shadow roots
            # must be traversed because widgets render inside their own roots.
            self.doc.js_on_event(DocumentReady, CustomJS(args=dict(
                widget=widget, pending=pending, button=self.load_button,
                status=self.status_text, source=self.data_source,
                css_class=css_class, label=label,
            ), code="""
                function findHost(root) {
                    const match = root.querySelector('.' + css_class);
                    if (match != null) return match;
                    for (const element of root.querySelectorAll('*')) {
                        if (element.shadowRoot != null) {
                            const found = findHost(element.shadowRoot);
                            if (found != null) return found;
                        }
                    }
                    return null;
                }
                const host = findHost(document);
                const input = host?.shadowRoot?.querySelector('input[type=file]')
                    ?? host?.querySelector('input[type=file]');
                if (input == null) return;
                input.addEventListener('change', () => {
                    const file = input.files?.[0];
                    // Clear the previous payload so even a replacement with
                    // identical contents produces a new server confirmation.
                    widget.setv({value: '', filename: ''});
                    pending.value = file?.name ?? '';
                    if (source.active === 0) {
                        button.disabled = true;
                        status.text = file == null ? 'Choose a file to continue.'
                            : 'Uploading ' + label + '… Please wait for server confirmation.';
                    }
                }, {capture: true});
            """))

        def reference_received(attr, old, new):
            if new:
                name = escape(self.reference_input.filename or "Single-cell reference")
                self.reference_upload_status.text = f"Reference uploaded: {name}. Received by server."
            else:
                self.reference_upload_status.text = "No single-cell reference received yet. Wait for upload confirmation."

        self.reference_input.on_change("value", reference_received)
        self.data_source.on_change(
            "active", lambda attr, old, new: self._update_source_mode()
        )
        self.example_select.on_change("value", lambda attr, old, new: self._update_source_mode())

        self.pseudospots_checkbox.on_change("active", lambda attr, old, new: self._reset_progress())
        for controls in self.segmentation_controls.values():
            def on_toggle(attr, old, new, box=controls["box"]):
                box.visible = new
                self._reset_progress()
            controls["checkbox"].on_change("active", on_toggle)

        # La sección Cells de index.html empieza oculta: mostrarla solo si se
        # han segmentado células.
        self.status_text.js_on_change("text", CustomJS(
            args=dict(checkboxes=[c["checkbox"] for c in self.segmentation_controls.values()]),
            code="""
            if (cb_obj.text === "Done") {
                showMainApp();
                const cells_section = document.getElementById('cells-section');
                if (cells_section && checkboxes.some((cb) => cb.active)) {
                    cells_section.style.display = 'block';
                }
            }
        """))

    def add_to_document(self):
        self.doc.add_root(self.layout)

    @property
    def _progress_divs(self):
        return (self.tiles_progress_bar, self.bins_progress, self.pseudospots_progress,
                *[c["progress"] for c in self.segmentation_controls.values()])

    def _progress_html(self, div, pct=0, status="Waiting"):
        label = self._task_label(next(task for task, d in self._task_divs.items() if d is div))
        pct = max(0, min(100, int(pct)))
        detail = escape(status) + f" · {pct}%"
        color = "#fb9292" if status.startswith("Failed") else "#5ee0c2"
        aria = f'aria-valuenow="{pct}"'
        return f"""
        <div style="width:{self.PROGRESS_WIDTH}px;max-width:100%;color:#e5edf0;font-size:14px;line-height:22px">
          <strong>{label}</strong>
          <div style="color:#b7cbd1;font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis">{detail}</div>
          <div role="progressbar" aria-label="{label}" aria-valuemin="0"
               aria-valuemax="100" {aria} aria-valuetext="{detail}"
               style="height:10px;margin-top:8px;background:#394c52;border-radius:5px;overflow:hidden">
            <div style="width:{pct}%;height:100%;background:{color};border-radius:5px"></div>
          </div>
        </div>"""

    def _reset_progress(self):
        with self._progress_lock:
            self._progress_cache.clear()
        self._progress_displayed.clear()

        enabled_by_div = [(self.tiles_progress_bar, True),
                          (self.bins_progress, True),
                          (self.pseudospots_progress, self.pseudospots_checkbox.active)]
        enabled_by_div += [(c["progress"], c["checkbox"].active)
                           for c in self.segmentation_controls.values()]

        for div, enabled in enabled_by_div:
            div.text = self._progress_html(div, status="Waiting" if enabled else "Skipped")

    def update_progress(self, div, pct=None, status="Processing"):
        with self._progress_lock:
            previous_pct = self._progress_cache.get(div.id, (0, ""))[0]
            pct = previous_pct if pct is None else max(previous_pct, min(100, int(pct)))
            self._progress_cache[div.id] = (pct, status)

    def _render_progress(self):
        with self._progress_lock:
            targets = dict(self._progress_cache)
        settled = True
        for div in self._progress_divs:
            if div.id not in targets:
                continue
            target, status = targets[div.id]
            shown = self._progress_displayed.get(div.id, 0)
            # Ease toward actual reported work, never ahead of it. At 20 fps,
            # large jumps settle in about a second without delaying computation.
            shown = target if status.startswith("Failed") else min(
                target, shown + max(1, (target - shown + 4) // 5)
            )
            self._progress_displayed[div.id] = shown
            settled = settled and shown == target
            label = "Finishing" if status == "Complete" and shown < 100 else status
            html = self._progress_html(div, shown, label)
            if div.text != html:
                div.text = html
        if settled and self._ready_to_finish:
            self._finish_loading()

    @property
    def _task_divs(self):
        divs = {"tiles": self.tiles_progress_bar, "bins": self.bins_progress,
                "spots": self.pseudospots_progress}
        divs.update({f"cells:{name}": c["progress"] for name, c in self.segmentation_controls.items()})
        return divs

    def _task_progress(self, task):
        return self._task_divs[task]

    TASK_LABELS = {"tiles": "Image tiles", "bins": "Bins", "spots": "Pseudospots"}

    def _task_label(self, task):
        if task and task.startswith("cells:"):
            return f"Cells · {self.segmentation_controls[task[6:]]['plugin'].title}"
        return self.TASK_LABELS.get(task, "Loading")

    def _report_failure(self, task, exc=None, fatal=True):
        """A fatal failure stops loading; others (one segmentation method) only
        mark their own progress bar so the remaining data still loads."""
        if fatal:
            self.failed = True
        name = self._task_label(task)
        reason = None
        if exc is not None:
            # Workers run in threads: without this the traceback is lost
            print(f"[startup] {name} failed:", file=sys.stderr)
            traceback.print_exception(type(exc), exc, exc.__traceback__, file=sys.stderr)
            sys.stderr.flush()
            reason = f"{type(exc).__name__}: {exc}"
        if task:
            self.update_progress(self._task_progress(task),
                                 status=f"Failed · {reason}" if reason else "Failed")
        if not fatal:
            return
        message = f"{name} failed: {reason}" if reason else "Loading failed."
        self.doc.add_next_tick_callback(lambda: setattr(
            self.status_text, "text", escape(message) + " Change the settings and try again."
        ))

    def _save_uploaded_file(self, file_input):
        self.SESSION_UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
        path = self.SESSION_UPLOAD_DIR / file_input.filename
        with open(path, "wb") as f:
            f.write(base64.b64decode(file_input.value))
        return str(path)

    def _validate_inputs(self):
        if not self.h5_input.value:
            self.status_text.text = "Please upload an .h5 file"
            return False

        if not self.positions_input.value:
            self.status_text.text = "Please upload a .parquet positions file"
            return False

        if not self.image_input.value:
            self.status_text.text = "Please upload a tissue image (.btf, .png, .tiff)"
            return False

        return True
    
    @property
    def _loading_controls(self):
        return (self.data_source, self.load_button, self.pseudospots_checkbox,
                self.example_select, self.h5_input,
                self.positions_input, self.image_input, self.resolution_select,
                self.reference_input, self.reference_annotations_input,
                *[widget for controls in self.segmentation_controls.values()
                  for widget in (controls["checkbox"], *controls["params"].values())])

    def _restore_controls_after_failure(self):
        if not self.failed or self.preparing or self._running_workers or not self.loading:
            return
        self.loading = False
        self._ready_to_finish = False
        self._render_progress()
        if self._progress_timer is not None:
            self.doc.remove_periodic_callback(self._progress_timer)
            self._progress_timer = None
        for widget in self._loading_controls:
            widget.disabled = False
        # Reapply the readiness rules while retaining the useful failure message.
        error = self.status_text.text
        self._update_source_mode()
        self.status_text.text = error

    def _check_all_done(self):
        if self.failed:
            if not self.preparing and not self._running_workers:
                self.doc.add_next_tick_callback(self._restore_controls_after_failure)
            return
        if not self.preparing and not self.failed and all(self.tasks.values()):
            self.doc.add_next_tick_callback(self._all_tasks_finished)
    
    def run_in_thread(self, func, sample_type=None, task=None, sample_kwargs=None, fatal=True):
        task = task or sample_type
        self._running_workers += 1

        def finish(error=None):
            if error is not None:
                self._report_failure(task, error, fatal=fatal)
            elif task in self.tasks:
                self.update_progress(self._task_progress(task), 100, "Complete")
            if task in self.tasks:
                self.tasks[task] = True
            self._running_workers -= 1
            self._check_all_done()

        def wrapper():
            try:
                output = func()
            except Exception as exc:
                finish(exc)
                return

            if not sample_type or output is None:
                finish()
                return

            # The task is only complete once the sample is in the viewer: adding it runs on the document thread and can fail too
            def update_doc():
                try:
                    self.app.add_sample(*output, self.b2s, sample_type=sample_type,
                                        **(sample_kwargs or {}))
                except Exception as exc:
                    finish(exc)
                else:
                    finish()

            self.doc.add_next_tick_callback(update_doc)
        
        threading.Thread(target=wrapper).start()

    def generate_tiles(self, base_dir):
        from .data.image.image import ImageProcessor
        from .data.image.transform import SpatialTransform
        from ..board.components.image import ImageComponent

        if self.dataset.df_tissue_positions is None:
            self.dataset.read_positions()

        self.spatial_transform = SpatialTransform.fit(
            self.dataset.df_tissue_positions,
            self.dataset.image_shape,
            self.dataset.image_scale,
            **grid_fit_kwargs(),
        )

        cache_key = f"padded-v2-{self._image_key()}-{self.spatial_transform.fingerprint}"
        dataset = self.dataset
        def load_image():
            if dataset.img is None:
                dataset.read_img()
            return dataset.img

        self.image_processor = ImageProcessor(
            dataset.img, 
            base_dir / 'tiles' / cache_key,
            spatial_transform=self.spatial_transform, 
            image_loader=load_image,
            **image_kwargs(),
        )

        single_image = self.image_processor.mem <= ImageComponent.TILE_MEM_THRESHOLD

        def worker():
            def on_progress(done, total):
                self.update_progress(self.tiles_progress_bar, int(100 * done / total) if total else 0,
                                     "Creating padded image tiles")
            return self.image_processor.create_tiles(overwrite=False, progress_callback=on_progress,
                                                     single_image=single_image)

        generated = worker()

        self.doc.add_next_tick_callback(lambda: self.app.add_image_component(self.image_processor))
        self.update_progress(self.tiles_progress_bar, 100, "Complete" if generated else "Using cached tiles")
        self.tasks['tiles'] = True


    def build_b2s(self, adata, df_tissue_positions, bins_size):
        from .operations.bins2spots import Bins2Spots

        def on_progress(step, done, total):
            # Three equally weighted stages: preparation, geometry, aggregation.
            pct = 33 + int(33 * done / total) if total else 33
            if self.pseudospots_checkbox.active:
                self.update_progress(self.pseudospots_progress, pct, "Building spot geometry")

        self.b2s = Bins2Spots(
            adata,
            df_tissue_positions,
            bins_size=bins_size,
            pixels_per_bin=self.dataset.pixels_per_bin,
            spatial_transform=self.spatial_transform,
            tile_scale_xfactor=self.image_processor.scale_xfactor * self.dataset.image_scale,
            tile_scale_yfactor=self.image_processor.scale_yfactor * self.dataset.image_scale,
            progress_callback=on_progress,
            **bins_to_spots_kwargs(),
        )

        self.update_progress(self.bins_progress, 50, "Building bin geometry")
        bin_squares = self.b2s.bin_squares
        self.update_progress(self.bins_progress, 75, "Adding bins to the viewer")

        def add_bins():
            try:
                self.app.add_sample(bin_squares, self.b2s.adata_hd, self.b2s, sample_type="bins")
                self.update_progress(self.bins_progress, 100, "Complete")
            except Exception as exc:
                self._report_failure("bins", exc)
            finally:
                self.tasks["bins"] = True
                self._check_all_done()

        self.doc.add_next_tick_callback(add_bins)

    def ensure_b2s(self):
        if self.b2s is None:

            if self.dataset.adata is None or self.dataset.df_tissue_positions is None:
                self.update_progress(self.bins_progress, status="Reading expression data")
                if self.dataset.adata is None:
                    self.dataset.read_adata()
                if self.dataset.df_tissue_positions is None:
                    self.dataset.read_positions()

            annotations_path = self._example_config().get("spatial_annotations_path")
            if annotations_path:
                from .data.singlecell import add_annotations
                self.update_progress(self.bins_progress, status="Reading spatial annotations")
                add_annotations(self.dataset.adata, annotations_path)

            self.update_progress(self.bins_progress, 25, "Preparing bins")
            self.build_b2s(self.dataset.adata, self.dataset.df_tissue_positions, self.dataset.bins_size)

    def generate_pseudospots(self):
        self.update_progress(self.pseudospots_progress, status="Preparing expression data")
        self.ensure_b2s()
        if self.b2s is None:
            return

        def worker():
            self.update_progress(self.pseudospots_progress, 33, "Building spot geometry")
            pseudospots_polygons = self.b2s.spot_polygons
            self.update_progress(self.pseudospots_progress, 66, "Aggregating expression")
            pseudospots_adata = self.b2s.convert(normalize=False)
            return pseudospots_polygons, pseudospots_adata

        self.run_in_thread(worker, sample_type="spots")

    def _image_key(self):
        """Identifies the image file, so segmentation results are cached per image."""
        path = Path(self.dataset.img_path).resolve()
        stat = path.stat()
        return hashlib.sha1(f"{path}|{stat.st_size}|{stat.st_mtime_ns}".encode()).hexdigest()[:12]

    def generate_segmentations(self, base_dir):
        from .operations.pooler import CellsPooler

        selected = [c for c in self.segmentation_controls.values() if c["checkbox"].active]
        if not selected:
            return
        if self.image_processor is None:
            raise RuntimeError("The tissue image must be loaded before segmenting cells")

        for controls in selected:
            self.update_progress(controls["progress"], status="Preparing image")
        self.ensure_b2s()

        work_dir = Path(base_dir) / "segmentation"
        image_key = self._image_key()
        # OpenCV (png/tif) reads BGR; plugins always receive RGB
        if self.dataset.img is None:
            self.dataset.read_img()

        image_path = self.segmentation.prepare_image(
            self.dataset.img, 
            image_key, 
            work_dir, 
            bgr=self.dataset.img_format in ("png", "tif", "tiff")
        )

        # Plugins see the given image: µm per pixel of that image (not of the full-res one)
        um_per_px = (float(self.dataset.bins_size) / (self.dataset.pixels_per_bin * self.dataset.image_scale)
                     if self.dataset.pixels_per_bin else None)
        # Segment the original image; transform only the resulting outlines.

        for controls in selected:
            plugin = controls["plugin"]

            def worker(plugin=plugin, div=controls["progress"]):
                def on_progress(fraction, message):
                    self.update_progress(div, 5 + int(80 * fraction), message or "Segmenting")

                cells, out_dir = self.segmentation.run(
                    plugin.name, image_path, image_key, work_dir,
                    params=self._segmentation_params[plugin.name],
                    um_per_px=um_per_px, on_progress=on_progress,
                )
                print(f"[segmentation] {plugin.title}: {len(cells)} cells ({out_dir})")
                if cells.empty:
                    raise RuntimeError(f"{plugin.title} found no cells")

                cells.geometry = self.spatial_transform.image_geometry(cells.geometry.values)
                cells["index"] = [f"{plugin.name}_{i + 1}" for i in range(len(cells))]
                cells.index = cells["index"].values

                self.update_progress(div, 85, f"Aggregating expression of {len(cells)} cells")
                cells_adata = CellsPooler(
                    self.b2s.bin_squares, self.b2s.adata_hd, cells, keep_cols=[]
                ).run(predicate="intersects")
                self.update_progress(div, 95, "Adding cells to the viewer")

                return cells, cells_adata

            self.run_in_thread(
                worker, task=f"cells:{plugin.name}", sample_type="cells",
                sample_kwargs=dict(method=plugin.name, title=plugin.title),
                fatal=False,
            )

    def set_tasks(self):
        # False = pending; tasks that were not selected count as done
        self.tasks = {
            "tiles": False,
            "bins": False,
            "spots": not self.pseudospots_checkbox.active,
        }
        for name, controls in self.segmentation_controls.items():
            self.tasks[f"cells:{name}"] = not controls["checkbox"].active

    def load_data_callback(self):
        if self.loading or self.source_mode == "url":
            return

        reference_error = self._reference_upload_error()
        if reference_error:
            # Leave the form enabled so the pending upload can finish and the user can retry.
            self.status_text.text = reference_error
            return

        if self.source_mode == "upload":
            self._update_upload_status()
            if not self._uploads_ready():
                return
            
        if self.source_mode == "example" and self.example_select.value not in AVAILABLE_DATASETS:
            self.status_text.text = "Please choose an example dataset"
            return

        try:
            self._segmentation_params = self._read_segmentation_params()
        except ValueError as exc:
            self.status_text.text = escape(str(exc))
            return

        self.restart()
        self.set_tasks()
        self.status_text.text = "Loading your dataset…"
        self.loading = True
        self.preparing = True
        self._reset_progress()
        self._progress_timer = self.doc.add_periodic_callback(self._render_progress, settings.STARTUP_PROGRESS_REFRESH_MS)
        for widget in self._loading_controls:
            widget.disabled = True

        self.update_progress(self.tiles_progress_bar, status="Checking image metadata and cached tiles")

        def prepare():
            try:
                self._prepare_data()
            except Exception as exc:
                self._report_failure(None, exc)
                for task, done in self.tasks.items():
                    if not done:
                        self.update_progress(self._task_progress(task), status="Failed")
            finally:
                self.preparing = False
                self._check_all_done()

        threading.Thread(target=prepare, daemon=True).start()

    def _reference_upload_error(self):
        if self.source_mode != "upload":
            return None
        if not self.reference_input.value and (self.reference_annotations_input.value or self.reference_input.filename):
            return ("The single-cell reference file has not reached the server yet. "
                    "Wait for the 'Reference uploaded' confirmation, then click Load data again. "
                    "If it never appears, check the Bokeh server's WebSocket upload-size limit.")
        return None

    def _prepare_deconvolution_reference(self):
        """Read the optional upload once, before any layer is constructed (worker thread)."""
        from .data.singlecell import read_reference

        reference = None
        name = None
        reference_error = self._reference_upload_error()
        if reference_error:
            raise ValueError(reference_error)
        config = self._example_config()
        if config.get("reference_path"):
            name = Path(config["reference_path"]).name
            try:
                reference = read_reference(config["reference_path"], config.get("reference_annotations_path"))
            except Exception as exc:
                raise ValueError(f"Cannot load example single-cell reference {name}: {exc}") from exc
        elif self.source_mode == "upload" and self.reference_input.value:
            name = Path(self.reference_input.filename or "reference.h5ad").name
            try:
                annotations = (base64.b64decode(self.reference_annotations_input.value, validate=True)
                               if self.reference_annotations_input.value else None)
                reference = read_reference(
                    base64.b64decode(self.reference_input.value, validate=True), annotations)
            except Exception as exc:
                raise ValueError(f"Cannot load single-cell reference {name}: {exc}") from exc
        self.app.deconvolution_reference = reference
        self.app.deconvolution_reference_name = name
        self.app.deconvolution_label_column = config.get("reference_label_column", "")

    def _prepare_data(self):
        self._prepare_deconvolution_reference()
        if self.source_mode == "upload":

            adata_path = self._save_uploaded_file(self.h5_input)
            positions_path = self._save_uploaded_file(self.positions_input)
            img_path = self._save_uploaded_file(self.image_input)

            bins_size = int(self.resolution_select.value)

        else:
            img_path = AVAILABLE_DATASETS[self.example_select.value]["image_path"]
            adata_path = AVAILABLE_DATASETS[self.example_select.value]["expression_path"]
            positions_path = AVAILABLE_DATASETS[self.example_select.value]["positions_path"]
            bins_size = AVAILABLE_DATASETS[self.example_select.value]["bins_size"]

        self.dataset = VisiumHD(
            img_path=img_path,
            adata_path=adata_path,
            positions_path=positions_path,
            bins_size=bins_size
        )

        base_dir = (
            self.SESSION_UPLOAD_DIR if self.source_mode == "upload"
            else Path(AVAILABLE_DATASETS[self.example_select.value]["image_path"]).resolve().parent
        )
        # Use the same dataset base directory as segmentation, before generators
        # add bins/spots/cell layers (including layers arriving asynchronously).
        self.app.deconvolution_work_dir = Path(base_dir) / 'deconvolution'
        self.generate_tiles(base_dir)
        
        if self.pseudospots_checkbox.active:
            self.generate_pseudospots()
        
        self.generate_segmentations(base_dir)

        # Bins are always displayed; the generators above only build them as a
        # dependency, so with both options unchecked they would never be added
        self.ensure_b2s()

    def _all_tasks_finished(self):
        if self.failed or self.finished:
            return
        self._ready_to_finish = True
        self._render_progress()

    def _finish_loading(self):
        if self.failed or self.finished:
            return
        if self._progress_timer is not None:
            self.doc.remove_periodic_callback(self._progress_timer)
            self._progress_timer = None
        try:
            self.app.build()
        except Exception as exc:
            self._report_failure(None, exc)
            self._check_all_done()
            return
        self.finished = True
        self.status_text.text = "Done"


def create_startup(doc):
    controller = StartupController(doc)
    controller.add_to_document()

def on_session_destroyed(session_context):
    upload_dir = settings.UPLOAD_DIR / session_context.id
    if upload_dir.exists():
        try:
            shutil.rmtree(upload_dir)
            print(f"[startup] Removed upload directory {upload_dir}")
        except Exception as exc:
            print(f"[startup] Failed to remove upload directory {upload_dir}: {exc}", file=sys.stderr)
