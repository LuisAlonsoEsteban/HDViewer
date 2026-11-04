
from ..settings.config import grid_fit_kwargs, image_kwargs, bins_to_spots_kwargs, segmentation_kwargs

from .. import settings
from ..presentation import configure_document
from bokeh.models import (
    Button, Div, Select,
    Checkbox, CustomJS, RadioButtonGroup, TextInput, NumericInput,
)
from bokeh.layouts import row, column

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
from ..plugins.jobs import JOBS
from .datasets import AVAILABLE_DATASETS
from .uploads import ChunkedUpload

from ..help import TUTORIAL_CARD_HTML, help_button


def __getattr__(name):
    # Board modules import startup.utils, which runs this package: a module-level
    # import of the board would make `import hd_viewer.board...` circular.
    if name == "App":
        from ..board.app import App
        return App
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


class StartupController:

    PROGRESS_WIDTH = settings.PROGRESS_WIDTH

    def __init__(self, doc):
        self.doc = doc
        configure_document(doc)
        session = doc.session_context
        session_id = session.id if session is not None else uuid.uuid4().hex
        self.session_id = session_id
        self.SESSION_UPLOAD_DIR = session_upload_dir(session_id)
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

        self.app = sys.modules[__name__].App(self.doc)
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
            labels=["Upload files", "Load example"],
            active=0,
            sizing_mode="stretch_width",
        )

        # Files stream to the session folder in chunks (see uploads.py), so
        # multi-GB uploads fit in Bokeh's default websocket message size.
        def upload_field(key, suffixes):
            return ChunkedUpload(key, suffixes, lambda: self.SESSION_UPLOAD_DIR,
                                 on_change=self._update_upload_status,
                                 can_upload=lambda: not self.loading)

        self.h5_input = upload_field("expression", (".h5",))
        self.h5_input_widget = column(
            Div(text="<b>1. Expression data (.h5)</b>"),
            self.h5_input.widget
        )

        self.positions_input = upload_field("positions", (".parquet",))
        self.positions_input_widget = column(
            Div(text="<b>2. Positions data (.parquet)</b>"),
            self.positions_input.widget
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

        self.image_input = upload_field("image", (".btf", ".png", ".tif", ".tiff"))
        self.tiles_progress_bar = Div(text="", width=1000)
        self.image_input_widget = column(
            Div(text="<b>4. Tissue image (.btf, .png, .tiff)</b>"),
            self.image_input.widget
        )

        self.reference_input = upload_field("reference", (".h5ad", ".h5"))
        self.reference_annotations_input = upload_field("reference_annotations", (".csv", ".parquet"))
        self.reference_upload_status = Div(text=
            "No single-cell reference received yet. After selecting a file, wait for the upload confirmation before loading data.")
        self.reference_input_widget = column(
            Div(text="<b>Single-cell reference (.h5ad or .h5, optional)</b>"
                "<p>Upload once to reuse it for deconvolution in every layer. "
                "AnnData annotations in .obs are loaded automatically. "
                "For a 10x .h5 counts file, supply cell-type annotations in a CSV or Parquet file below.</p>"),
            self.reference_input.widget,
            self.reference_upload_status,
            Div(text="<b>Reference annotations (.csv or .parquet, optional)</b>"
                "<p>First column: cell barcode matching the reference. Parquet may also store "
                "barcodes in its index or a barcode column. "
                "Remaining columns: cell-type labels or other annotations.</p>"),
            self.reference_annotations_input.widget,
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

        # A file counts only once every chunk has reached the server; choosing a
        # replacement discards the previous file and blocks loading until it arrives.
        self._uploads = [
            (self.h5_input, "Expression data", True),
            (self.positions_input, "Positions data", True),
            (self.image_input, "Tissue image", True),
            (self.reference_input, "Single-cell reference", False),
            (self.reference_annotations_input, "Reference annotations", False),
        ]
        self.upload_status = Div(text="", styles={"overflow_wrap": "anywhere"})

    def _build_segmentation_widgets(self):
        self.segmentation_controls = {}
        try:
            self.segmentation = segmentation_service()
        except Exception as exc:
            traceback.print_exc()
            self.segmentation = None
            self.segmentation_info = Div(
                text=f'<span style="color:#fb9292">Segmentation is unavailable: {escape(str(exc))}</span>')
            return

        self.segmentation_info = Div(text=self._segmentation_info(self.segmentation.runner_name))
        if self.segmentation.runner_name is None:
            # Detection may probe the Docker daemon for seconds: never on the IO loop
            def detected(runner, info=self.segmentation_info, doc=self.doc):
                text = self._segmentation_info(runner)
                try:
                    doc.add_next_tick_callback(lambda: setattr(info, "text", text))
                except Exception:
                    pass  # session already closed
            self.segmentation.detect_runner(detected)

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
    def _segmentation_info(runner):
        """`runner`: its name, a Runner, the detection error, or None while detecting."""
        if isinstance(runner, Exception):
            return f'<span style="color:#fb9292">Segmentation is unavailable: {escape(str(runner))}</span>'
        runner = getattr(runner, "name", runner)
        hint = ("" if runner != "local" else
                " No container runtime found: plugins use this Python environment "
                "(see segmentation/config.yaml).")
        return ('<p style="color:#b7cbd1;margin:0">Select one or more methods to compare; '
                f'each becomes its own cells layer. Runner: <b>{runner or "detecting…"}</b>.{escape(hint)}</p>')

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
        return ("upload", "example")[self.data_source.active]

    def _build_layout(self):
        widget_css = self.doc._visium_stylesheet
        card_style = dict(
            background="#192c28", border="1px solid #344d45",
            border_radius="14px", padding="24px", box_sizing="border-box",
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
            Div(text='<div style="color:#a5dbc3;letter-spacing:2px;font-size:11px;'
                'font-weight:700">YOUR SPATIAL WORKSPACE</div>'
                '<h1 style="font-size:32px;letter-spacing:-1px;margin:10px 0">Start with a tissue.</h1>'
                '<p style="color:#b7cbd1;font-size:15px">Choose your data, configure the layers, '
                'and open the interactive workspace.</p>'),
            card(Div(text=TUTORIAL_CARD_HTML)),
            card(
                heading('<h2 style="margin:0;font-size:19px">Choose your data source</h2>', "data_source"),
                self.data_source, self.upload_panel, self.example_panel,
                self.reference_input_widget,
            ),
            self.preprocessing_panel,
            self.upload_status, self.load_button, self.status_text, self.progress_panel,
            sizing_mode="stretch_width", spacing=20,
            styles=dict(color="#e5edf0"), name="startup",
        )

        for model in self.layout.references():
            if isinstance(model, (Select, TextInput, NumericInput, Checkbox,
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
        self.upload_status.visible = mode == "upload"
        self.load_button.disabled = self.loading

        bundled = self._example_config().get("reference_path")
        self.reference_input_widget.visible = mode == "upload"
        self.status_text.text = {
            "upload": "Choose your files to get started.",
            "example": "This example will load without a single-cell reference.",
        }[mode]
        if bundled:
            self.status_text.text = ("The example's single-cell reference and spatial metadata will load from disk. "
                                     "No upload is required.")
            if self._example_config().get("reference_annotations_path"):
                self.status_text.text += " Single-cell annotations are included for deconvolution."
        self._update_upload_status()

    @staticmethod
    def _upload_state(field, required):
        ready = field.complete or (field.status == field.EMPTY and not required)
        return ready, field.describe(required)

    def _uploads_ready(self):
        return (all(self._upload_state(field, required)[0]
                    for field, _, required in self._uploads)
                and not self._reference_upload_error())

    def _update_upload_status(self):
        self._update_reference_status()
        rows = []
        for field, label, required in self._uploads:
            _, message = self._upload_state(field, required)
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

    def _update_reference_status(self):
        field = self.reference_input
        if field.complete:
            text = f"Reference uploaded: {escape(field.filename)}. Received by server."
        elif field.status == field.UPLOADING:
            text = escape(field.describe(False))
        else:
            text = "No single-cell reference received yet. Wait for upload confirmation."
        self.reference_upload_status.text = text

    def _add_callbacks(self):
        for field, *_ in self._uploads:
            field.bind_js(self.load_button, self.data_source)

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

    @staticmethod
    def _uploaded_path(field):
        """Path of a fully received upload (already assembled in the session folder)."""
        if not field.complete:
            raise ValueError(f"{field.filename or 'A required file'} has not finished uploading")
        return str(field.path)
    
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
                    um_per_px=um_per_px, on_progress=on_progress, owner=self.session_id,
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
        if self.loading:
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
        reference = self.reference_input
        if not reference.complete and (self.reference_annotations_input.complete or reference.filename):
            return ("The single-cell reference file has not reached the server yet. "
                    "Wait for the 'Reference uploaded' confirmation, then click Load data again.")
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
        elif self.source_mode == "upload" and self.reference_input.complete:
            name = self.reference_input.filename
            try:
                annotations = (self._uploaded_path(self.reference_annotations_input)
                               if self.reference_annotations_input.complete else None)
                reference = read_reference(self._uploaded_path(self.reference_input), annotations)
            except Exception as exc:
                raise ValueError(f"Cannot load single-cell reference {name}: {exc}") from exc
        self.app.deconvolution_reference = reference
        self.app.deconvolution_reference_name = name
        self.app.deconvolution_label_column = config.get("reference_label_column", "")

    def _prepare_data(self):
        self._prepare_deconvolution_reference()
        if self.source_mode == "upload":

            adata_path = self._uploaded_path(self.h5_input)
            positions_path = self._uploaded_path(self.positions_input)
            img_path = self._uploaded_path(self.image_input)

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
        self.app.deconvolution_bin_size_um = bins_size
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

def session_upload_dir(session_id):
    # Session ids can be chosen by the client (?bokeh-session-id=...), so never use
    # one as a path: a hash cannot escape UPLOAD_DIR.
    return settings.UPLOAD_DIR / hashlib.sha256(session_id.encode()).hexdigest()[:32]


_segmentation_service = None
_segmentation_lock = threading.Lock()


def segmentation_service():
    """One queue (max_parallel) and runner per viewer process, shared by all sessions."""
    global _segmentation_service
    with _segmentation_lock:
        if _segmentation_service is None:
            _segmentation_service = SegmentationService(**segmentation_kwargs())
        return _segmentation_service


def on_session_destroyed(session_context):
    # Stopping containers (docker rm -f) may take seconds: not on the IO loop
    threading.Thread(target=_cleanup_session, args=(session_context.id,), daemon=True).start()


def _cleanup_session(session_id):
    # Stop its segmentation/deconvolution runs, unless another session awaits the same one,
    # before removing the folders they write to
    cancelled = JOBS.cancel_owner(session_id)
    if cancelled:
        print(f"[startup] Cancelled {cancelled} analysis job(s) of a closed session")
    upload_dir = session_upload_dir(session_id)
    if upload_dir.exists():
        try:
            shutil.rmtree(upload_dir)
            print(f"[startup] Removed upload directory {upload_dir}")
        except Exception as exc:
            print(f"[startup] Failed to remove upload directory {upload_dir}: {exc}", file=sys.stderr)
