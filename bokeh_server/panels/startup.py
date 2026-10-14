from bokeh.models import (
    FileInput, Button, Div, Select,
    Checkbox, CustomJS, RadioButtonGroup, TextInput,
)
from bokeh.layouts import row, column

import base64
import os
import sys
from pathlib import Path
from html import escape

import threading
import traceback

from .app.app import App
from ..datasets import AVAILABLE_DATASETS

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.load_data import DataSet


class StartupController:

    UPLOAD_DIR = Path("./tmp")

    TILE_MEM_THRESHOLD = 0.3 # GB
    PROGRESS_WIDTH = 560  # px; fixed so the bars do not resize as their text changes

    def __init__(self, doc):
        self.doc = doc

        self.loading = False
        self.preparing = False
        self.failed = False
        self.finished = False
        self._progress_cache = {}
        self._progress_displayed = {}
        self._progress_lock = threading.Lock()
        self._progress_timer = None
        self._ready_to_finish = False
        self.dataset = None

        self.image_processor = None
        self.b2s = None

        self._build_widgets()
        self._build_layout()
        self._add_callbacks()
        self._update_source_mode()

        self.tasks = {
            "tiles": False,
            "bins": False,
            "spots": False,
            "cells": False
        }

        self.app = App(self.doc)

    def restart(self):
        self.app.reset_geometries()

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

        self.h5_input = FileInput(accept=".h5", value="")
        self.h5_input_widget = column(
            Div(text="<b>1. Expression data (.h5)</b>"),
            self.h5_input
        )

        self.positions_input = FileInput(accept=".parquet", value="")
        self.positions_input_widget = column(
            Div(text="<b>2. Positions data (.parquet)</b>"),
            self.positions_input
        )

        self.resolution_select = Select(
            title="Resolution (µm)",
            options=["2", "8", "16"],
            value="16"
        )
        self.resolution_select_widget = column(
            Div(text="<b>3. Resolution</b>"),
            self.resolution_select
        )

        self.image_input = FileInput(accept=".btf, .png, .tiff", value="")
        self.tiles_progress_bar = Div(text="", width=1000)
        self.image_input_widget = column(
            Div(text="<b>4. Tissue image (.btf, .png, .tiff)</b>"),
            self.image_input
        )

        self.pseudospots_checkbox = Checkbox(
            label="Generate Visium v2 pseudospots",
            active=True
        )
        self.cells_checkbox = Checkbox(
            label="Segment cells",
            active=True
        )

        self.example_select = Select(
            title="Choose example dataset",
            options=[(key, key) for key in AVAILABLE_DATASETS.keys()],
            value=list(AVAILABLE_DATASETS.keys())[0]
        )

        self.bins_progress = Div(text="", width=200)
        self.pseudospots_progress = Div(text="", width=200)
        self.cells_progress = Div(text="", width=200)
        self.status_text = Div(text="Upload your data")

        self.load_button = Button(label="Load Data", button_type="success")
        self.load_button.on_click(self.load_data_callback)

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
        self.preprocessing_panel = card(
            Div(text='<h2 style="margin:0;font-size:19px">Preprocessing</h2>'
                '<p style="color:#b7cbd1;margin-bottom:0">Shared across all data sources. '
                'Choose the representations to generate.</p>'),
            self.pseudospots_checkbox,
            self.cells_checkbox,
        )
        self.progress_panel = card(
            Div(text='<h2 style="margin:0;font-size:19px">Loading progress</h2>'
                '<p style="color:#b7cbd1;margin-bottom:0">Tiles report completed tiles. '
                'Preprocessing percentages track stages completed, not time remaining.</p>'),
            *self._progress_divs,
        )
        self.layout = column(
            Div(text='<div style="color:#88f5dc;letter-spacing:2px;font-size:12px;'
                'font-weight:700">SPATIAL TRANSCRIPTOMICS</div>'
                '<h1 style="font-size:32px;margin:10px 0">Visium HD Viewer</h1>'
                '<p style="color:#b7cbd1;font-size:15px">Start with your data. '
                'Bring spatial expression into focus.</p>'),
            card(
                Div(text='<h2 style="margin:0;font-size:19px">Choose your data source</h2>'),
                self.data_source, self.upload_panel, self.example_panel, self.url_panel,
            ),
            self.preprocessing_panel,
            self.load_button, self.status_text, self.progress_panel,
            sizing_mode="stretch_width", spacing=20,
            styles=dict(color="#e5edf0"), name="startup",
        )
        for model in self.layout.references():
            if isinstance(model, (FileInput, Select, TextInput, Checkbox, Button,
                                  RadioButtonGroup, Div)):
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
        self.load_button.disabled = mode == "url"
        self.status_text.text = {
            "upload": "Choose your files to get started.",
            "example": "Choose an example and review preprocessing options.",
            "url": "URL loading is coming soon.",
        }[mode]

    def _add_callbacks(self):
        self.data_source.on_change(
            "active", lambda attr, old, new: self._update_source_mode()
        )

        self.pseudospots_checkbox.on_change("active", lambda attr, old, new: self._reset_progress())
        self.cells_checkbox.on_change("active", lambda attr, old, new: self._reset_progress())

        # La sección Cells de index.html empieza oculta: mostrarla solo si se
        # han segmentado células.
        self.status_text.js_on_change("text", CustomJS(
            args=dict(cells_checkbox=self.cells_checkbox),
            code="""
            if (cb_obj.text === "Done") {
                showMainApp();
                const cells_section = document.getElementById('cells-section');
                if (cells_section && cells_checkbox.active) {
                    cells_section.style.display = 'block';
                }
            }
        """))

    def add_to_document(self):
        self.doc.add_root(self.layout)

    @property
    def _progress_divs(self):
        return (self.tiles_progress_bar, self.bins_progress,
                self.pseudospots_progress, self.cells_progress)

    def _progress_html(self, div, pct=0, status="Waiting"):
        label = {
            self.tiles_progress_bar: "Image tiles",
            self.bins_progress: "Bins",
            self.pseudospots_progress: "Pseudospots",
            self.cells_progress: "Cell segmentation",
        }[div]
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
        for div, enabled in ((self.tiles_progress_bar, True),
                             (self.bins_progress, True),
                             (self.pseudospots_progress, self.pseudospots_checkbox.active),
                             (self.cells_progress, self.cells_checkbox.active)):
            div.text = self._progress_html(div, status="Waiting" if enabled else "Skipped")

    def update_progress(self, div, pct=None, status="Processing"):
        # Workers publish only the latest target; do not queue one UI patch per event.
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

    def _task_progress(self, task):
        return {"tiles": self.tiles_progress_bar, "bins": self.bins_progress,
                "spots": self.pseudospots_progress,
                "cells": self.cells_progress}[task]

    TASK_LABELS = {"tiles": "Image tiles", "bins": "Bins",
                   "spots": "Pseudospots", "cells": "Cell segmentation"}

    def _report_failure(self, task, exc=None):
        self.failed = True
        name = self.TASK_LABELS.get(task, "Loading")
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
        message = f"{name} failed: {reason}" if reason else "Loading failed."
        self.doc.add_next_tick_callback(lambda: setattr(
            self.status_text, "text", escape(message) + " Reload to try again."
        ))

    def _save_uploaded_file(self, file_input):
        self.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
        path = self.UPLOAD_DIR / file_input.filename
        with open(path, "wb") as f:
            f.write(base64.b64decode(file_input.value))
        return path

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
    
    def _check_all_done(self):
        if not self.preparing and not self.failed and all(self.tasks.values()):
            self.doc.add_next_tick_callback(self._all_tasks_finished)
    
    def run_in_thread(self, func, sample_type=None):
        def finish(error=None):
            if error is not None:
                self._report_failure(sample_type, error)
            elif sample_type in self.tasks:
                self.update_progress(self._task_progress(sample_type), 100, "Complete")
            if sample_type in self.tasks:
                self.tasks[sample_type] = True
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

            # The task is only complete once the sample is in the viewer: adding
            # it runs on the document thread and can fail too
            def update_doc():
                try:
                    self.app.add_sample(*output, b2s=self.b2s, sample_type=sample_type)
                except Exception as exc:
                    finish(exc)
                else:
                    finish()
            self.doc.add_next_tick_callback(update_doc)
        
        threading.Thread(target=wrapper).start()

    def generate_tiles(self, base_dir):
        from src.image import ImageProcessor

        if self.dataset:
            if self.dataset.img is None:
                self.dataset.read_img()

            self.image_processor = ImageProcessor(
                self.dataset.img,
                base_dir / "tiles",
                base_dir / "StarDist",
                # pixels_per_bin is the full-res width of one bin of the loaded size
                um_per_px = float(self.dataset.bins_size) / self.dataset.pixels_per_bin
            )
        
            self.doc.add_next_tick_callback(
                lambda: self.app.add_image_component(self.image_processor))

            def worker():
                def on_progress(done, total):
                    pct = int(100 * done / total) if total else 0
                    self.update_progress(self.tiles_progress_bar, pct, "Creating tiles")

                self.image_processor.create_tiles(overwrite=False, progress_callback=on_progress)

            if self.image_processor.mem > self.TILE_MEM_THRESHOLD:
                self.run_in_thread(worker, sample_type="tiles")
            else:
                from PIL import Image
                import numpy as np

                im = Image.fromarray(np.flipud(self.dataset.img))
                self.image_processor.tiles_path.mkdir(parents=True, exist_ok=True)
                im.save(self.image_processor.tiles_path / "z0_x0_y0.png")
                
                self.update_progress(self.tiles_progress_bar, 100, "Complete")
                self.tasks["tiles"] = True
                self._check_all_done()
                
        else:
            raise ValueError

    def build_b2s(self, adata, df_tissue_positions, bins_size):
        from src.utils.bins2spots import Bins2Spots

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
            tile_scale_xfactor=self.image_processor.scale_xfactor,
            tile_scale_yfactor=self.image_processor.scale_yfactor,
            progress_callback=on_progress
        )
        self.update_progress(self.bins_progress, 50, "Building bin geometry")
        bin_squares = self.b2s.bin_squares
        self.update_progress(self.bins_progress, 75, "Adding bins to the viewer")

        def add_bins():
            try:
                self.app.add_sample(bin_squares, self.b2s.adata_hd, b2s=self.b2s, sample_type="bins")
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
                self.dataset.load()
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

    def generate_cells(self):
        from src.utils.operations import CellsPooler

        if self.image_processor is None:
            raise RuntimeError("The tissue image must be loaded before segmenting cells")

        self.update_progress(self.cells_progress, status="Preparing expression data")
        self.ensure_b2s()

        def worker():
            def on_progress(done, total):
                # Four stages: preparation, inference, boundaries, aggregation.
                pct = 50 + int(25 * done / total) if total else 50
                self.update_progress(self.cells_progress, pct, "Extracting cell boundaries")

            self.update_progress(self.cells_progress, 25, "Segmenting tissue")
            polys = self.image_processor.segment_image(min_overlap=320)[1]
            self.update_progress(self.cells_progress, 50, "Extracting cell boundaries")
            cells_gdf = self.image_processor.get_cells_gdf(polys, progress_callback=on_progress)

            cells_pooler = CellsPooler(
                self.b2s.bin_squares,
                self.b2s.adata_hd,
                cells_gdf,
                keep_cols=[]
            )

            self.update_progress(self.cells_progress, 75, "Aggregating cell expression")
            cells_adata = cells_pooler.run(predicate="intersects")

            return cells_gdf, cells_adata

        self.run_in_thread(worker, sample_type="cells")

    def set_tasks(self):
        self.tasks["tiles"] = False
        self.tasks["bins"] = False
        
        if self.pseudospots_checkbox.active:
            self.tasks["spots"] = False
        else:
            self.tasks["spots"] = True
        
        if self.cells_checkbox.active:
            self.tasks["cells"] = False
        else:
            self.tasks["cells"] = True

    def load_data_callback(self):
        if self.loading or self.source_mode == "url":
            return
        if self.source_mode == "upload" and not self._validate_inputs():
            return
        if self.source_mode == "example" and self.example_select.value not in AVAILABLE_DATASETS:
            self.status_text.text = "Please choose an example dataset"
            return

        self.restart()
        self.set_tasks()
        self.status_text.text = "Loading your dataset…"
        self.loading = True
        self.preparing = True
        self._reset_progress()
        self._progress_timer = self.doc.add_periodic_callback(self._render_progress, 50)
        for widget in (self.data_source, self.load_button, self.pseudospots_checkbox,
                       self.cells_checkbox, self.example_select, self.h5_input,
                       self.positions_input, self.image_input, self.resolution_select):
            widget.disabled = True
        self.update_progress(self.tiles_progress_bar, status="Reading image and preparing tiles")

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

    def _prepare_data(self):
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

        self.dataset = DataSet(
            img_path=img_path,
            adata_path=adata_path,
            positions_path=positions_path,
            bins_size=bins_size
        )

        self.generate_tiles(
            self.UPLOAD_DIR if self.source_mode == "upload"
            else Path(AVAILABLE_DATASETS[self.example_select.value]["image_path"]).resolve().parent
        )
        
        if self.pseudospots_checkbox.active:
            self.generate_pseudospots()
        
        if self.cells_checkbox.active:
            self.generate_cells()

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
        self.finished = True
        self.app.build()
        self.status_text.text = "Done"


def create_startup(doc):
    controller = StartupController(doc)
    controller.add_to_document()
