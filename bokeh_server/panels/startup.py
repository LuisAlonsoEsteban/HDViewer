from bokeh.models import (
    FileInput, Button, Div, Select,
    Checkbox, CustomJS,
    ColumnDataSource
)
from bokeh.layouts import row, column

import base64
import os
import sys
from pathlib import Path

import threading
from time import sleep

from .app import build_main_app
from ..datasets import AVAILABLE_DATASETS

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.load_data import DataSet


class StartupController:

    UPLOAD_DIR = Path("./tmp")

    def __init__(self, doc):
        self.doc = doc

        self.dataset = None

        self.image_processor = None
        self.b2s = None

        self.pseudospots_adata = None
        self.cells_adata = None
        self.cells_gdf = None

        self._build_widgets()
        self.add_callbacks()
        self._build_layout()

        self.pending_tasks = 0
        self.lock = threading.Lock()

    def _build_widgets(self):
        self.data_source = Select(
            title="Data Source:",
            options=[
                ("upload", "Upload data"),
                ("example", "Use Example Dataset"),
            ],
            value="upload",
        )

        self.h5_input = FileInput(accept=".h5")
        self.h5_input_widget = column(
            Div(text="<b>1. Expression data (.h5)</b>"),
            self.h5_input
        )

        self.positions_input = FileInput(accept=".parquet")
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

        self.image_input = FileInput(accept=".btf")
        self.tiles_progress_bar = Div(text="", width=1000)
        self.image_input_widget = column(
            Div(text="<b>4. Tissue image (.btf)</b>"),
            self.image_input
        )

        self.pseudospots_checkbox = Checkbox(
            label="Generate Visium v2 pseudospots",
            active=False
        )
        self.cells_checkbox = Checkbox(
            label="Segment cells",
            active=False
        )

        self.example_select = Select(
            title="Choose example dataset",
            options=[(key.replace(" ", "_"), key) for key in AVAILABLE_DATASETS.keys()],
            value=list(AVAILABLE_DATASETS.keys())[0]
        )

        self.pseudospots_progress = Div(text="", width=200)
        self.cells_progress = Div(text="", width=200)
        self.status_text = Div(text="Upload your data")

        self.load_button = Button(label="Load Data", button_type="success")
        self.load_button.on_click(self.load_data_callback)

    def _build_layout(self):
        pseudospots_row = row(
            self.pseudospots_checkbox,
            self.pseudospots_progress
        )

        cells_row = row(
            self.cells_checkbox,
            self.cells_progress
        )

        processing_options = column(
            pseudospots_row,
            cells_row
        )

        self.layout = column(
            Div(text="<h1>Visium HD Viewer</h1>"),
            self.data_source,

            self.h5_input_widget,
            self.positions_input_widget,
            self.resolution_select_widget,
            self.image_input_widget,

            row(self.example_select, self.tiles_progress_bar),

            Div(text="<b>5. Processing options</b>"),
            processing_options,

            self.load_button,
            self.status_text,

            width=500,
            name="startup"
        )

    def add_callbacks(self):
        def _update_source_mode(attr, old, new):
            is_upload = (new == "upload")

            self.h5_input_widget.visible = is_upload
            self.positions_input_widget.visible = is_upload
            self.resolution_select_widget.visible = is_upload
            self.image_input_widget.visible = is_upload

            self.example_select.visible = not is_upload

        self.data_source.on_change("value", _update_source_mode)

        # JS trigger model
        self.js_trigger = ColumnDataSource(data=dict(done=[0]))
        self.js_trigger.js_on_change("data", CustomJS(code="""
            showMainApp();
        """))

    def add_to_document(self):
        self.doc.add_root(self.layout)

    def update_progress(self, div, pct):
        def _update():
            div.text = f"""
            <div style="width:200px; height:15px; background:#555;">
                <div style="width:{pct}%; height:15px; background:#4CAF50;"></div>
            </div>
            """
        self.doc.add_next_tick_callback(_update)

    def _save_uploaded_file(self, file_input):
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
            self.status_text.text = "Please upload a .btf image file"
            return False

        return True
    
    def _task_done(self):
        with self.lock:
            self.pending_tasks -= 1
            done = self.pending_tasks == 0

        if done:
            self.doc.add_next_tick_callback(self._all_tasks_finished)
    
    def run_in_thread(self, func):
        def wrapper():
            try:
                func()
            finally:
                self._task_done()

        with self.lock:
            self.pending_tasks += 1

        threading.Thread(target=wrapper).start()

    def generate_tiles(self, base_dir):
        from src.image import ImageProcessor

        if self.dataset:
            if self.dataset.img is None:
                self.dataset.read_img()

            self.image_processor = ImageProcessor(
                self.dataset.img,
                base_dir / "new_tiles",
                base_dir / "StarDist"
            )

            def worker():
                def on_progress(done, total):
                    pct = int(100 * done / total) if total else 0
                    self.update_progress(self.tiles_progress_bar, pct)

                self.image_processor.create_tiles(progress_callback=on_progress)

            self.run_in_thread(worker)

        else:
            raise ValueError

    def build_b2s(self, adata, df_tissue_positions, bins_size):
        from src.utils.bins2spots import Bins2Spots

        def on_progress(step, done, total):
            pct = int((100 / 2) * (step + done / total)) if total else 0
            self.update_progress(self.pseudospots_progress, pct)

        self.b2s = Bins2Spots(
            adata,
            df_tissue_positions,
            bins_size=bins_size,
            tile_scale_xfactor=self.image_processor.scale_xfactor,
            tile_scale_yfactor=self.image_processor.scale_yfactor,
            progress_callback=on_progress
        )

    def ensure_b2s(self):
        if self.b2s is None:
            if not self.dataset.adata or self.dataset.df_tissue_positions:
                self.dataset.load()
            self.build_b2s(self.dataset.adata, self.dataset.df_tissue_positions, self.dataset.bins_size)

    def generate_pseudospots(self):
        self.ensure_b2s()

        def worker():
            self.b2s.convert(normalize=False)

        self.run_in_thread(worker)

    def generate_cells(self):
        from src.utils.operations import CellsPooler

        if self.image_processor is None:
            return
        
        self.ensure_b2s()

        def worker():
            def on_progress(done, total):
                pct = int(100 * done / total) if total else 0
                self.update_progress(self.cells_progress, pct)

            self.cells_gdf = self.image_processor.get_cells_gdf(
                                    self.image_processor.segment_image()[1], 
                                    progress_callback=on_progress
                                )


            cells_pooler = CellsPooler(
                self.b2s.bin_squares,
                self.b2s.adata_hd,
                self.cells_gdf,
                keep_cols=[]
            )

            self.cells_adata = cells_pooler.run(predicate="intersects")

        self.run_in_thread(worker)

    def load_data_callback(self):
        if self.data_source.value == "upload":
            if not self._validate_inputs():
                return

            self.status_text.text = "Loading..."

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
            self.UPLOAD_DIR if self.data_source.value == "upload"
            else Path(AVAILABLE_DATASETS[self.example_select.value]["image_path"]).resolve().parent
        )

        if self.pseudospots_checkbox.active:
            sleep(1) # Ensure the Bins2Spots has been created
            self.generate_pseudospots()

        if self.cells_checkbox.active:
            sleep(5) # Ensure the ImageProcessor has been created
            self.generate_cells()

    def _all_tasks_finished(self):
        build_main_app(
            self.doc, 
            self.image_processor,
            self.b2s,
            self.cells_adata,
            self.cells_gdf
        )

        self.js_trigger.data = dict(done=[self.js_trigger.data["done"][0] + 1])


def create_startup(doc):
    controller = StartupController(doc)
    controller.add_to_document()
