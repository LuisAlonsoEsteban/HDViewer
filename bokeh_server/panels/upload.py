
from bokeh.models import Div, FileInput, NumericInput, MultiChoice

from . import Panel
from ..uploader import DataInputPaths

import os
import base64
import io

class UploadPanel(Panel):

    def __init__(self, on_data_loaded_callback, on_image_loaded_callback):
        super().__init__(name="metadata")

        self.data_loader = DataInputPaths(on_data_loaded_callback, on_image_loaded_callback)

    def get_placeholder(self):
        self.adata_input = FileInput(title="Select the filtered counts matrix (.h5)", accept=".h5")
        def on_load_adata(attr, old, new):
            self.data_loader.adata = new
        self.adata_input.on_change('value', on_load_adata)

        self.positions_input = FileInput(title="Select the positions parquet file", accept=".parquet")
        def on_load_positions(attr, old, new):
            self.data_loader.positions = new
        self.positions_input.on_change('value', on_load_positions)

        self.bins_size_input = NumericInput(low=2, high=16, title="Enter the micrometers of a bin side")
        def on_load_bins_size(attr, old, new):
            self.data_loader.bins_size = new
        self.bins_size_input.on_change('value', on_load_bins_size)

        self.image_input = FileInput(title="Select the image file", accept="*")
        def on_load_img(attr, old, new):
            print(new)
            print("----------")
            decoded = base64.b64decode(new)
            print(decoded)
            print("----------")
            data = io.BytesIO(decoded)
            
            self.data_loader.img = data

        self.image_input.on_change('value', on_load_img)

        self.technology_selector = MultiChoice(value=["Visium HD"], options=["Visium HD", "Visium v2"])
        def on_change_technology(attr, old, new):
            if new == ["Visium HD"]:
                self.positions_input.disabled = False
                self.bins_size_input.disabled = False
            elif new == ["Visium v2"]:
                self.positions_input.disabled = True
                self.bins_size_input.disabled = True
        self.technology_selector.on_change('value', on_change_technology)

        return [
            self.technology_selector,
            self.adata_input,
            self.positions_input,
            self.bins_size_input,
            self.image_input
        ]