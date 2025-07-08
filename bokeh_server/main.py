
from bokeh. plotting import curdoc

from .panels.upload import UploadPanel
from .panels.thumbnail import ThumbnailPanel
from .panels.figure import FigurePanel

from .components.image import ImageComponent

class VisHDApp:

    def __init__(self):
        self.geometries = []

        self.upload_panel = UploadPanel(self.on_image_loaded, self.on_image_loaded)
        self.thumbnail_panel = ThumbnailPanel()
        self.figure_panel = FigurePanel()


    def on_image_loaded(self, image_tiles):
        print("Here")
        # image_component = ImageComponent(image_tiles)
        # self.figure_panel.populate(image_component)
        # self.thumbnail_panel.populate(self.figure_panel.p)

    def on_data_loaded(self, adata, df_tissue_positions=None):
        if "Visium HD" in self.upload_panel.technology_selector.value:
            pass
        elif "Visium v2" in self.upload_panel.technology_selector.value:
            pass

    def add_geometry(self, geom):
        self.geometries.append(geom)
        doc.add_root(geom.get_widgets_layout(doc))
        doc.add_root(geom.get_table_layout())

    def add_roots(self, doc):
        doc.add_root(self.upload_panel.layout)
        doc.add_root(self.thumbnail_panel.layout)
        doc.add_root(self.figure_panel.layout)

        for geom in self.geometries:
            doc.add_root(geom.get_widgets_layout(doc))
            doc.add_root(geom.get_table_layout())
                    

doc = curdoc()

server = VisHDApp()
server.add_roots(doc)

from bokeh.models import Div
doc.add_root(Div(name="spots_widgets"))
doc.add_root(Div(name="bins_widgets"))
doc.add_root(Div(name="cells_widgets"))
doc.add_root(Div(name="spots_params"))
doc.add_root(Div(name="bins_params"))
doc.add_root(Div(name="cells_params"))