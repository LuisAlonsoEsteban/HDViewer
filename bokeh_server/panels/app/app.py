
import os

# from .upload import UploadPanel
from .thumbnail import ThumbnailPanel
from .figure import FigurePanel
from .metadata import MetadataPanel

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

from ...components.image import ImageComponent

from src.sample import Bins, Spots, Cells
from ...components.glyphs import Squares, Circles, Patches

class App:

    INITIAL_ATTRS = {
        "color": "#6baed6",
        "alpha": 0.5,
    }
    NOT_SEQUENCED_ATTRS = {
        "color": "#b8b8b8",
        "alpha": 0.5,
    }

    def __init__(self, doc, image_processor=None):
        self.doc = doc

        self.geometries = []

        # self.upload_panel = UploadPanel(self.on_data_loaded)
        self.thumbnail_panel = ThumbnailPanel()
        self.figure_panel = FigurePanel()
        self.metadata_panel = MetadataPanel()

        self.image_component = None
        if image_processor is not None:
            self.image_component = ImageComponent(image_processor)

    def reset_geometries(self):
        self.geometries.clear()

    def on_data_loaded(self, adata, df_tissue_positions=None):
        pass

    def add_image_component(self, image_processor):
        self.image_component = ImageComponent(image_processor)

    def add_sample(self, polygons, adata, b2s=None, sample_type="spots"):
        if sample_type not in ["spots", "bins", "cells"]:
            return
        
        if self.image_component is None:
            raise ValueError("The ImageComponent must be loaded first with the 'add_image_component' function.")
        
        if self.figure_panel.p is None:
            self.figure_panel.populate(self.image_component)

        if sample_type == "spots":
            sample = Spots(polygons, adata)
            radius = b2s.converter.get_pixels(b2s.converter.VisiumV2.SPOT_DIAMETER) / 2

            geom = Circles(radius, self.figure_panel.p, sample, hover_color='green')
            geom.add_parameter(
                attr='alpha',
                static=False,
                min_val=0.0,
                max_val=1.0,
                colname=geom.get_colnames()[0],
                static_not_sequenced=True,
                not_sequenced_val=self.NOT_SEQUENCED_ATTRS["alpha"],
                colname_hovered=True
            )
            geom.add_parameter(
                attr='color',
                static=True,
                val=self.INITIAL_ATTRS["color"],
                not_sequenced_val=self.NOT_SEQUENCED_ATTRS["color"]
            )

            geom.manager.fetch_data()

        elif sample_type == "bins":
            sample = Bins(
                polygons, 
                adata, 
                on='index', 
                pixels_per_bin=b2s.converter.visiumHD.pixels_per_bin
            )

            width = sample.pixels_per_bin * b2s.tile_scale_xfactor
            height = sample.pixels_per_bin * b2s.tile_scale_yfactor
            angle = b2s.get_rotation(rad=True)

            geom = Squares(width, height, angle, self.figure_panel.p, sample, lazy_loading=True, alias='Bin', hover_color='green')

            geom.add_parameter(
                attr='alpha',
                static=False,
                min_val=0.0,
                max_val=1.0,
                colname=geom.get_colnames()[0],
                static_not_sequenced=True,
                not_sequenced_val=self.NOT_SEQUENCED_ATTRS["alpha"],
                colname_hovered=True
            )
            geom.add_parameter(
                attr='color',
                static=True,
                val=self.INITIAL_ATTRS["color"],
                not_sequenced_val=self.NOT_SEQUENCED_ATTRS["color"]
            )

            geom.attach_loader(self.doc)

        elif sample_type == "cells":
            sample = Cells(polygons, adata)
            geom = Patches(self.figure_panel.p, sample, lazy_loading=True, alias='Cell', hover_color='blue')

            geom.add_parameter(
                attr='alpha',
                static=False,
                min_val=0.0,
                max_val=1.0,
                colname=geom.get_colnames()[0],
                static_not_sequenced=True,
                not_sequenced_val=0.5,
                colname_hovered=True
            )
            geom.add_parameter(
                attr='color',
                static=True,
                val=self.INITIAL_ATTRS["color"],
                not_sequenced_val=self.NOT_SEQUENCED_ATTRS["color"]
            )

            geom.attach_loader(self.doc)

        self.add_geometry(geom)

    def add_geometry(self, geom):
        self.geometries.append(geom)

    def add_root(self, layout):
        # for a_root in self.doc.select(dict(name=layout.name)):
        #     self.doc.remove_root(a_root)
        # self.doc.add_rot(layout)

        placeholder = self.doc.select_one({"name": layout.name})
        if hasattr(placeholder, "children"):
            placeholder.children = [layout]

    def build(self):
        if self.image_component is None:
            raise ValueError("The ImageComponent must be loaded first with the 'add_image_component' function.")
        
        # self.upload_panel.populate()
        # doc.add_root(self.upload_panel.layout)

        self.figure_panel.populate(self.image_component)
        self.add_root(self.figure_panel.layout)

        self.figure_panel.p.circle([0, 100], [0, 100], size=20, color="red")

        self.thumbnail_panel.populate(self.image_component, self.figure_panel.p)
        self.add_root(self.thumbnail_panel.layout)

        self.metadata_panel.populate(self.geometries)
        self.add_root(self.metadata_panel.layout)

        for geom in self.geometries:
            self.add_root(geom.get_widgets_layout(self.doc))
            self.add_root(geom.get_table_layout())

