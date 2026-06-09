
import os

from bokeh.models import Range1d

# from .upload import UploadPanel
from .thumbnail import ThumbnailPanel
from .figure import FigurePanel
from .metadata import MetadataPanel
from .embedding import EmbeddingPanel

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

        self.thumbnail_panel = ThumbnailPanel()
        self.figure_panel = FigurePanel()
        self.metadata_panel = MetadataPanel()
        self.embedding_panel = EmbeddingPanel()

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
            self.figure_panel.create_figure(self.image_component.image_tiles.max_dim)

        if not self.figure_panel.image_attached:
            self.figure_panel.attach_image(self.image_component)

        if not self.figure_panel.alpha_setup:
            self.figure_panel.setup_image_alpha(self.image_component)

        if sample_type == "spots":
            sample = Spots(polygons, adata)
            radius = b2s.converter.get_pixels(b2s.converter.VisiumV2.SPOT_DIAMETER) / 2

            geom = Circles(
                radius, 
                self.figure_panel.p, 
                sample, 
                hover_color='green', 
                alias='Spot',
            )
            geom.sample_type = sample_type
            
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

            # width = sample.pixels_per_bin * b2s.tile_scale_xfactor
            # height = sample.pixels_per_bin * b2s.tile_scale_yfactor
            angle = b2s.get_rotation(rad=True)

            geom = Squares(
                angle, 
                self.figure_panel.p, 
                sample, 
                lazy_loading=True, 
                alias='Bin', 
                hover_color='green'
            )
            geom.sample_type = sample_type

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
            geom = Patches(
                self.figure_panel.p, 
                sample, 
                lazy_loading=True, 
                alias='Cell', 
                hover_color='blue'
            )
            geom.sample_type = sample_type

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
        geom._metadata_panel = self.metadata_panel #referencia metadata panel

        # Reconstruir metadata si la app ya está construida
        if hasattr(self, '_app_built') and self._app_built:
            self.doc.add_next_tick_callback(
                lambda: self.metadata_panel.populate(self.geometries)
            )
        self.embedding_panel.set_components(self.geometries)

        # If the image is lowres, adjust the width and height to the max of geometries
        if self.image_component:
            if not self.image_component.are_tiles_needed:
                for r in self.figure_panel.p.renderers:
                    if hasattr(r, "glyph") and hasattr(r.glyph, "url"):
                        x = min([geom.sample.gdf.total_bounds[0] for geom in self.geometries])
                        y = min([geom.sample.gdf.total_bounds[1] for geom in self.geometries])
                        w = max([geom.sample.gdf.total_bounds[2] - geom.sample.gdf.total_bounds[0] for geom in self.geometries])
                        h = max([geom.sample.gdf.total_bounds[3] - geom.sample.gdf.total_bounds[1] for geom in self.geometries])

                        r.glyph.x = x
                        r.glyph.y = y
                        r.glyph.w = w
                        r.glyph.h = h

                        self.figure_panel.p.x_range = Range1d(x, x+w)
                        self.figure_panel.p.y_range = Range1d(y, y+h)

    def add_root(self, layout):
        placeholder = self.doc.select_one({"name": layout.name})
        if hasattr(placeholder, "children"):
            if hasattr(layout, "children"):
                placeholder.children = layout.children
            else:
                placeholder.children = [layout]

    def build(self):
        if self.image_component is None:
            raise ValueError("The ImageComponent must be loaded first with the 'add_image_component' function.")

        self.figure_panel.populate(self.image_component)

        for geom in self.geometries:
            geom._figure_de_div = self.figure_panel.de_div

        self.add_root(self.figure_panel.layout)

        self.thumbnail_panel.populate(self.image_component, self.figure_panel.p)
        self.add_root(self.thumbnail_panel.layout)

        self.metadata_panel.populate(self.geometries)
        self.add_root(self.metadata_panel.layout)

        self.add_root(self.embedding_panel.layout)

        for geom in self.geometries:
            self.add_root(geom.get_widgets_layout(self.doc))
        self._app_built = True  # ← añadir esta línea al final
