
from ..settings.config import sample_kwargs

from .. import settings

import os

from bokeh.layouts import column
from bokeh.models import Range1d, Tabs, TabPanel

from .layout import (
    ThumbnailPanel,
    FigurePanel,
    MetadataPanel,
    EmbeddingPanel,
    ComparisonPanel
)

from .components.image import ImageComponent

from ..startup.sample import Bins, Spots, Cells
from .components.layer import Squares, Circles, Patches

class App:

    INITIAL_ATTRS = settings.INITIAL_ATTRS
    NOT_SEQUENCED_ATTRS = settings.NOT_SEQUENCED_ATTRS

    def __init__(self, doc, image_processor=None):
        self.doc = doc

        self.geometries = []
        self.deconvolution_reference = None
        self.deconvolution_reference_name = None
        self.deconvolution_label_column = ""

        self.thumbnail_panel = ThumbnailPanel()
        self.figure_panel = FigurePanel()
        self.metadata_panel = MetadataPanel()
        self.embedding_panel = EmbeddingPanel(doc)
        self.comparison_panel = ComparisonPanel(doc)

        self.image_component = None
        if image_processor is not None:
            self.image_component = ImageComponent(image_processor)

    # Outline colours that tell segmentation methods apart on the figure
    CELL_LAYER_COLORS = settings.CELL_LAYER_COLORS

    def reset_geometries(self):
        for geometry in self.geometries:
            geometry.roi_manager.dispose()
        self.geometries.clear()
        self.cells_tabs = None

    def on_data_loaded(self, adata, df_tissue_positions=None):
        pass

    def add_image_component(self, image_processor):
        self.image_component = ImageComponent(image_processor)

    def add_sample(self, polygons, adata, b2s, sample_type="spots", method=None, title=None):
        """For cells, `method`/`title` identify the segmentation method; each
        method becomes its own layer."""
        if sample_type not in ["spots", "bins", "cells"]:
            return
        
        if self.image_component is None:
            raise ValueError("The ImageComponent must be loaded first with the 'add_image_component' function.")
        
        if self.figure_panel.p is None:
            self.figure_panel.create_figure(self.image_component.image_processor.max_dim)

        if not self.figure_panel.image_attached:
            self.figure_panel.attach_image(self.image_component)

        if not self.figure_panel.alpha_setup:
            self.figure_panel.setup_image_alpha(self.image_component)

        spatial_transform = getattr(b2s, "spatial_transform", None)

        if sample_type == "spots":
            radius = (spatial_transform.pitch * b2s.converter.VisiumV2.SPOT_DIAMETER / (2 * b2s.bins_size)
                      if spatial_transform is not None else
                      b2s.converter.get_pixels(b2s.converter.VisiumV2.SPOT_DIAMETER) / 2 * b2s.tile_scale_xfactor)
            sample = Spots(polygons, adata, radius=radius, **sample_kwargs())
            sample.spatial_transform = spatial_transform
            geom = Circles('radius', self.figure_panel.p, sample, hover_color='green', alias='Spot')
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
                **sample_kwargs("bins"),
                pixels_per_bin=b2s.converter.visiumHD.pixels_per_bin,
                spatial_transform=spatial_transform,
                bin_size=b2s.bins_size,
                tile_scale_xfactor=b2s.tile_scale_xfactor,
                tile_scale_yfactor=b2s.tile_scale_yfactor,
            )
            um_per_px = (float(sample.bin_size) / sample.pixels_per_bin
                         if sample.pixels_per_bin and sample.bin_size > 0 else None)
            if spatial_transform is not None:
                if self.figure_panel.measure_tool is not None:
                    self.figure_panel.measure_tool.set_transform(spatial_transform, um_per_px)
                if self.figure_panel.roi_readout is not None:
                    self.figure_panel.roi_readout.set_transform(spatial_transform, um_per_px)
            else:
                if self.figure_panel.measure_tool is not None and sample.native_width > 0:
                    self.figure_panel.measure_tool.set_scale(sample.bin_size / sample.native_width,
                                                            sample.bin_size / sample.native_height)
                if self.figure_panel.roi_readout is not None:
                    self.figure_panel.roi_readout.set_scale(1 / sample.tile_scale_xfactor,
                                                          1 / sample.tile_scale_yfactor, um_per_px)

            geom = Squares(
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
            sample = Cells(polygons, adata, **sample_kwargs("cells"))
            sample.spatial_transform = spatial_transform
            n_cell_layers = sum(g.sample_type == "cells" for g in self.geometries)
            geom = Patches(
                self.figure_panel.p, 
                sample, 
                lazy_loading=True, 
                alias='Cell', 
                hover_color='blue',
                key=f"cells_{method}" if method else None,
                name=f"Cells · {title}" if title else None,
            )
            geom.method_title = title or "Cells"
            cell_color = self.CELL_LAYER_COLORS[n_cell_layers % len(self.CELL_LAYER_COLORS)]
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
                val=cell_color,
                not_sequenced_val=self.NOT_SEQUENCED_ATTRS["color"]
            )
            geom.attach_loader(self.doc)
        self.add_geometry(geom)

    def add_geometry(self, geom):
        if self.figure_panel.roi_readout is not None:
            geom.roi_manager.coordinate_scale = self.figure_panel.roi_readout.scale
        geom.deconvolution_manager.reference = self.deconvolution_reference
        geom.deconvolution_manager.reference_name = self.deconvolution_reference_name
        geom.deconvolution_manager.label_column = self.deconvolution_label_column
        if hasattr(self.figure_panel, "de_layout"):
            geom._figure_de_layout = self.figure_panel.de_layout
        self.geometries.append(geom)
        geom._metadata_panel = self.metadata_panel #referencia metadata panel

        # Reconstruir metadata si la app ya está construida
        if hasattr(self, '_app_built') and self._app_built:
            def _update(g=geom):
                self.add_widgets(g)
                self.metadata_panel.populate(self.geometries)
            self.doc.add_next_tick_callback(_update)
        self.embedding_panel.set_components(self.geometries)
        self.comparison_panel.set_components(self.geometries)

        # The image (tiled or a single PNG) always covers the square 0..max_dim of the
        # figure, like the geometries' coordinates: no need to refit it to them.

    def add_root(self, layout):
        placeholder = self.doc.select_one({"name": layout.name})
        if hasattr(placeholder, "children"):
            if hasattr(layout, "children"):
                placeholder.children = layout.children
            else:
                placeholder.children = [layout]
        return placeholder

    def add_widgets(self, geom):
        layout = geom.get_widgets_layout(self.doc)
        if geom.sample_type != "cells":
            self.add_root(layout)
            return
        # Every segmentation method shares the Cells section: one tab per method
        panel = TabPanel(child=layout, title=geom.method_title)
        if getattr(self, "cells_tabs", None) is None:
            self.cells_tabs = Tabs(tabs=[panel], sizing_mode="stretch_width")
            self.add_root(column(self.cells_tabs, name="cells_widgets"))
        else:
            self.cells_tabs.tabs = [*self.cells_tabs.tabs, panel]

    def build(self):
        if self.image_component is None:
            raise ValueError("The ImageComponent must be loaded first with the 'add_image_component' function.")

        self.figure_panel.populate(self.image_component)

        de_root = self.add_root(self.figure_panel.de_layout)
        if de_root is not None:
            self.figure_panel.de_layout = de_root

        for geom in self.geometries:
            geom._figure_de_layout = self.figure_panel.de_layout

        self.add_root(self.figure_panel.layout)

        self.thumbnail_panel.populate(self.image_component, self.figure_panel.p,
                                      roi_readout=self.figure_panel.roi_readout)
        self.add_root(self.thumbnail_panel.layout)

        self.metadata_panel.populate(self.geometries)
        self.add_root(self.metadata_panel.layout)

        self.add_root(self.embedding_panel.layout)
        self.add_root(self.comparison_panel.layout)

        for geom in self.geometries:
            self.add_widgets(geom)
        self._app_built = True
