
from bokeh.layouts import column
from bokeh.plotting import figure
from bokeh.models import Div

import os

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))


def update_component(doc, comp):
    placeholder = doc.select_one({"name": comp.name})
    if placeholder:
        doc.remove_root(placeholder)
    doc.add_root(comp)

def build_main_app(
        doc, 
        loaded_data,
        image_processor,
        b2s,
        cells_adata=None,
        cells_gdf=None,
        initial_color="#6baed6"
    ):

    # Figure
    p = figure(
        sizing_mode='stretch_width', 
        aspect_ratio=1, 
        x_range=(0, image_processor.max_dim), 
        y_range=(0, image_processor.max_dim),
        background_fill_color='#7e7e7e',
        tools='pan,wheel_zoom,save,reset')
    p.grid.visible = False
    p.xaxis.visible = False
    p.yaxis.visible = False
    p.toolbar.logo = None
    p.toolbar.autohide = True

    ## Image
    from ..components.image import ImageComponent
    image_component = ImageComponent(image_processor)
    image_component.add_tile_source(p)
    image_component.add_img_overlay(p)
    img_alpha_slider = image_component.get_alpha_slider()

    thumb_p = image_component.generate_thumbnail()
    image_component.add_lowres_image(thumb_p)
    image_component.add_viewport(thumb_p)
    image_component.add_viewport_callback(p)

    ## Spots
    if loaded_data["pseudospots"]:
        if b2s.spot_polygons is not None and b2s.adata_v2:
            from src.sample import Spots
            from ..components.glyphs import Circles

            spots = Spots(
                b2s.spot_polygons, 
                b2s.adata_v2)
            radius = b2s.converter.get_pixels(b2s.converter.VisiumV2.SPOT_DIAMETER) / 2

            spots_component = Circles(radius, p, spots, hover_color='green')
            spots_component.add_parameter(
                attr='alpha',
                static=False,
                min_val=0.0,
                max_val=1.0,
                colname=spots_component.get_colnames()[0],
                static_not_sequenced=True,
                not_sequenced_val=0.5,
                colname_hovered=True
            )
            spots_component.add_parameter(
                attr='color',
                static=True,
                val=initial_color,
                not_sequenced_val='#b8b8b8'
            )

            spots_component.manager.fetch_data()

    ## Bins
    if loaded_data["bins"]:
        from src.sample import Bins
        from ..components.glyphs import Squares

        bins = Bins(
            b2s.bin_squares, 
            b2s.adata_hd, 
            on='index', 
            pixels_per_bin=b2s.converter.visiumHD.pixels_per_bin
        )

        width = bins.pixels_per_bin * b2s.tile_scale_xfactor
        height = bins.pixels_per_bin * b2s.tile_scale_yfactor
        angle = b2s.get_rotation(rad=True)

        bins_component = Squares(width, height, angle, p, bins, lazy_loading=True, alias='Bin', hover_color='green')

        bins_component.add_parameter(
            attr='alpha',
            static=False,
            min_val=0.0,
            max_val=1.0,
            colname=bins_component.get_colnames()[0],
            static_not_sequenced=True,
            not_sequenced_val=0.5,
            colname_hovered=True
        )
        bins_component.add_parameter(
            attr='color',
            static=True,
            val=initial_color,
            not_sequenced_val='#b8b8b8'
        )

        bins_component.attach_loader(doc)

    ## Cells
    if loaded_data["cells"]:
        if cells_adata and cells_gdf is not None:
            from src.sample import Cells
            from ..components.glyphs import Patches

            cells = Cells(cells_gdf, cells_adata)
            cells_component = Patches(p, cells, lazy_loading=True, alias='Cell', hover_color='blue')

            cells_component.add_parameter(
                attr='alpha',
                static=False,
                min_val=0.0,
                max_val=1.0,
                colname=cells_component.get_colnames()[0],
                static_not_sequenced=True,
                not_sequenced_val=0.5,
                colname_hovered=True
            )
            cells_component.add_parameter(
                attr='color',
                static=True,
                val=initial_color,
                not_sequenced_val='#b8b8b8'
            )

            cells_component.attach_loader(doc)


    ## Layout
    def get_metadata(sample, name):
        return f"""
        <div>
            <p><strong>{name}:</strong><br>
            Total {name}: {sample.n_geometries}<br>
            Not Sequenced: {sample.n_not_sequenced}</p>
        </div>
        """
    
    text = ""
    if loaded_data["bins"]:
        text += get_metadata(bins, "Bins")

    if loaded_data["pseudospots"]:
        text += get_metadata(spots, "Pseudospots")

    if loaded_data["cells"]:
        text += get_metadata(cells, "Cells")

    metadata = Div(text=f"<div><h3>Metadata</h3>{text}</div>", name='metadata')

    if loaded_data["bins"]:
        bins_widgets = bins_component.get_widgets_layout(doc)
        bins_params = bins_component.get_table_layout()

    if loaded_data["pseudospots"]:
        spots_widgets = spots_component.get_widgets_layout(doc)
        spots_params = spots_component.get_table_layout()

    if loaded_data["cells"]:
        cells_widgets = cells_component.get_widgets_layout(doc)
        cells_params = cells_component.get_table_layout()

    plot_column = column(
        p, img_alpha_slider,
        name='plot_column',
        )
    plot_column.sizing_mode='stretch_width' # Set separately to avoid assign to children 

    thumbnail_column = column(
        thumb_p,
        name='thumbnail_plot',
        sizing_mode='stretch_both'
    )

    # doc.add_root(metadata)

    # doc.add_root(plot_column)
    # doc.add_root(thumbnail_column)

    # if loaded_data["bins"]:
    #     doc.add_root(bins_widgets)
    #     doc.add_root(bins_params)

    # if loaded_data["pseudospots"]:
    #     doc.add_root(spots_widgets)
    #     doc.add_root(spots_params)

    # if loaded_data["cells"]:
    #     doc.add_root(cells_widgets)
    #     doc.add_root(cells_params)

    update_component(doc, metadata)

    update_component(doc, plot_column)
    update_component(doc, thumbnail_column)

    if loaded_data["bins"]:
        update_component(doc, bins_widgets)
        update_component(doc, bins_params)

    if loaded_data["pseudospots"]:
        update_component(doc, spots_widgets)
        update_component(doc, spots_params)

    if loaded_data["cells"]:
        update_component(doc, cells_widgets)
        update_component(doc, cells_params)