
from bokeh.layouts import column
from bokeh.plotting import figure, curdoc

import os

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.load_data import DataSet
from src.tiles import ImageTiles
from src.utils.bins2spots import Bins2Spots
from src.sample import Bins, Spots, Cells
from src.utils.operations import CellsPooler

# Constants
INITIAL_COLOR = '#6baed6'

dir_base = '/home/lalonso/PhD/VisHD/data/'
dataset = DataSet(dir_base=dir_base, img_file="Visium_HD_Mouse_Small_Intestine_tissue_image.btf", bins_size=16)
img, adata, df_tissue_positions = dataset.load()
print("Data loaded succesfully")

# Deploy image tiles server
image_tiles = ImageTiles(img, os.path.join(dir_base, "./tiles/"), segm_path=os.path.join(dir_base, "./StarDist/"))

# Figure
p = figure(
    sizing_mode='stretch_width', 
    aspect_ratio=1, 
    x_range=(0, image_tiles.max_dim), 
    y_range=(0, image_tiles.max_dim),
    background_fill_color='#7e7e7e',
    tools='pan,wheel_zoom,save,reset')
p.grid.visible = False
p.xaxis.visible = False
p.yaxis.visible = False
p.toolbar.logo = None
p.toolbar.autohide = True

doc = curdoc()

## Image
from .components.image import ImageComponent
image_component = ImageComponent(image_tiles)
image_component.add_tile_source(p)
image_component.add_img_overlay(p)
img_alpha_slider = image_component.get_alpha_slider()

thumb_p = image_component.generate_thumbnail()
image_component.add_lowres_image(thumb_p)
image_component.add_viewport(thumb_p)
image_component.add_viewport_callback(p)


# Bins2Spots
b2s = Bins2Spots(
    adata, 
    df_tissue_positions, 
    bins_size=dataset.bins_size, 
    tile_scale_xfactor=image_tiles.scale_xfactor,
    tile_scale_yfactor=image_tiles.scale_yfactor
)

## Spots
from .components.glyphs import Circles

spots = Spots(
    b2s.spot_polygons, 
    b2s.convert(normalize=False))
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
    val=INITIAL_COLOR,
    not_sequenced_val='#b8b8b8'
)

spots_component.manager.fetch_data()

## Bins
from .components.glyphs import Squares

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
    val=INITIAL_COLOR,
    not_sequenced_val='#b8b8b8'
)

bins_component.attach_loader(doc)

## Cells
from .components.glyphs import Patches

cells_gdf = image_tiles.get_cells_gdf(image_tiles.segment_image()[1])
cells_pooler = CellsPooler(b2s.bin_squares, b2s.adata_hd, cells_gdf, keep_cols=[])

cells_adata = cells_pooler.run(predicate='intersects')
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
    val=INITIAL_COLOR,
    not_sequenced_val='#b8b8b8'
)

cells_component.attach_loader(doc)


## Layout
metadata = column(
    Div(text=f"""
        <div>
            <h3>Metadata</h3>
            <div>
                <p><strong>Spots:</strong><br>
                Total Spots: {spots.n_geometries}<br>
                Not Sequenced: {spots.n_not_sequenced}</p>
            </div>
            <div>
                <p><strong>Bins:</strong><br>
                Total Bins: {bins.n_geometries}<br>
                Not Sequenced: {bins.n_not_sequenced}</p>
            </div>
            <div>
                <p><strong>Cells:</strong><br>
                Total Cells: {cells.n_geometries}<br>
                Not Sequenced: {cells.n_not_sequenced}</p>
            </div>
        </div>
        """),
    name="metadata"
)

spots_widgets = spots_component.get_widgets_layout(doc)
bins_widgets = bins_component.get_widgets_layout(doc)
cells_widgets = cells_component.get_widgets_layout(doc)

spots_params = spots_component.get_table_layout()
bins_params = bins_component.get_table_layout()
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

doc.add_root(metadata)

doc.add_root(spots_widgets)
doc.add_root(bins_widgets)
doc.add_root(cells_widgets)

doc.add_root(plot_column)

doc.add_root(thumbnail_column)

doc.add_root(spots_params)
doc.add_root(bins_params)
doc.add_root(cells_params)

