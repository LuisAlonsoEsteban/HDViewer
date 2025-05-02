
from bokeh.layouts import row, column
from bokeh.models import Button, WMTSTileSource, ColumnDataSource, Slider, CustomJS, HoverTool, CheckboxGroup, ColorPicker, TextInput, Div, AutocompleteInput
from bokeh.plotting import figure, curdoc

import os


from .loader import SourceManager, StaticSourceManager

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.load_data import DataSet
from src.tiles import ImageTiles
from src.utils.bins2spots import Bins2Spots
from src.sample import Bins, Spots

# Constants
INITIAL_COLOR = '#6baed6'
INITIAL_ALPHA = 0.5

# Preload data: Replace 'path_to_file' with actual file paths
dir_base = '/home/lalonso/PhD/VisHD/data/'
dataset = DataSet(dir_base=dir_base, img_file="Visium_HD_Mouse_Small_Intestine_tissue_image.btf", bins_size=8)
img, adata, df_tissue_positions = dataset.load()
print("Data loaded succesfully")

doc = curdoc()

# Deploy image tiles server
# image_tiles = ImageTiles(img, os.path.join(dir_base + "./tiles/"))
# port = 4000
# image_tiles.deploy_server(port=port)
# print(f"Tiles server deployed on port {port}")

# Calculate max dimension for plot range
max_dim = max(img.shape[:2])
background_color = 'white'
p = figure(x_range=(0, max_dim), y_range=(0, max_dim), background_fill_color=background_color)
p.grid.visible = False
# p.xaxis.visible = False
# p.yaxis.visible = False
# p.x_range.bounds = (0, max_dim)
# p.y_range.bounds = (0, max_dim)

# Configure tile source to point to the local tile server
# custom_tile = WMTSTileSource(
#     url= f"http://localhost:{port}" + "/{Z}/{X}/{Y}",  # Ensure port matches deployed tile server
#     tile_size=256,
#     max_zoom=image_tiles.max_zoom-1,
#     min_zoom=0,
#     wrap_around=False,
#     x_origin_offset=0,
#     y_origin_offset=0,
#     initial_resolution=max_dim / image_tiles.tile_size
# )
# p.add_tile(custom_tile)

# # Widgets for image
# img_overlay_source = ColumnDataSource({
#     'x': [max_dim / 2],
#     'y': [max_dim / 2],
#     'width': [max_dim],
#     'height': [max_dim],
#     'alpha': [0.0]  # Initial transparent state
# })
# overlay = p.rect('x', 'y', 'width', 'height', source=img_overlay_source, color=background_color, alpha='alpha')

# # Create a slider to control the alpha of the image overlay
# img_alpha_slider = Slider(start=0.0, end=1.0, value=1.0, step=0.01, title="Image alpha")

# # JavaScript to update the alpha of the overlay rectangle
# update_img_alpha = CustomJS(args=dict(source=img_overlay_source), code="""
#     const data = source.data;
#     const alpha = cb_obj.value;
#     data['alpha'][0] = 1 - alpha;
#     source.change.emit();
# """)
# img_alpha_slider.js_on_change('value', update_img_alpha)


b2s = Bins2Spots(adata, df_tissue_positions, bins_size=dataset.bins_size)
INITIAL_GENE = adata.var_names.tolist()[2]

## Gene input
gene_pool = adata.var_names.tolist() + adata.obs.columns.tolist()
gene_input = AutocompleteInput(title="Gene/Metadata", value=INITIAL_GENE, completions=gene_pool, min_characters=2, case_sensitive=False, restrict=True)

# Spots
from .parameter import Scaler, StaticParameter

# spots = Spots(b2s.spot_polygons, b2s.convert(normalize=False))
# # spots_manager = SourceManager(spots)
# spots_manager = SourceManager(spots)
# spots_manager.add_parameter(Scaler(0, 1, 
#                                    INITIAL_GENE, 
#                                    attr='alpha', 
#                                    static_not_sequenced=True,
#                                    not_sequenced_val=1.0
#                                 ))
# spots_manager.add_parameter(Scaler(0.5, 2, 
#                                    INITIAL_GENE, 
#                                    attr='line_width', 
#                                    static_not_sequenced=True,
#                                    not_sequenced_val=0
#                                 ))
# spots_manager.add_parameter(StaticParameter(1.0, attr='line_alpha', not_sequenced_val=0.0))
# spots_manager.add_parameter(StaticParameter(INITIAL_COLOR, attr='color'))
# # spots_manager.fetch_data()
# spots_manager.add_data()

# radius = b2s.converter.get_pixels(b2s.converter.VisiumV2.SPOT_DIAMETER) / 2

# circle_glyph = p.circle(
#     x='x', 
#     y='y', 
#     color='color', 
#     radius=radius, 
#     fill_alpha='alpha', 
#     hover_color="pink", 
#     line_width='line_width',
#     source=spots_manager.source
# )

# ## Hover 
# hover = HoverTool(tooltips=[("Spot", "x: @x, y: @y, a: @alpha")], renderers=[circle_glyph])
# p.add_tools(hover)

# ## Alpha slider
# spots_alpha_slider = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title="Spots alpha")

# def update_spots_alpha(attr, old, new):
#     print("Changing alpha")
#     def weight_alpha():
#         spots_manager.update_param('alpha', new_max=new)
#     doc.add_next_tick_callback(weight_alpha)

# spots_alpha_slider.on_change('value', update_spots_alpha)

# ## Color picker
# spots_color_picker = ColorPicker(title="Select Spots Color", color=INITIAL_COLOR)
# def update_spots_color(attr, old, new):
#     print("Changing color")
#     def update_color():
#         spots_manager.update_param('color', new_val=new)
#     doc.add_next_tick_callback(update_color)

# spots_color_picker.on_change('color', update_spots_color)

# ## Visibility checkbox
# checkbox = CheckboxGroup(labels=["Show spots"], active=[0])

# def on_checkbox_change(attr, old, new):
#     visible = 0 in new
#     circle_glyph.visible = visible
#     spots_alpha_slider.disabled = not visible
#     spots_color_picker.disabled = not visible

# checkbox.on_change('active', on_checkbox_change)


# def on_input_change(attr, old, new):
#     print("Changing gene")
#     print(spots_manager.get_param("alpha").colname)
#     def update_gene():
#         spots_manager.update_colname(old, new)
#         spots_manager.add_data()
#     doc.add_next_tick_callback(update_gene)

# gene_input.on_change('value', on_input_change)

# layout = row(p, column(img_alpha_slider, checkbox, spots_alpha_slider, spots_color_picker, gene_input))


## Bins
bins = Bins(b2s.bin_points, adata, on='index')
bins_manager = SourceManager(bins)
bins_manager.add_parameter(Scaler(0, 0.8, 
                                   INITIAL_GENE, 
                                   attr='alpha', 
                                   static_not_sequenced=True,
                                   not_sequenced_val=0.0
                                ))
# bins_manager.add_parameter(StaticParameter(1.0, attr='line_alpha', not_sequenced_val=0.0))
# bins_manager.add_parameter(StaticParameter(INITIAL_COLOR, attr='color'))

pixels_per_bin = b2s.converter.visiumHD.pixels_per_bin
angle = b2s.get_rotation(rad=True)
rect_glyph = p.rect(
    x='x', 
    y='y', 
    width=pixels_per_bin, 
    height=pixels_per_bin, 
    angle=-angle, 
    hover_color="green", 
    fill_alpha="alpha", 
    source=bins_manager.source
)
# rect_glyph.visible = True

def update_visible_squares(attr, old, new):
    left = p.x_range.start
    right = p.x_range.end
    bottom = p.y_range.start
    top = p.y_range.end

    if (right - left < 1000) and (top - bottom < 1000):
        # current_time = time.time()
        print("Visualizing bins")
        bbox = [left - 50, right + 50, bottom - 50, top + 50]
        
        async def add_bbox():
            bins_manager.add_data(bbox)
            # rect_glyph.visible = True

        doc.add_next_tick_callback(add_bbox)

p.x_range.on_change('start', update_visible_squares)
# p.x_range.on_change('end', update_visible_squares)
# p.y_range.on_change('start', update_visible_squares)
# p.y_range.on_change('end', update_visible_squares)

layout = row(p, column(gene_input))


# # Bins
# bins = Bins(b2s.bin_points, adata, on='index')
# # bins.deploy_dynamic_loading_server(port=5002)
# bins_source = ColumnDataSource(data={'x': [], 'y': []})
# pixels_per_bin = b2s.converter.visiumHD.pixels_per_bin
# angle = b2s.get_rotation(rad=True)
# rect_glyph = p.rect(
#     x='x', 
#     y='y', 
#     width=pixels_per_bin, 
#     height=pixels_per_bin, 
#     angle=angle, 
#     hover_color="green", 
#     fill_alpha=INITIAL_ALPHA, 
#     source=bins_source
# )
# rect_glyph.visible = False
# hover = HoverTool(tooltips=[("Spot", "x: @x, y: @y")], renderers=[rect_glyph])
# p.add_tools(hover)

# def update_visible_squares(attr, old, new):
#     left = p.x_range.start
#     right = p.x_range.end
#     bottom = p.y_range.start
#     top = p.y_range.end

#     if (right - left < 1000) and (top - bottom < 1000):
#         bbox = [left - 50, right + 50, bottom - 50, top + 50]
#         new_data = bins.get_repr_dict(colname=None, bbox=bbox)  # Assuming 'bins' can handle this call
        
#         async def set_visible():
#             rect_glyph.visible = True
#             bins_source.data = new_data
        
#         doc.add_next_tick_callback(set_visible)

#     else:
#         async def set_invisible():
#             rect_glyph.visible = False

#         doc.add_next_tick_callback(set_invisible)

# # Attach to range changes
# p.x_range.on_change('start', update_visible_squares)
# p.x_range.on_change('end', update_visible_squares)
# p.y_range.on_change('start', update_visible_squares)
# p.y_range.on_change('end', update_visible_squares)

# ## Alpha slider
# bins_alpha_slider = Slider(start=0.0, end=1.0, value=INITIAL_ALPHA, step=0.01, title="Bins alpha")

# def update_bins_alpha(attr, old, new):
#     rect_glyph.glyph.fill_alpha = new
#     rect_glyph.glyph.line_alpha = new

# bins_alpha_slider.on_change('value', update_bins_alpha)

# ## Color picker
# bins_color_picker = ColorPicker(title="Select Bins Color", color=INITIAL_COLOR)

# def update_bins_color(attr, old, new):
#     rect_glyph.glyph.fill_color = new
#     rect_glyph.glyph.line_color = new

# bins_color_picker.on_change('color', update_bins_color)


# Create layout
# layout = row(p, column(img_alpha_slider, checkbox, spots_alpha_slider, spots_color_picker, bins_alpha_slider, bins_color_picker, gene_input))

doc.add_root(layout)
