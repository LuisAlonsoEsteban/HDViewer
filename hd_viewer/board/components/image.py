
from bokeh.models import (
    WMTSTileSource,
    ColumnDataSource, 
    Slider,
    CustomJS,
    BoxAnnotation
)
from bokeh.plotting import figure

from ... import settings

class ImageComponent:

    TILE_MEM_THRESHOLD = settings.TILE_MEM_THRESHOLD

    def __init__(self, image_processor, port=settings.TILE_SERVER_PORT):
        self.image_processor = image_processor
        self.port = port

        self.server_deployed = False

    @property
    def are_tiles_needed(self):
        return self.image_processor.mem > self.TILE_MEM_THRESHOLD

    def deploy_server(self):
        self.tile_key = self.image_processor.deploy_server(host=settings.SERVER_HOST, port=self.port)
        self.server_deployed = True

    def generate_tile_source(self):
        try:
            self.deploy_server()
        except:
            raise ValueError("An error has occurred deploying the server providing the tiles")

        self.tile_source = WMTSTileSource(
                                url=f"http://localhost:{self.port}/{self.tile_key}" + "/{Z}/{X}/{Y}",
                                tile_size=self.image_processor.tile_size,
                                max_zoom=self.image_processor.max_zoom,
                                min_zoom=0,
                                wrap_around=False,
                                x_origin_offset=0,
                                y_origin_offset=0,
                                initial_resolution=self.image_processor.max_dim / self.image_processor.tile_size
                            )
        
        return self.tile_source
    
    def add_tile_source(self, p):
        tile_source = self.generate_tile_source()
        if tile_source is not None:
            p.add_tile(tile_source)

    def attach_image(self, p, w=None, h=None):
        if self.are_tiles_needed:
            self.add_tile_source(p)
        else:
            self.deploy_server()
            self.add_lowres_image(p, w=w, h=h)

    def get_overlay_source(self):
        self.overlay_source = ColumnDataSource({
                    'x': [self.image_processor.max_dim / 2],
                    'y': [self.image_processor.max_dim / 2],
                    'width': [self.image_processor.max_dim],
                    'height': [self.image_processor.max_dim],
                    'alpha': [0.0]  # Initial transparent state
                })
        return self.overlay_source
    
    def add_img_overlay(self, p):
        return p.rect(
                    'x', 
                    'y', 
                    'width', 
                    'height', 
                    source=self.get_overlay_source(), 
                    color=p.background_fill_color, 
                    alpha='alpha')
    
    def get_alpha_slider(self, width=300, margin=(10,50,10,50)):
        self.alpha_slider = Slider(start=0.0, end=1.0, value=1.0, step=0.01, title="Image alpha", width=width, margin=margin)
        self.update_img_alpha = CustomJS(args=dict(source=self.overlay_source), code="""
                                        const data = source.data;
                                        const alpha = cb_obj.value;
                                        data['alpha'][0] = 1 - alpha;
                                        source.change.emit();
                                    """)
        self.alpha_slider.js_on_change('value', self.update_img_alpha)
        return self.alpha_slider

    def generate_thumbnail(self):
        thumb_p = figure(
            sizing_mode='scale_width',
            tooltips="Thumbnail overview", 
            x_range=(0, self.image_processor.max_dim), 
            y_range=(0, self.image_processor.max_dim),
            match_aspect=True
        )

        thumb_p.toolbar.active_drag = None
        thumb_p.toolbar.active_scroll = None
        thumb_p.toolbar.active_tap = None
        
        thumb_p.grid.visible = False
        thumb_p.xaxis.visible = False
        thumb_p.yaxis.visible = False

        thumb_p.toolbar.logo = None
        thumb_p.toolbar_location = None

        return thumb_p
    
    def add_lowres_image(self, p, w=None, h=None):
        if self.server_deployed:
            p.image_url(
                url=[f"http://localhost:{self.port}/{self.tile_key}/0/0/0"], 
                x=0, 
                y=h or self.image_processor.max_dim, 
                w=w or self.image_processor.max_dim, 
                h=h or self.image_processor.max_dim)
        else:
            pass
    
    def add_viewport(self, thumb_p, fill_color='black', line_color='black'):
        self.viewport = BoxAnnotation(left=0, right=0, top=0, bottom=0, fill_alpha=0.1, fill_color=fill_color, line_width=2, line_color=line_color)
        thumb_p.add_layout(self.viewport)

    def add_viewport_callback(self, p):
        def update_viewport(attr, old, new):
            self.viewport.left = p.x_range.start
            self.viewport.right = p.x_range.end
            self.viewport.bottom = p.y_range.start
            self.viewport.top = p.y_range.end

        p.x_range.on_change('start', update_viewport)
        p.x_range.on_change('end', update_viewport)
        p.y_range.on_change('start', update_viewport)
        p.y_range.on_change('end', update_viewport)
