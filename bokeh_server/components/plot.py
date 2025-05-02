
from ..data.tiles import ImageTiles
from bokeh.models import (
    WMTSTileSource,
    ColumnDataSource,
    Slider,
    HoverTool
)

from .widget import Widget

class Component:

    def __init__(self):
        self.layout = None
        self.source = {
            'x': [],
            'y': []
        }
        self.widgets = []

    @property
    def layout(self):
        if self.layout is None:
            self.add_layout()
        return self.layout

    def add_layout(self):
        pass

    def add_hover(self, tooltips, p):
        hover = HoverTool(tooltips=[tooltips], renderers=[self.layout])
        p.add_tools(hover)

    def add_widget(self, widget):
        self.widgets.append(widget)

    def connect_widgets(self):
        for widget in self.widgets:
            widget.connect(self.source)

class ImageComponent:

    def __init__(self, *args):
        self.img = ImageTiles(*args)
        self.max_dim = max(self.img.image.shape[:2])

    def get_tile_source(self):
        url = self.img.deploy_server()
        return WMTSTileSource(
            url = url,
            tile_size = self.img.tile_size,
            min_xoom = 0,
            max_zoom = self.img.max_zoom,
            wrap_around = False,
            x_origin_offset = 0,
            y_origin_offset = 0,
            initial_resolution = self.max_dim / self.img.tile_size
        )

    class Overlay(Component):

        INITIAL_ALPHA = 0.0

        def __init__(self, max_dim):
            super().__init__()

            self.source = ColumnDataSource(data={
                'x': [max_dim / 2],
                'y': [max_dim / 2],
                'width': [max_dim],
                'height': [max_dim],
            })

        def add_layout(self, p):
            if self.layout is None:
                self.layout = p.rect(
                                'x', 
                                'y', 
                                'width', 
                                'height', 
                                source=self.source, 
                                color=p.background_fill_color, 
                                alpha=self.INITIAL_ALPHA
                            )
                
        def add_widgets(self):
            component = Slider(start=0.0, end=1.0, value=1-self.INITIAL_ALPHA, step=0.01, title="Image alpha")
            def update_alpha(attr, old, new):
                self.layout.fill_alpha = 1- new
                self.layout.line_alpha = 1- new
            self.add_widget(Widget(component, update_alpha, "value"))

            self.connect_widgets()
        
    def add_layout(self, p):
        tile_source = p.add_tile(self.get_tile_source())
        overlay = self.Overlay(self.max_dim).add_layout(p)
        return tile_source, overlay
        
class SpotsComponents(Component):

    INITIAL_COLOR = '#6baed6'
    INITIAL_ALPHA = 0.5 

    def __init__(self, spots, radius, hover_color):
        super().__init__()

        self.source = ColumnDataSource(data=spots.get_repr_dict(colname=None))
        self.radius = radius
        self.hover_color = hover_color

    def add_layout(self, p):
        if self.layout is None:
            self.layout = p.circle(
                                x = 'x', 
                                y = 'y', 
                                color = self.INITIAL_COLOR, 
                                radius = self.radius, 
                                fill_alpha = self.INITIAL_ALPHA, 
                                hover_color = self.hover_color, 
                                source = self.source
                            )
            self.add_hover(("Spot", "x: @x, y: @y"),p)
            
class BinsComponent(Component):

    INITIAL_COLOR = '#6baed6'
    INITIAL_ALPHA = 0.5 

    def __init__(self, bins, pixels_per_bin, angle, hover_color):
        super().__init__()

        self.source = ColumnDataSource(data=bins.get_repr_dict(colname=None))
        self.pixels_per_bin = pixels_per_bin
        self.angle = angle
        self.hover_color = hover_color

    def add_layout(self, p):
        if self.layout is None:
            self.layout = p.rect(
                                x='x', 
                                y='y', 
                                width=self.pixels_per_bin, 
                                height=self.pixels_per_bin, 
                                angle=self.angle, 
                                hover_color=self.hover_color, 
                                fill_alpha=self.INITIAL_ALPHA, 
                                source=self.source
                            )
            self.layout.visible = False
            self.add_hover(("Bin", "x: @x, y: @y"), p)
            
    