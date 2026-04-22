
from . import Panel

from bokeh. plotting import figure

class FigurePanel(Panel):

    p, img_alpha_slider = None, None

    def __init__(self):
        super().__init__(name="plot_column")

        self.image_attached = False
        self.alpha_setup = False
    
    def create_figure(self, max_dim):
        self.p = figure(
            sizing_mode='stretch_width', 
            aspect_ratio=1, 
            x_range=(0, max_dim), 
            y_range=(0, max_dim),
            background_fill_color='#7e7e7e',
            tools='pan,wheel_zoom,save,reset')
        
        self.p.grid.visible = False
        self.p.xaxis.visible = False
        self.p.yaxis.visible = False
        self.p.toolbar.logo = None
        self.p.toolbar.autohide = True

    def attach_image(self, image_component):
        if self.p:
            image_component.attach_image(self.p)
            self.image_attached = True

    def setup_image_alpha(self, image_component):
        if self.p:
            image_component.add_img_overlay(self.p)
            self.img_alpha_slider = image_component.get_alpha_slider()
            self.alpha_setup = True
    
    def populate(self, image_component=None):
        if image_component is not None:
            if self.p is None:
                self.create_figure(image_component.image_tiles.max_dim)
            self.set_children([self.p])

            if not self.image_attached:
                self.attach_image(image_component)
            if not self.alpha_setup:
                self.setup_image_alpha(image_component)

            self.add_children(self.img_alpha_slider)



        