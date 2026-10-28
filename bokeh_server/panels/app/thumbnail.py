
from . import Panel

class ThumbnailPanel(Panel):

    thumb_p = None

    def __init__(self):
        super().__init__(name="thumbnail_plot")

    def populate(self, image_component, p, roi_readout=None):
        self.thumb_p = image_component.generate_thumbnail()

        image_component.add_lowres_image(self.thumb_p)
        image_component.add_viewport(self.thumb_p)
        image_component.add_viewport_callback(p)
        
        if roi_readout is not None and roi_readout.panel is None:
            roi_readout.place_on(self.thumb_p)

        self.set_children([self.thumb_p])