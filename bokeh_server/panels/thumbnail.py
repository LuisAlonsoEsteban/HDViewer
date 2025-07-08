
from . import Panel

class ThumbnailPanel(Panel):

    def __init__(self):
        super().__init__(name="thumbnail_plot")

    def populate(self, image_component, p):
        self.thumb_p = image_component.generate_thumbnail()
        image_component.add_lowres_image(self.thumb_p)
        image_component.add_viewport(self.thumb_p)
        image_component.add_viewport_callback(p)

        self.set_children([self.thumb_p])
