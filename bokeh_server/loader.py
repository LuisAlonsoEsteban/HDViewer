
from time import sleep

class LazyLoader:

    STATIC_RANGE = 0.6
    OFFSET = 0.7

    assert STATIC_RANGE <= OFFSET, "Offset must be larger than the static range for efficient loading"

    def __init__(self, dynamic_source_manager, p):
        self.manager = dynamic_source_manager

        self.x_range = p.x_range
        self.y_range = p.y_range

        self.fixed_static_bbox = None

    @property
    def window_size(self):
        return self.x_range.end - self.x_range.start, self.y_range.end - self.y_range.start
    
    @property
    def npixels(self):
        return self.window_size[0] * self.window_size[1]
    
    @property
    def x(self):
        return (self.x_range.start + self.x_range.end) / 2
    
    @property
    def y(self):
        return (self.y_range.start + self.y_range.end) / 2
    
    @property
    def static_bbox(self):
        return [
            self.x - self.window_size[0] * self.STATIC_RANGE / 2,
            self.x + self.window_size[0] * self.STATIC_RANGE / 2,
            self.y - self.window_size[1] * self.STATIC_RANGE / 2,
            self.y + self.window_size[1] * self.STATIC_RANGE / 2,
        ]
    
    @property
    def compute_bbox(self):
        return [
            self.x - self.window_size[0] * (1 + self.OFFSET) / 2,
            self.x + self.window_size[0] * (1 + self.OFFSET) / 2,
            self.y - self.window_size[1] * (1 + self.OFFSET) / 2,
            self.y + self.window_size[1] * (1 + self.OFFSET) / 2,
        ]
    
    @staticmethod
    def is_inside(x, y, bbox):
        return x > bbox[0] and x < bbox[1] and y > bbox[2] and y < bbox[3]

    def add_range_callback(self, component, doc):
        
        def callback(attr, old, new):

            def update_bbox():
                if self.npixels < 1e6:
                    if self.manager.bbox is None:
                        self.manager.bbox = self.compute_bbox
                        self.fixed_static_bbox = self.static_bbox
                        self.manager.fetch_data()
                    else:
                        x, y = self.x, self.y
                        if not LazyLoader.is_inside(x, y, self.fixed_static_bbox):
                            self.fixed_static_bbox = self.static_bbox
                            self.manager.bbox = self.compute_bbox
                            self.manager.fetch_data()
                            self.manager.remove_further_data(x, y, distance=1000)
                        
                    if 0 in component.checkbox.active:
                        component.set_visibility(True, disable_checkbox=True)

                else:
                    component.set_visibility(False, disable_checkbox=True)

            doc.add_next_tick_callback(update_bbox)

        self.x_range.on_change('start', callback)
        self.x_range.on_change('end', callback)
        self.y_range.on_change('start', callback)
        self.y_range.on_change('end', callback)