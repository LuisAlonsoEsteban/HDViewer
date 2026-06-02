
class LazyLoader:

    STATIC_RANGE = 0.6
    OFFSET = 0.7

    HYSTERESIS_RATIO = 0.15

    assert STATIC_RANGE <= OFFSET, "Offset must be larger than the static range for efficient loading"

    def __init__(self, dynamic_source_manager, p):
        self.manager = dynamic_source_manager

        self.x_range = p.x_range
        self.y_range = p.y_range

        self.fixed_static_bbox = None
        self._last_factor = 1

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

    def _get_hysteresis_bounds(self, factor):
        next_factor = None
        prev_factor = None
        factors = getattr(self.manager.sample, "AGGREGATION_FACTORS", (1,))
        if factor in factors:
            idx = factors.index(factor)
            if idx > 0:
                prev_factor = factors[idx - 1]
            if idx < len(factors) - 1:
                next_factor = factors[idx + 1]
        return prev_factor, next_factor

    def get_resolution_level(self):
        sample = self.manager.sample
        if not hasattr(sample, "pixels_per_bin") or not sample.pixels_per_bin:
            return None

        max_span_in_bins = max(self.window_size) / sample.pixels_per_bin
        factors = tuple(getattr(sample, "AGGREGATION_FACTORS", (1,)))
        factor = factors[0]

        for candidate in factors:
            if max_span_in_bins / candidate <= 180:
                factor = candidate
                break
        else:
            factor = factors[-1]

        prev_factor, next_factor = self._get_hysteresis_bounds(self._last_factor)
        if factor != self._last_factor:
            if factor > self._last_factor and next_factor is not None:
                promote_at = 180 * next_factor * (1 + self.HYSTERESIS_RATIO)
                if max_span_in_bins < promote_at:
                    factor = self._last_factor
            elif factor < self._last_factor and prev_factor is not None:
                demote_at = 180 * self._last_factor * (1 - self.HYSTERESIS_RATIO)
                if max_span_in_bins > demote_at:
                    factor = self._last_factor

        self._last_factor = factor
        return sample.get_resolution_level_for_factor(factor)

    def add_range_callback(self, component, doc):
        
        def callback(attr, old, new):

            def update_bbox():
                resolution_level = self.get_resolution_level()
                current_fetch_bbox = self.compute_bbox
                if self.manager.bbox is None:
                    self.manager.set_fetch_context(
                        bbox=current_fetch_bbox,
                        resolution_level=resolution_level,
                    )
                    self.fixed_static_bbox = self.static_bbox
                    self.manager.fetch_data()
                else:
                    x, y = self.x, self.y
                    context_changed = self.manager.set_fetch_context(
                        bbox=self.manager.bbox,
                        resolution_level=resolution_level,
                    )
                    if context_changed or not LazyLoader.is_inside(x, y, self.fixed_static_bbox):
                        self.manager.bbox = current_fetch_bbox
                        self.fixed_static_bbox = self.static_bbox
                        self.manager.fetch_data()

                    if not LazyLoader.is_inside(x, y, self.fixed_static_bbox):
                        self.manager.remove_further_data(x, y, distance=1000)

            doc.add_next_tick_callback(update_bbox)

        self.x_range.on_change('start', callback)
        self.x_range.on_change('end', callback)
        self.y_range.on_change('start', callback)
        self.y_range.on_change('end', callback)
