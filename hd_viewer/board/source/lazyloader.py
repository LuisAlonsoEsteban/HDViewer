
import math

from bokeh.core.property.descriptors import UnsetValueError
from bokeh.document import Document

from ... import settings
from ...core import viewport

class LazyLoader:

    STATIC_RANGE = settings.STATIC_RANGE
    OFFSET = settings.OFFSET

    HYSTERESIS_RATIO = settings.HYSTERESIS_RATIO

    assert STATIC_RANGE <= OFFSET, "Offset must be larger than the static range for efficient loading"

    def __init__(self, dynamic_source_manager, p):
        self.manager = dynamic_source_manager

        self.p = p
        self.x_range = p.x_range
        self.y_range = p.y_range

        self.fixed_static_bbox = None
        self._last_factor = 1
        self._last_lod = None
        self._last_valid_ranges = None
        self._update_pending = False
        # Pan/zoom loads run in a worker thread (only with a real Bokeh document,
        # whose next-tick callbacks run on the server loop; see add_range_callback)
        self.doc = None
        self.background = False

    def _load(self):
        if self.background:
            self.manager.fetch_data_async(self.doc)
        else:
            self.manager.fetch_data()

    @property
    def _ranges(self):
        return self.x_range.start, self.x_range.end, self.y_range.start, self.y_range.end

    def _ranges_valid(self):
        return viewport.ranges_valid(*self._ranges)

    def _check_ranges(self):
        """True si los rangos son válidos. Si BokehJS manda NaN (p.ej. mientras
        se recoloca el layout tras el clustering y el plot mide 0 px), se
        restauran los últimos rangos válidos y se devuelve False."""
        if self._ranges_valid():
            self._last_valid_ranges = (
                self.x_range.start, self.x_range.end,
                self.y_range.start, self.y_range.end,
            )
            return True
        if self._last_valid_ranges is not None:
            xs, xe, ys, ye = self._last_valid_ranges
            self.x_range.update(start=xs, end=xe)
            self.y_range.update(start=ys, end=ye)
        return False

    @property
    def window_size(self):
        return viewport.window_size(*self._ranges)

    @property
    def npixels(self):
        return self.window_size[0] * self.window_size[1]

    @property
    def x(self):
        return viewport.center(*self._ranges)[0]

    @property
    def y(self):
        return viewport.center(*self._ranges)[1]

    @property
    def static_bbox(self):
        return viewport.static_bbox(*self._ranges, self.STATIC_RANGE)

    @property
    def compute_bbox(self):
        return viewport.compute_bbox(*self._ranges, self.OFFSET)

    is_inside = staticmethod(viewport.is_inside)

    # Canvas width (px) assumed before the browser has reported the real one
    DEFAULT_CANVAS_PX = settings.DEFAULT_CANVAS_PX

    @staticmethod
    def _canvas_px(*getters):
        """First positive size among `getters` (inner_* are unset until rendered)."""
        for get in getters:
            try:
                value = get()
            except UnsetValueError:
                continue
            if value:
                return value
        return LazyLoader.DEFAULT_CANVAS_PX

    @property
    def pixel_size(self):
        """Figure units spanned by one screen pixel (the coarser axis)."""
        width = self._canvas_px(lambda: self.p.inner_width, lambda: self.p.width)
        height = self._canvas_px(lambda: self.p.inner_height, lambda: self.p.height, lambda: width)
        return viewport.pixel_size(*self._ranges, width, height)

    def get_resolution_level(self):
        sample = self.manager.sample
        # Cells: polygons simplified to about one screen pixel
        get_lod_level = getattr(sample, "get_lod_level", None)
        if get_lod_level is not None:
            self._last_lod = get_lod_level(self.pixel_size, current=self._last_lod)
            return self._last_lod
        # native_width/height ya vienen en el espacio (escalado) de la figura,
        # así que window_size y tamaño de bin son dimensionalmente consistentes.
        bin_w = getattr(sample, "native_width", None)
        bin_h = getattr(sample, "native_height", None)
        if not bin_w or not bin_h:
            return None

        max_span_in_bins = max(
            self.window_size[0] / bin_w,
            self.window_size[1] / bin_h,
        )
        factor, self._last_factor = viewport.aggregation_factor(
            max_span_in_bins, getattr(sample, "AGGREGATION_FACTORS", (1,)), self._last_factor,
            target_bins_per_axis=settings.TARGET_BINS_PER_AXIS, hysteresis_ratio=self.HYSTERESIS_RATIO,
            minimum=getattr(self.manager, 'minimum_aggregation_factor', None),
        )
        return sample.get_resolution_level_for_factor(factor)

    def add_range_callback(self, component, doc):
        self.doc = doc
        self.background = isinstance(doc, Document)

        def callback(attr, old, new):
            # A pan or zoom changes start and end of both ranges: update once for all
            if self._update_pending:
                return
            self._update_pending = True

            def update_bbox():
                self._update_pending = False
                # Hidden layer: nothing to load; showing it reloads the viewport
                if not self.manager.active or not self._check_ranges():
                    return
                resolution_level = self.get_resolution_level()
                current_fetch_bbox = self.compute_bbox
                n_before = len(self.manager.source.data.get('x', []))
                if self.manager.bbox is None:
                    self.manager.set_fetch_context(
                        bbox=current_fetch_bbox,
                        resolution_level=resolution_level,
                    )
                    self.fixed_static_bbox = self.static_bbox
                    self._load()
                else:
                    x, y = self.x, self.y
                    context_changed = self.manager.set_fetch_context(
                        bbox=self.manager.bbox,
                        resolution_level=resolution_level,
                    )
                    inside = LazyLoader.is_inside(x, y, self.fixed_static_bbox)
                    if context_changed or not inside:
                        self.manager.bbox = current_fetch_bbox
                        self.fixed_static_bbox = self.static_bbox
                        self._load()

                    if not inside and self.manager.needs_pruning():
                        # El radio debe cubrir TODA la región fetcheada (compute_bbox),
                        # cuyas esquinas están a ~1.2*window del centro. Usamos
                        # (1+OFFSET)*window para no recortar nunca bins visibles
                        # y solo podar los de un paneo lejano anterior.
                        prune_distance = viewport.prune_distance(*self._ranges, self.OFFSET)
                        self.manager.remove_further_data(x, y, distance=prune_distance)

            doc.add_next_tick_callback(update_bbox)

        self.x_range.on_change('start', callback)
        self.x_range.on_change('end', callback)
        self.y_range.on_change('start', callback)
        self.y_range.on_change('end', callback)

    def force_reload(self):
        """Reload current viewport data after source has been cleared."""
        if not self.manager.active or not self._check_ranges():
            return
        resolution_level = self.get_resolution_level()
        self.manager.set_fetch_context(
            bbox=self.compute_bbox,
            resolution_level=resolution_level,
        )
        self.manager.bbox = self.compute_bbox
        self.fixed_static_bbox = self.static_bbox
        self.manager.fetch_data()
