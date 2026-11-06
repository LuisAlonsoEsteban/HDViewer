"""Encoding of displayed values: alpha scaling and colour codes + palettes.

Bokeh-free: a colour is sent as a small integer code per row plus a palette
(list of hex colours); the Bokeh UI and the API each turn the palette into
their own colour mapper.
"""
import numpy as np
import pandas as pd

from .. import settings


class Parameter:

    def __init__(self, is_static):
        self.is_static = is_static


class ColorCodes:
    """Colours travel to the browser as small integer codes plus `palette`:
    changing a colour then sends the palette, not one hex string per row.
    Code 0 is the not-sequenced colour."""

    palette = None

    def _sync_palette(self, palette):
        palette = list(palette)
        if self.palette != palette:
            self.palette = palette
            self._palette_changed()

    def _palette_changed(self):
        """Hook for renderers that mirror the palette (e.g. a Bokeh colour mapper)."""

    def decode(self, codes):
        """Hex colours of `codes` (for the server side and tests)."""
        return np.asarray(self.palette, dtype=object)[np.asarray(codes, dtype=int)].tolist()


class StaticParameter(Parameter, ColorCodes):

    is_sequenced = None

    def __init__(self, value, attr='color', not_sequenced_val=settings.NOT_SEQUENCED_ATTRS["color"]):
        super().__init__(is_static=True)

        self.attr = attr

        self.value = value
        if self.is_color and not isinstance(not_sequenced_val, str):
            not_sequenced_val = value          # e.g. a numeric default for line_color
        self.not_sequenced_val = not_sequenced_val
        if self.is_color:
            self._sync_palette([not_sequenced_val, value])

    @property
    def is_color(self):
        return "color" in self.attr

    def set_is_sequenced(self, is_sequenced):
        self.is_sequenced = is_sequenced

    def _build_values(self, is_sequenced_values):
        sequenced_array = np.asarray(is_sequenced_values, dtype=bool)
        if sequenced_array.size == 0:
            return []
        if self.is_color:
            return sequenced_array.astype(np.uint8)
        return np.where(sequenced_array, self.value, self.not_sequenced_val)

    def get_values(self, is_sequenced_values:list=None):
        if is_sequenced_values is not None and len(is_sequenced_values) > 0:
            return self._build_values(is_sequenced_values)

        elif self.is_sequenced is not None and len(self.is_sequenced) > 0:
            return self._build_values(self.is_sequenced)

        else:
            return 1 if self.is_color else self.value

    def update(self, new_val=None, new_not_sequenced_val=None):
        if new_val is not None:
            self.value = new_val
        if new_not_sequenced_val is not None:
            self.not_sequenced_val = new_not_sequenced_val
        if self.is_color:
            self._sync_palette([self.not_sequenced_val, self.value])

        return self.get_values()


class DynamicParameter(Parameter):

    def __init__(
            self, 
            colname, 
            attr='color', 
            min_counts=0.0, 
            max_counts=None, 
            static_not_sequenced=True, 
            not_sequenced_val=settings.PARAMETER_NOT_SEQUENCED_COLOR,
            colname_hovered=False
            ):
        
        super().__init__(is_static=False)

        self.colname = colname
        self.attr = attr

        self.min_counts = min_counts
        self.max_counts = max_counts

        self.static_not_sequenced = static_not_sequenced
        self.not_sequenced_val = not_sequenced_val

        self.colname_hovered = colname_hovered

    def get_values(self, series):
        raise NotImplementedError


class Scaler(DynamicParameter):

    def __init__(self, minv, maxv, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.minv = minv
        self.maxv = maxv
        self.avoid_scaling = False

    def get_values(self, series):
        """Vectorised: values scaled to [minv, maxv] as float32; missing values
        (not sequenced) take `not_sequenced_val`."""
        series = pd.Series(series, copy=False)
        sequenced = series.notna().to_numpy()
        if self.avoid_scaling:
            # Uniform mode (categorical colouring): every value takes the slider max
            scaled = np.full(len(series), self.maxv, dtype=float)
        else:
            try:
                values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float)
            except (ValueError, TypeError):
                values = np.arange(len(series), dtype=float)
            range_val = self.max_counts - self.min_counts if self.max_counts != self.min_counts else 1
            scaled = self.minv + (values - self.min_counts) * (self.maxv - self.minv) / range_val
        if self.static_not_sequenced:
            scaled = np.where(sequenced, scaled, self.not_sequenced_val)
        return scaled.astype(np.float32)

    def update(self, series, new_minv=None, new_maxv=None, avoid_scaling=None):
        # avoid_scaling=None keeps the current mode, so moving the alpha slider
        # while clusters are shown does not bring back the old gene scaling
        if avoid_scaling is not None:
            self.avoid_scaling = avoid_scaling
        if new_minv is not None:
            self.minv = new_minv
            if new_minv > self.maxv:
                self.maxv = new_minv
        if new_maxv is not None:
            self.maxv = new_maxv
            if new_maxv < self.minv:
                self.minv = new_maxv
        return self.get_values(series)


class CategoricalColor(DynamicParameter, ColorCodes):

    PALETTE = settings.CATEGORICAL_PALETTE

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Codes index this order (+1, after the not-sequenced colour). It only
        # grows, so codes already sent to the browser stay valid.
        self._order = []
        self._colors = {}
        self._sync()

    @property
    def category_colors(self):
        return dict(self._colors)

    @category_colors.setter
    def category_colors(self, colors):
        for category, color in dict(colors).items():
            if category not in self._colors:
                self._order.append(category)
            self._colors[category] = color
        self._sync()

    def _sync(self):
        self._sync_palette([self.not_sequenced_val, *(self._colors[c] for c in self._order)])

    def ensure_colors(self, categories):
        missing = pd.Index(categories).difference(self._colors.keys())
        for category in missing:
            self._order.append(category)
            self._colors[category] = self.PALETTE[(len(self._order) - 1) % len(self.PALETTE)]
        if len(missing):
            self._sync()

    def get_values(self, series):
        """Vectorised: one uint16 code per row (0 = not sequenced)."""
        series = pd.Series(series, copy=False)
        sequenced = series.notna().to_numpy()
        labels = series[sequenced].astype(str)
        self.ensure_colors(pd.unique(labels))
        codes = np.zeros(len(series), dtype=np.uint16)
        codes[sequenced] = pd.Index(self._order).get_indexer(labels) + 1
        return codes
