"""Extra style parameters of a layer (the Parameters tab), without any UI.

Besides its colour and opacity (CORE_ATTRS), a layer's glyphs can bind
`EXTRA_ATTRS` to a fixed value or to data: colours to a categorical metadata
column (one colour per category), widths and opacities to a gene or numeric
metadata column scaled to an output range.
"""
import pandas as pd

from .. import settings
from .encoding import CategoricalColor, Scaler, StaticParameter

CORE_ATTRS = ("color", "alpha")
EXTRA_ATTRS = ("fill_alpha", "line_width", "line_color", "line_alpha")


def is_color(attr):
    return "color" in attr


def attr_spec(attr):
    """Input limits and defaults per kind of attribute."""
    if is_color(attr):
        return dict(settings.PARAMETER_COLOR_SPEC)
    if "alpha" in attr:
        return dict(settings.PARAMETER_ALPHA_SPEC)
    return dict(settings.PARAMETER_WIDTH_SPEC)  # widths


def column_choices(sample, attr):
    """Columns an attribute can follow: categorical metadata for colours; genes
    and numeric metadata (what Sample.is_continous accepts) otherwise."""
    adata = sample.adata
    numeric = {c: pd.api.types.is_numeric_dtype(adata.obs[c]) for c in adata.obs.columns}
    if is_color(attr):
        return [c for c, is_num in numeric.items() if not is_num]
    return list(adata.var_names) + [c for c, is_num in numeric.items() if is_num]


def make_parameter(sample, attr, *, mode, value=None, column=None, minv=None, maxv=None,
                   not_sequenced=None, hovered=False):
    """The parameter for one attribute; raises ValueError on invalid input."""
    if attr not in EXTRA_ATTRS:
        raise ValueError(f"Choose one of {', '.join(EXTRA_ATTRS)}.")
    color = is_color(attr)
    spec = attr_spec(attr)
    if not_sequenced is None:
        not_sequenced = spec["not_sequenced"]
    if mode == "static":
        if value is None or value == "":
            raise ValueError("Enter a value.")
        if not color:
            value = float(value)
            if not spec.get("low", value) <= value <= spec.get("high", value):
                raise ValueError(f"The value must be between {spec['low']} and {spec['high']}.")
        return StaticParameter(value, attr=attr, not_sequenced_val=not_sequenced)
    if mode != "data":
        raise ValueError("mode must be 'static' or 'data'.")
    if column not in column_choices(sample, attr):
        raise ValueError("Choose a categorical metadata column." if color
                         else "Choose a gene or a numeric metadata column.")
    if color:
        return CategoricalColor(column, attr=attr, static_not_sequenced=True,
                                not_sequenced_val=not_sequenced, colname_hovered=hovered)
    minv = spec["minv"] if minv is None else float(minv)
    maxv = spec["maxv"] if maxv is None else float(maxv)
    return Scaler(minv, maxv, column, attr=attr, static_not_sequenced=True,
                  not_sequenced_val=float(not_sequenced), colname_hovered=hovered)


def describe(param):
    entry = dict(attr=param.attr, not_sequenced=param.not_sequenced_val)
    if isinstance(param, StaticParameter):
        entry.update(mode="static", value=param.value)
    else:
        entry.update(mode="data", column=param.colname, hovered=bool(param.colname_hovered))
        if isinstance(param, Scaler):
            entry.update(min=param.minv, max=param.maxv)
    return entry
