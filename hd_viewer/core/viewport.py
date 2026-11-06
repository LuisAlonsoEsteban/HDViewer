"""Viewport geometry for lazy loading, independent of any plotting library.

Ranges are (start, end) in figure units. The Bokeh lazy loader and the HTTP API
both decide what to fetch with these functions.
"""
import math


def ranges_valid(x0, x1, y0, y1):
    values = (x0, x1, y0, y1)
    if any(v is None or not math.isfinite(v) for v in values):
        return False
    return x1 > x0 and y1 > y0


def window_size(x0, x1, y0, y1):
    return x1 - x0, y1 - y0


def center(x0, x1, y0, y1):
    return (x0 + x1) / 2, (y0 + y1) / 2


def scaled_bbox(x0, x1, y0, y1, scale):
    """[left, right, bottom, top] of the window scaled by `scale` around its centre."""
    (width, height), (x, y) = window_size(x0, x1, y0, y1), center(x0, x1, y0, y1)
    return [x - width * scale / 2, x + width * scale / 2, y - height * scale / 2, y + height * scale / 2]


def static_bbox(x0, x1, y0, y1, static_range):
    """While the view centre stays inside this box, no new fetch is needed."""
    return scaled_bbox(x0, x1, y0, y1, static_range)


def compute_bbox(x0, x1, y0, y1, offset):
    """Region fetched for a view: the window plus a margin of `offset` windows."""
    return scaled_bbox(x0, x1, y0, y1, 1 + offset)


def is_inside(x, y, bbox):
    return bbox[0] < x < bbox[1] and bbox[2] < y < bbox[3]


def prune_distance(x0, x1, y0, y1, offset):
    """Rows farther than this from the view centre may be dropped: it covers the
    whole fetched region (whose corners are ~1.2 windows away)."""
    return max(window_size(x0, x1, y0, y1)) * (1 + offset)


def pixel_size(x0, x1, y0, y1, width_px, height_px):
    """Figure units spanned by one screen pixel (the coarser axis)."""
    width, height = window_size(x0, x1, y0, y1)
    return max(width / width_px, height / height_px)


def neighbour_factors(factor, factors):
    """(previous, next) aggregation factors around `factor`."""
    factors = tuple(factors)
    if factor not in factors:
        return None, None
    i = factors.index(factor)
    return (factors[i - 1] if i > 0 else None), (factors[i + 1] if i < len(factors) - 1 else None)


def aggregation_factor(max_span_in_bins, factors, last_factor, *, target_bins_per_axis,
                       hysteresis_ratio, minimum=None):
    """Bin aggregation factor for a view spanning `max_span_in_bins` native bins.

    The smallest factor keeping at most `target_bins_per_axis` bins per axis, with
    hysteresis around `last_factor` so zooming near a threshold does not flip levels.
    Factors below `minimum` (e.g. clusters computed on aggregated bins) are skipped.
    Returns (factor, last_factor to remember).
    """
    factors = all_factors = tuple(factors)
    if minimum is not None:
        factors = tuple(value for value in factors if value >= minimum) or (minimum,)
        if last_factor < factors[0]:
            last_factor = factors[0]

    factor = factors[-1]
    for candidate in factors:
        if max_span_in_bins / candidate <= target_bins_per_axis:
            factor = candidate
            break

    prev_factor, next_factor = neighbour_factors(last_factor, all_factors)
    if factor != last_factor:
        if factor > last_factor and next_factor is not None:
            promote_at = target_bins_per_axis * next_factor * (1 + hysteresis_ratio)
            if max_span_in_bins < promote_at:
                factor = last_factor
        elif factor < last_factor and prev_factor is not None:
            demote_at = target_bins_per_axis * last_factor * (1 - hysteresis_ratio)
            if max_span_in_bins > demote_at:
                factor = last_factor
    return factor, factor
