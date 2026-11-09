"""One map layer as seen by one browser session.

`LayerView` wraps the shared core.layer_state.LayerRows (with a recording sink)
and adds the per-layer interface state: the view geometry
decisions of the lazy loader, the colouring choices of the Color tab and the
native selection. Every change returns the recorded operations, which the API
sends as a binary frame (api.frames).
"""
import math
import threading

import numpy as np
import pandas as pd

from .. import settings
from ..core import viewport
from ..core.clustering import ClusteringState
from ..core.encoding import CategoricalColor, Scaler, StaticParameter
from ..core.groups import GROUPS_COL, GroupsState
from ..core import parameters as extra_params
from ..core.layer_state import LayerRows
from ..core.roi import ROIGeometry

KIND_DEFAULTS = {
    "bins": dict(not_sequenced_alpha=settings.NOT_SEQUENCED_ATTRS["alpha"], hover="green"),
    "spots": dict(not_sequenced_alpha=settings.NOT_SEQUENCED_ATTRS["alpha"], hover="green"),
    "cells": dict(not_sequenced_alpha=0.5, hover="blue"),
}


class LayerView:
    def __init__(self, layer_id, kind, sample, title, colour):
        self.id, self.kind, self.title = layer_id, kind, title
        self.sample = sample
        self.lock = threading.RLock()
        self.rows = LayerRows(sample)
        self.visible = True
        self.colour = colour
        self.alpha_max = 1.0
        self.colour_by = dict(kind="gene", name=sample.get_colnames()[0])
        self.seq = -1
        self.version = 0                 # data version: bumped whenever ops are sent
        self._last_factor = 1
        self._last_lod = None
        self._static_bbox = None
        self.selected = set()
        self.clustering = ClusteringState(sample)
        self.groups = GroupsState(sample)
        self.rois = ROIGeometry(sample)
        self.roi_preview = dict(generation=0, key=None, result=None)   # last group-seed preview
        self.de = dict(header=None, merged=None, plot_data=None, names=None, colours=None, generation=0,
                       running=None, compared=None)      # (source, names) being / last compared
        self.job = dict(name=None, state="idle", message="")      # clustering or DE (they exclude each other)
        self.deconvolution_job = dict(name="deconvolution", state="idle", message="")   # may take hours
        self.deconvolution_run = None   # the running deconvolution's handle (api.routers.deconvolution)
        # Last metadata comparison (Bokeh json_item); `generation` changes with groups and clusters,
        # so a comparison computed meanwhile is discarded
        self.metadata = dict(item=None, compared=None, generation=0)
        # Last embedding: method, embedded native obs names and coordinates; `version`
        # changes with the hover values (recolouring)
        self.embedding = dict(method=None, names=None, x=None, y=None, values=None, version=0)

        defaults = KIND_DEFAULTS[kind]
        self.rows.add_parameter(Scaler(0.0, 1.0, self.colour_by["name"], attr="alpha",
                                       not_sequenced_val=defaults["not_sequenced_alpha"], colname_hovered=True))
        self.rows.add_parameter(StaticParameter(colour, attr="color",
                                                not_sequenced_val=settings.NOT_SEQUENCED_ATTRS["color"]))
        if kind == "spots":
            # The circle glyph reads a per-row radius: declare it before the first
            # fetch, so even an empty state matches the figure's source
            self.rows.coordinate_columns = ("x", "y", "radius")
            if "radius" not in self.rows.sink.data:
                self.rows.sink.data["radius"] = []
        self.rows.sink.take()

    # -- description -----------------------------------------------------------

    def describe(self):
        return dict(id=self.id, kind=self.kind, title=self.title, visible=self.visible,
                    colour=self.colour, alpha=self.alpha_max, colour_by=self.colour_by,
                    scale=self.sample.expression_scale,
                    hover_colour=KIND_DEFAULTS[self.kind]["hover"],
                    aggregation_factors=list(getattr(self.sample, "AGGREGATION_FACTORS", ())),
                    fixed_aggregation=self.rows.fixed_aggregation_factor)

    def _value_label(self, name):
        """A column's tooltip label; genes name the layer's expression scale."""
        scale = self.sample.expression_scale
        if name in self.sample.adata.var_names and scale != "counts":
            return name + (" (log-normalized)" if scale == "log1p" else " (normalized)")
        return name

    def tooltips(self):
        rows = [["Barcode", "$index{custom}"]]
        shown = set()
        name = self.colour_by.get("name")
        if self.colour_by["kind"] != "static" and name in self.rows.sink.data:
            rows.append([self._value_label(name), f"@{{{name}}}"])
            shown.add(name)
        for param in self.extra_parameters():         # Parameters tab columns shown on hover
            column = getattr(param, "colname", None)
            if (not param.is_static and param.colname_hovered and column in self.rows.sink.data
                    and column not in shown):
                rows.append([self._value_label(column), f"@{{{column}}}"])
                shown.add(column)
        return rows

    # -- extra glyph parameters (Parameters tab) -------------------------------------

    def extra_parameters(self):
        return [param for param in self.rows.parameters if param.attr in extra_params.EXTRA_ATTRS]

    def parameter_palettes(self):
        return {param.attr: list(param.palette) for param in self.extra_parameters()
                if extra_params.is_color(param.attr) and getattr(param, "palette", None) is not None}

    def describe_parameters(self):
        current = {param.attr: extra_params.describe(param) for param in self.extra_parameters()}
        return [dict(attr=attr, color=extra_params.is_color(attr), spec=extra_params.attr_spec(attr),
                     current=current.get(attr)) for attr in extra_params.EXTRA_ATTRS]

    def set_parameter(self, attr, **options):
        """Bind `attr` to a fixed value or a column (raises ValueError). The glyph
        changes, so the browser rebuilds the map figure afterwards."""
        param = extra_params.make_parameter(self.sample, attr, **options)
        self.rows.add_parameter(param)
        self.rows.sink.take()
        return extra_params.describe(param)

    def remove_parameter(self, attr):
        if attr not in extra_params.EXTRA_ATTRS:
            raise ValueError(f"Choose one of {', '.join(extra_params.EXTRA_ATTRS)}.")
        self.rows.remove_parameter(attr)
        self.rows.sink.take()

    def header(self, ops):
        base = self.version
        if ops:
            self.version += 1
        rows = self.rows
        colour = rows.get_param("color")
        return dict(layer=self.id, base=base, version=self.version, row_version=rows.row_version,
                    rows=len(rows.indexes), level=rows.effective_resolution_level,
                    constants=dict(rows.glyph_constants), palette=list(colour.palette),
                    palettes=self.parameter_palettes(),
                    tooltips=self.tooltips(), visible=self.visible, bins=self.bin_info())

    # -- bin resolution readout (the map's overlay) ------------------------------------

    def _native_geometry(self):
        """(width, height, rotation°) of a native bin, measured on its footprint."""
        if not hasattr(self, "_native"):
            self._native = None
            gdf = self.sample.gdf
            if not gdf.empty:
                corners = list(gdf.geometry.iloc[0].exterior.coords)
                dx, dy = corners[1][0] - corners[0][0], corners[1][1] - corners[0][1]
                hx, hy = corners[2][0] - corners[1][0], corners[2][1] - corners[1][1]
                angle = (math.degrees(math.atan2(dy, dx)) + 90) % 180 - 90
                self._native = (math.hypot(dx, dy), math.hypot(hx, hy), 0.0 if abs(angle) < 0.005 else angle)
        return self._native

    def bin_info(self):
        """What the bins overlay shows: the displayed aggregation and its geometry."""
        if self.kind != "bins":
            return None
        rows, sample = self.rows, self.sample
        factor = sample.get_aggregation_factor(rows.effective_resolution_level)
        native = self._native_geometry()
        if factor > 1:
            rotation = native[2] if native and getattr(sample, "spatial_transform", None) else 0.0
            geometry = (sample.native_width * factor, sample.native_height * factor, rotation)
        else:
            geometry = native
        return dict(factor=int(factor), width=geometry and float(geometry[0]), height=geometry and float(geometry[1]),
                    rotation=geometry and float(geometry[2]), bin_size_um=getattr(sample, "bin_size", None),
                    fixed=rows.fixed_aggregation_factor, minimum=rows.minimum_aggregation_factor,
                    factors=[int(f) for f in getattr(sample, "AGGREGATION_FACTORS", ())])

    def visibility(self):
        """Rows to draw: those of the shown clusters or groups (None: every row)."""
        rows = self.rows
        mask = np.asarray(rows.view_filter, dtype=bool)
        kind, column = self.colour_by.get("kind"), self.colour_by.get("name")
        labels = rows.sink.data.get(column)
        if kind in ("clusters", "groups") and labels is not None and len(labels) == len(mask):
            labels = pd.Series(labels, dtype=object)
            if kind == "clusters":
                shown = self.clustering.visible_mask(labels)
            else:
                shown = labels.isin(self.groups._active_groups).to_numpy()
            mask = mask & shown
        return None if mask.all() else mask.astype(np.uint8)

    def frame_parts(self, ops):
        """(header, ops, extra arrays) of a response."""
        mask = self.visibility()
        header = self.header(ops)
        header["filtered"] = mask is not None
        return header, ops, ({"visible": mask} if mask is not None else None)

    def full_state(self):
        """Everything the browser needs to rebuild this layer's data source."""
        return self.frame_parts([("replace", dict(self.rows.sink.data))])

    # -- viewport ------------------------------------------------------------------

    def _resolution_level(self, x0, x1, y0, y1, width_px, height_px):
        sample = self.sample
        get_lod_level = getattr(sample, "get_lod_level", None)
        if get_lod_level is not None:
            self._last_lod = get_lod_level(viewport.pixel_size(x0, x1, y0, y1, width_px, height_px),
                                           current=self._last_lod)
            return self._last_lod
        bin_w, bin_h = getattr(sample, "native_width", None), getattr(sample, "native_height", None)
        if not bin_w or not bin_h:
            return None
        width, height = viewport.window_size(x0, x1, y0, y1)
        factor, self._last_factor = viewport.aggregation_factor(
            max(width / bin_w, height / bin_h), getattr(sample, "AGGREGATION_FACTORS", (1,)), self._last_factor,
            target_bins_per_axis=settings.TARGET_BINS_PER_AXIS, hysteresis_ratio=settings.HYSTERESIS_RATIO,
            minimum=self.rows.minimum_aggregation_factor)
        return sample.get_resolution_level_for_factor(factor)

    def update_viewport(self, x0, x1, y0, y1, width_px, height_px):
        """Load what the view needs (same decisions as the Bokeh lazy loader)."""
        rows = self.rows
        if not self.visible or not viewport.ranges_valid(x0, x1, y0, y1):
            return rows.sink.take()
        level = self._resolution_level(x0, x1, y0, y1, width_px, height_px)
        x, y = viewport.center(x0, x1, y0, y1)
        compute_bbox = viewport.compute_bbox(x0, x1, y0, y1, settings.OFFSET)
        if rows.bbox is None:
            rows.set_fetch_context(bbox=compute_bbox, resolution_level=level)
            self._static_bbox = viewport.static_bbox(x0, x1, y0, y1, settings.STATIC_RANGE)
            rows.fetch_data()
        else:
            changed = rows.set_fetch_context(bbox=rows.bbox, resolution_level=level)
            inside = self._static_bbox is not None and viewport.is_inside(x, y, self._static_bbox)
            if changed or not inside:
                rows.bbox = compute_bbox
                self._static_bbox = viewport.static_bbox(x0, x1, y0, y1, settings.STATIC_RANGE)
                rows.fetch_data()
            if not inside and rows.needs_pruning():
                rows.remove_further_data(x, y, distance=viewport.prune_distance(x0, x1, y0, y1, settings.OFFSET))
        return rows.sink.take()

    # -- style -------------------------------------------------------------------------

    def set_visible(self, visible):
        self.visible = visible
        self.rows.set_active(visible)
        if not visible:
            self._static_bbox = None
        return self.rows.sink.take()

    def set_scale(self, scale):
        with self.sample.fetch_lock:
            self.sample.expression_scale = scale
        self.rows.refresh_values()
        return self.rows.sink.take()

    def set_aggregation(self, factor):
        self.rows.set_aggregation_factor(factor)
        self._static_bbox = None
        return self.rows.sink.take()

    def set_colour(self, colour):
        self.colour = colour
        param = self.rows.get_param("color")
        if isinstance(param, StaticParameter):
            self.rows.update_param("color", new_val=colour)
        return self.rows.sink.take()

    def set_alpha(self, alpha):
        self.alpha_max = alpha
        self.rows.update_param("alpha", new_maxv=alpha)
        return self.rows.sink.take()

    def colour_by_column(self, name):
        """Gene or metadata column: continuous values scale the opacity of one
        colour; categories get one colour each (as in the Bokeh Color tab)."""
        rows, sample = self.rows, self.sample
        if name not in sample.get_colnames():
            raise KeyError(name)
        alpha = rows.get_param("alpha")
        if sample.is_continous(name):
            colour = rows.get_param("color")
            if not isinstance(colour, StaticParameter):
                rows.add_parameter(StaticParameter(self.colour, attr="color",
                                                   not_sequenced_val=settings.NOT_SEQUENCED_ATTRS["color"]))
            if alpha.colname != name:
                rows.update_colname(alpha.colname, name, attrs=("alpha",))
            rows.update_param("alpha", new_minv=0, new_maxv=self.alpha_max, avoid_scaling=False)
            kind = "gene" if name in sample.adata.var_names else "metadata"
        else:
            colour = rows.get_param("color")
            if not isinstance(colour, CategoricalColor) or colour.colname != name:
                rows.add_parameter(CategoricalColor(name, attr="color", colname_hovered=True,
                                                    not_sequenced_val=settings.NOT_SEQUENCED_ATTRS["color"]))
            rows.update_param("alpha", avoid_scaling=True)
            kind = "category"
        self.colour_by = dict(kind=kind, name=name)
        return rows.sink.take()

    def colour_by_categories(self, column, colours, kind):
        """Colour by a cluster or group column with its assigned colours."""
        rows, sample = self.rows, self.sample
        if hasattr(sample, "clear_multiscale_cache"):
            sample.clear_multiscale_cache()
        param = CategoricalColor(column, attr="color", colname_hovered=True,
                                 not_sequenced_val=settings.NOT_SEQUENCED_ATTRS["color"])
        param.category_colors = dict(colours)
        rows.add_parameter(param)
        rows.update_param("alpha", avoid_scaling=True)
        self.colour_by = dict(kind=kind, name=column)
        self._static_bbox = None
        return rows.sink.take()

    def recolour_categories(self, colours):
        param = self.rows.get_param("color")
        if isinstance(param, CategoricalColor):
            param.category_colors = dict(colours)
        return self.rows.sink.take()

    # -- differential expression -----------------------------------------------------

    def drop_comparisons(self, clusters=False, groups=(), groups_changed=False):
        """Forget the DE results (and any running DE) and the metadata comparison
        that compared re-clustered clusters or one of the removed `groups`; a
        metadata comparison against The rest also depends on every group
        (`groups_changed`). True if the DE results were dropped."""
        if clusters or groups or groups_changed:
            self.metadata["generation"] += 1
        compared = self.metadata["compared"]
        if compared and ((clusters and any(source == "clusters" for source, _ in compared["chosen"]))
                         or any(name in groups for _, name in compared["chosen"])
                         or (compared["include_rest"] and (groups or groups_changed))):
            self.metadata.update(item=None, compared=None)
        de = self.de

        def affected(compared):
            if not compared:
                return False
            source, names = compared
            return (clusters and source in ("clusters", "mixed")) or any(name in names for name in groups)

        if affected(de.get("running")):        # its result would be stale on arrival
            de["generation"] += 1
            de["running"] = None
        if not affected(de.get("compared")):
            return False
        de.update(header=None, merged=None, plot_data=None, names=None, colours=None, compared=None)
        return True

    # -- names and selection -------------------------------------------------------

    def names(self, row_version, positions):
        if row_version != self.rows.row_version:
            return []
        names = self.rows.loaded_obs_names
        return [[int(i), str(names[i])] for i in positions if isinstance(i, int) and 0 <= i < len(names)]

    def select(self, row_version, positions):
        from ..core.selection import native_selection
        if row_version != self.rows.row_version:
            raise ValueError("stale rows")
        names = self.rows.loaded_obs_names
        picked = [names[i] for i in positions if isinstance(i, int) and 0 <= i < len(names)]
        self.selected = native_selection(self.sample, picked)
        return len(self.selected)

    def selected_positions(self):
        from ..core.selection import selection_indices
        return selection_indices(self.sample, self.selected, self.rows.loaded_obs_names)

    def column_names(self, query="", limit=50):
        query = query.lower()
        names = self.sample.get_colnames()
        matches = [name for name in names if name.lower().startswith(query)] if query else names
        return matches[:limit]


def build_views(loaded, um_per_px=None):
    """One view per loaded layer; `um_per_px` calibrates ROI radii in µm."""
    from types import SimpleNamespace
    views, cells = {}, 0
    for layer_id, (kind, sample, title) in loaded.samples.items():
        colour = settings.INITIAL_ATTRS["color"]
        if kind == "cells":
            colour = settings.CELL_LAYER_COLORS[cells % len(settings.CELL_LAYER_COLORS)]
            cells += 1
        view = views[layer_id] = LayerView(layer_id, kind, sample, title, colour)
        if um_per_px:
            view.rois.coordinate_scale = SimpleNamespace(data=dict(px_x=[1.0], px_y=[1.0], um_per_px=[um_per_px]))
    return views

