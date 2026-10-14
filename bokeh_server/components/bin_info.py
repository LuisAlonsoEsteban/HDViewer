"""A screen-anchored readout of the bins currently displayed on the main plot."""

import math
from html import escape

from bokeh.layouts import column
from bokeh.models import Div, Node, Panel, Select

from ..loader import LazyLoader


class BinInfoOverlay:
    def __init__(self, component):
        self.component = component
        self._renderer = None
        self._native_geometry = self._measure_native_geometry()
        self._syncing = False
        self.mode_select = Select(
            title="Bin resolution", value="dynamic", sizing_mode="stretch_width",
            options=[("dynamic", "Dynamic · follow zoom"), ("1", "Fixed · raw/native")]
                    + [(str(factor), f"Fixed · {factor} × {factor}")
                       for factor in component.sample.AGGREGATION_FACTORS if factor > 1],
            stylesheets=["""
                :host { color: #e5edf0; font-size: 12px; }
                label { color: #a9c7c9; }
                .bk-input { background: #23383e; color: #e5edf0;
                    border: 1px solid #5b8585; border-radius: 6px; }
                .bk-input:focus-visible { outline: 2px solid #5ee0c2; outline-offset: 2px; }
            """],
        )
        self.label = Div(name="bin_resolution_details", sizing_mode="stretch_width")
        self.lock_notice = Div(
            text="Clustering fixes the resolution. Display another column to change it.",
            visible=False, sizing_mode="stretch_width",
            styles=dict(color="#a9c7c9", font_size="11px"),
        )
        self.panel = Panel(
            name="bin_resolution_info",
            position=Node(target="frame", symbol="top_right"),
            anchor="top_right", width=270,
            elements=[column(self.mode_select, self.label, self.lock_notice,
                             width=244, spacing=6, margin=0)],
            styles=dict(background="rgba(24,39,45,0.94)", color="#e5edf0",
                        border="1px solid #5b8585", border_radius="8px",
                        padding="12px", width="270px", margin_left="-12px", margin_top="12px",
                        box_sizing="border-box", pointer_events="auto"),
        )
        component.figure.elements.append(self.panel)
        self.mode_select.on_change("value", self._mode_changed)
        component.manager.source.on_change("data", self._source_changed)
        # Empty sources may remain identical when resolution or cluster locking
        # changes, so source events alone are not sufficient.
        component.manager.on_context_change(self.update)
        self.update()

    def _measure_native_geometry(self):
        sample = self.component.sample
        if sample.gdf.empty:
            return None
        # All native footprints use the same fitted grid basis. Measure an
        # actual edge after scaling, rather than the obsolete rectangle angle.
        corners = list(sample.gdf.geometry.iloc[0].exterior.coords)
        dx, dy = corners[1][0] - corners[0][0], corners[1][1] - corners[0][1]
        hx, hy = corners[2][0] - corners[1][0], corners[2][1] - corners[1][1]
        angle = (math.degrees(math.atan2(dy, dx)) + 90) % 180 - 90
        if abs(angle) < 0.005:
            angle = 0.0
        return math.hypot(dx, dy), math.hypot(hx, hy), angle

    def metadata(self):
        manager = self.component.manager
        sample = self.component.sample
        factor = sample.get_aggregation_factor(manager.effective_resolution_level)
        if factor > 1:
            geometry = sample.native_width * factor, sample.native_height * factor, 0.0
        else:
            geometry = self._native_geometry
        return {
            "factor": factor,
            "geometry": geometry,
            "bin_size_um": sample.bin_size_um,
            "locked": sample.get_forced_aggregation_factor(manager.get_colnames()) is not None,
        }

    def update(self):
        info = self.metadata()
        factor = info["factor"]
        title = f"BINS · {factor} × {factor}"
        if factor == 1:
            title += " · native"
        lines = [title]
        if info["geometry"] is None:
            lines.append("No bin geometry available")
        else:
            width, height, angle = info["geometry"]
            lines.extend((f"Width      {width:,.2f} plot units",
                          f"Height     {height:,.2f} plot units",
                          f"Rotation   {angle:.2f}°"))
        if info["bin_size_um"] is not None:
            lines.append(f"Native bin {info['bin_size_um']:g} µm")
        fixed = self.component.manager.fixed_aggregation_factor
        lines.append("Resolution Cluster-locked" if info["locked"] else
                     "Resolution Fixed" if fixed is not None else "Resolution Zoom-driven")
        self._syncing = True
        try:
            self.mode_select.value = (str(factor) if info["locked"] else
                                      "dynamic" if fixed is None else str(fixed))
            self.mode_select.disabled = info["locked"]
            self.lock_notice.visible = info["locked"]
        finally:
            self._syncing = False
        text = ('<div style="white-space:pre-wrap;font:11px/1.5 monospace;color:#e5edf0">'
                + escape("\n".join(lines)) + '</div>')
        if self.label.text != text:
            self.label.text = text

    def _mode_changed(self, attr, old, new):
        if self._syncing:
            return
        if self.metadata()["locked"]:
            self.update()
            return
        factor = None if new == "dynamic" else int(new)
        self.component.manager.set_aggregation_factor(factor)
        # Reload immediately, using current ranges even if the user has not
        # panned. The persistent loader keeps tracking zoom while fixed.
        loader = getattr(self.component, "loader", None)
        if loader is None:
            loader = LazyLoader(self.component.manager, self.component.figure)
        loader.force_reload()

    def _source_changed(self, attr, old, new):
        self.update()

    def _visibility_changed(self, attr, old, new):
        self.panel.visible = new

    def bind_renderer(self, renderer):
        if self._renderer is not None:
            self._renderer.remove_on_change("visible", self._visibility_changed)
        self._renderer = renderer
        renderer.on_change("visible", self._visibility_changed)
        self.panel.visible = renderer.visible
        self.update()
