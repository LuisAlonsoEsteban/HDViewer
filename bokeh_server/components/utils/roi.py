"""Per-layer, session-only spatial ROIs and shared map drawing tools."""

from dataclasses import dataclass, field
from html import escape
import uuid
import threading

import numpy as np
import shapely
from shapely.geometry import Point, Polygon, box
from shapely.affinity import scale
from bokeh.events import SelectionGeometry, Tap
from bokeh.layouts import column, row
from bokeh.models import BoxSelectTool, Button, Checkbox, ColumnDataSource, Div, LassoSelectTool, NumericInput, Select, TapTool, TextInput

from .loading import LoadingLogger
from .palettes import GROUP_PALETTE, palette_color
from src.utils.progress import report


@dataclass(eq=False)
class ROI:
    name: str
    geometry: object
    kind: str
    color: str
    id: str = field(default_factory=lambda: uuid.uuid4().hex)
    mask: np.ndarray | None = None
    visible: bool = True
    error: str | None = None
    seed_group: str | None = None
    seed_count: int | None = None

    @property
    def count(self):
        return None if self.mask is None else int(self.mask.sum())


class _DrawingController:
    """One set of drawing tools/handlers per shared figure, with one active owner."""

    def __init__(self, figure):
        self.figure = figure
        self.owner = None
        self.box_tool = BoxSelectTool(renderers=[], continuous=False, description="Draw ROI bounding box")
        self.lasso_tool = LassoSelectTool(renderers=[], continuous=False, description="Draw ROI lasso")
        self.seed_tool = TapTool(renderers=[], description="Set ROI seed")
        figure.add_tools(self.box_tool, self.lasso_tool, self.seed_tool)
        figure.on_event(SelectionGeometry, self._geometry)
        figure.on_event(Tap, self._tap)

    def arm(self, owner, mode):
        if self.owner is not None:
            self.owner.cancel()
        self.owner, self.mode = owner, mode
        toolbar = self.figure.toolbar
        self.previous = toolbar.active_drag, toolbar.active_tap
        toolbar.active_drag = {"box": self.box_tool, "lasso": self.lasso_tool}.get(mode)
        toolbar.active_tap = self.seed_tool if mode == "radius" else None

    def release(self, owner):
        if self.owner is owner:
            self.figure.toolbar.active_drag, self.figure.toolbar.active_tap = self.previous
            self.owner = None

    def _geometry(self, event):
        owner = self.owner
        if owner is None or not event.final or not event.geometry:
            return
        geometry = event.geometry
        expected = "rect" if self.mode == "box" else "poly"
        if self.mode == "radius" or geometry.get("type") != expected:
            return
        if self.figure.toolbar.active_drag is not (self.box_tool if self.mode == "box" else self.lasso_tool):
            return
        try:
            if self.mode == "box":
                shape = owner.box_geometry(geometry['x0'], geometry['y0'], geometry['x1'], geometry['y1'])
            else:
                shape = owner.lasso_geometry(geometry['x'], geometry['y'])
            owner.set_draft(shape, self.mode)
        except (ValueError, KeyError, TypeError) as exc:
            owner._error(exc)
        finally:
            self.release(owner)

    def _tap(self, event):
        owner = self.owner
        if owner is None or self.mode != "radius" or self.figure.toolbar.active_tap is not self.seed_tool:
            return
        try:
            factors, _ = owner.radius_scale()
            shape = owner.seed_radius_geometry(event.x, event.y, owner.radius.value, *factors)
            owner.set_draft(shape, "radius")
        except (ValueError, TypeError) as exc:
            owner._error(exc)
        finally:
            self.release(owner)


class ROIManager:
    def __init__(self, component):
        self.component = component
        self.rois = []
        self.coordinate_scale = None
        self.doc = None
        self.logger = None
        self.draft = None
        self.draft_kind = None
        self.visible = True
        self.disposed = False
        self._layout = None
        self._counter = 1
        self._preview_generation = 0
        self._preview_cancel = None
        self._preview_timeout = None
        self._pending_preview = None
        self._draft_seed_info = None
        self._requested_seed_mask = None
        figure = component.figure
        if not hasattr(figure, '_roi_drawing_controller'):
            figure._roi_drawing_controller = _DrawingController(figure)
        self.drawing = figure._roi_drawing_controller
        self.source = ColumnDataSource(dict(xs=[], ys=[], color=[]))
        self.renderer = figure.multi_polygons('xs', 'ys', source=self.source, fill_alpha=.06,
                                       line_color='color', fill_color='color', line_width=2,
                                       selection_fill_alpha=.06, nonselection_fill_alpha=.06,
                                       level='overlay', name=f'{component.key}_rois')
        component.groups_manager.add_change_listener(self._groups_changed)

    @property
    def sample(self):
        return self.component.sample

    @staticmethod
    def validate_geometry(geometry):
        if geometry.is_empty or geometry.geom_type not in ('Polygon', 'MultiPolygon') or not geometry.is_valid or geometry.area <= 0:
            raise ValueError("Draw a nonempty, non-self-intersecting region with positive area.")
        if not np.isfinite(shapely.get_coordinates(geometry)).all():
            raise ValueError("ROI coordinates must be finite.")
        return geometry

    @classmethod
    def box_geometry(cls, x0, y0, x1, y1):
        if not np.isfinite([x0, y0, x1, y1]).all() or x0 == x1 or y0 == y1:
            raise ValueError("A bounding box needs finite coordinates and positive width and height.")
        return cls.validate_geometry(box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))

    @classmethod
    def lasso_geometry(cls, xs, ys):
        if len(xs) != len(ys) or len(xs) < 3 or not np.isfinite([xs, ys]).all():
            raise ValueError("A lasso needs at least three finite vertices.")
        return cls.validate_geometry(Polygon(zip(xs, ys)))

    @classmethod
    def radius_geometry(cls, x, y, radius, units_x=1., units_y=1.):
        if not np.isfinite([x, y, radius, units_x, units_y]).all() or min(radius, units_x, units_y) <= 0:
            raise ValueError("Seed coordinates must be finite and the radius and calibration must be positive.")
        circle = Point(x, y).buffer(radius, quad_segs=64)
        return cls.validate_geometry(scale(circle, xfact=1 / units_x, yfact=1 / units_y, origin=(x, y)))

    def seed_radius_geometry(self, x, y, radius, units_x=1., units_y=1.):
        transform = getattr(self.sample, 'spatial_transform', None)
        if transform is None:
            return self.radius_geometry(x, y, radius, units_x, units_y)
        seed = transform.points([x, y], inverse=True)
        circle = self.radius_geometry(*seed, radius, units_x, units_y)
        return self.validate_geometry(transform.geometry(circle))

    def radius_scale(self):
        data = None if self.coordinate_scale is None else self.coordinate_scale.data
        fx = (data['px_x'][0] if data else None) or 1 / getattr(self.sample, 'tile_scale_xfactor', 1.)
        fy = (data['px_y'][0] if data else None) or 1 / getattr(self.sample, 'tile_scale_yfactor', 1.)
        um = data['um_per_px'][0] if data else None
        if getattr(self.sample, 'spatial_transform', None) is not None:
            fx = fy = 1.
        return ((fx * um, fy * um), 'µm') if um is not None and um > 0 else ((fx, fy), 'full-res px')

    def compute_membership(self, geometry):
        """Intersect full native geometries, then align to sequenced observations."""
        with self.sample.fetch_lock:
            gdf = self.sample.gdf
            candidates = gdf.sindex.query(geometry, predicate='intersects')
            names = gdf.iloc[candidates][self.sample.on].astype(str)
            return np.asarray(self.sample.adata.obs_names.isin(names), dtype=bool)

    def compute_group_radius(self, group, radius, factors, mask=None, cancelled=None):
        """Union radius neighbourhoods of every native bin in a group (worker-safe)."""
        if not hasattr(self.sample, 'AGGREGATION_FACTORS'):
            raise ValueError("Group-bin seeds are available on bins layers.")
        self.radius_geometry(0., 0., radius, *factors)  # Validate radius and calibration.
        with self.sample.fetch_lock:
            key = f'grp_{group}'
            if key not in self.sample.adata.obs:
                raise ValueError("The seed group no longer exists.")
            if mask is None:
                mask = self.sample.adata.obs[key].fillna(False).to_numpy(dtype=bool)
            names = self.sample.adata.obs_names[np.asarray(mask, dtype=bool)]
            if not len(names):
                raise ValueError("The seed group contains no sequenced bins.")
            gdf = self.sample.gdf
            rows = gdf.loc[gdf[self.sample.on].astype(str).isin(names)].drop_duplicates(self.sample.on)
            centers = shapely.centroid(rows.geometry.to_numpy())
            coordinates = np.column_stack((shapely.get_x(centers), shapely.get_y(centers)))
            if len(rows) != len(names) or not np.isfinite(coordinates).all():
                raise ValueError("Every grouped bin must have valid geometry for its seed.")
        transform = getattr(self.sample, "spatial_transform", None)
        if transform is not None:
            coordinates = transform.points(coordinates, inverse=True)
        coordinates *= np.asarray(factors)
        levels = []
        # Bound temporary geometry memory; nearby circles merge within each chunk.
        for start in range(0, len(coordinates), 256):
            if cancelled is not None and cancelled():
                raise InterruptedError("Group preview cancelled.")
            points = shapely.points(coordinates[start:start + 256])
            region = shapely.union_all(shapely.buffer(points, radius, quad_segs=64))
            level = 0
            # Merge similarly sized batches instead of repeatedly unioning an
            # ever-growing region, especially costly for disjoint group bins.
            while level < len(levels) and levels[level] is not None:
                if cancelled is not None and cancelled():
                    raise InterruptedError("Group preview cancelled.")
                region = shapely.union_all([levels[level], region])
                levels[level] = None
                level += 1
            if level == len(levels):
                levels.append(region)
            else:
                levels[level] = region
            report(f"Buffered {min(start + 256, len(coordinates)):,} / {len(coordinates):,} group-bin seeds")
        region = shapely.union_all([region for region in levels if region is not None])
        region = scale(region, xfact=1 / factors[0], yfact=1 / factors[1], origin=(0, 0))
        if transform is not None:
            region = transform.geometry(region)
        return self.validate_geometry(region), len(names)

    def _group_mode(self):
        return (self._layout is not None and self.method.value == 'radius'
                and self.seed_source.value == 'group' and hasattr(self.sample, 'AGGREGATION_FACTORS'))

    def _groups_changed(self):
        self._refresh_list()
        if self._layout is not None:
            names = list(self.component.groups_manager._clusters)
            current = self.seed_group.value
            self.seed_group.options = [('', 'Choose a group')] + [(name, name) for name in names]
            if current and current not in names:
                self.seed_group.value = ''
            elif current and self._group_mode() and self._requested_seed_mask is not None:
                group, previous = self._requested_seed_mask
                values = self.sample.adata.obs.get(f'grp_{current}')
                if (group == current and values is not None
                        and not np.array_equal(previous, values.fillna(False).to_numpy(dtype=bool))):
                    self.preview_group_radius()

    def preview_group_radius(self, immediate=False):
        """Queue the latest group/radius preview; obsolete geometry work is cancelled."""
        if not self._group_mode() or self.disposed or not self.visible:
            return
        self.cancel()
        group = self.seed_group.value
        if not group:
            self.seed_hint.text = ''
            self.logger.set_message("Choose a seed group to preview one combined ROI.")
            return
        try:
            factors, unit = self.radius_scale()
            radius = self.radius.value
            self.radius_geometry(0., 0., radius, *factors)
            key = f'grp_{group}'
            if key not in self.sample.adata.obs:
                raise ValueError("The seed group no longer exists.")
            mask = self.sample.adata.obs[key].fillna(False).to_numpy(dtype=bool, copy=True)
            if not mask.any():
                raise ValueError("The seed group contains no sequenced bins.")
        except (ValueError, TypeError) as exc:
            self._error(exc)
            return
        self.radius.title = f'Radius ({unit})'
        self._requested_seed_mask = (group, mask.copy())
        self._pending_preview = (self._preview_generation, group, radius, factors, mask,
                                 self.sample.adata, self.sample.adata.obs_names.copy())
        self.seed_hint.text = f"All {int(mask.sum()):,} native bins are seeds; their neighbourhoods form one ROI."
        if self.logger.running:
            return
        if immediate:
            self._start_pending_preview()
        else:
            self._preview_timeout = self.doc.add_timeout_callback(self._start_pending_preview, 200)

    def _start_pending_preview(self):
        self._preview_timeout = None
        if self.logger.running or self._pending_preview is None or self.disposed or not self._group_mode():
            return
        generation, group, radius, factors, mask, adata, names = self._pending_preview
        self._pending_preview = None
        cancel = self._preview_cancel = threading.Event()

        def done(result):
            if generation != self._preview_generation or self.disposed:
                if not self.disposed:
                    self.logger.set_message("Preview cancelled.")
                    self._start_pending_preview()
                return
            key = f'grp_{group}'
            if (self.sample.adata is not adata or not names.equals(adata.obs_names)
                    or key not in adata.obs or not np.array_equal(mask, adata.obs[key].fillna(False).to_numpy(dtype=bool))):
                self.preview_group_radius(immediate=True)
                return
            geometry, count = result
            self.set_draft(geometry, 'radius')
            self._draft_seed_info = (group, count)
            self.logger.success(f"Preview ready: one ROI from all {count:,} group-bin seeds. Save ROI to keep it.")

        def failed(exc):
            if self.disposed:
                return
            if generation == self._preview_generation:
                self._error(exc)
            else:
                self.logger.set_message("Preview cancelled.")
            self._start_pending_preview()

        self.logger.run('Group-radius preview',
                        lambda: self.compute_group_radius(group, radius, factors, mask, cancel.is_set), done, failed)

    def set_draft(self, geometry, kind):
        self.draft = self.validate_geometry(geometry)
        self.draft_kind = kind
        self._draft_seed_info = None
        self._refresh_overlay()
        if self._layout is not None:
            self.save_button.disabled = self.logger.running
            self.logger.set_message("Draft ready. Save ROI or cancel.")

    def begin(self):
        if self._group_mode():
            self.preview_group_radius(immediate=True)
            return
        if self.disposed or not self.visible or self.logger.running:
            return
        self.cancel()
        self.radius.title = f"Radius ({self.radius_scale()[1]})"
        if self.method.value == 'radius' and (self.radius.value is None or not np.isfinite(self.radius.value) or self.radius.value <= 0):
            self._error("Enter a positive radius.")
            return
        self.drawing.arm(self, self.method.value)
        self.logger.set_message(f"Drawing for {escape(self.component.name)}: " +
                                ("click the seed on the map." if self.method.value == 'radius' else "drag the region on the map."))

    def cancel(self):
        self._preview_generation += 1
        self._pending_preview = None
        if self._preview_cancel is not None:
            self._preview_cancel.set()
        if self._preview_timeout is not None:
            try:
                self.doc.remove_timeout_callback(self._preview_timeout)
            except ValueError:
                pass
            self._preview_timeout = None
        self.drawing.release(self)
        self.draft = None
        self.draft_kind = None
        self._draft_seed_info = None
        self._refresh_overlay()
        if self._layout is not None:
            self.save_button.disabled = True
            if not self.logger.running:
                self.logger.set_message("Choose a method and draw an ROI.")

    def add_roi(self, roi):
        self.validate_geometry(roi.geometry)
        if roi.mask is not None:
            roi.mask = np.asarray(roi.mask, dtype=bool).copy()
            if roi.mask.shape != (self.sample.adata.n_obs,):
                raise ValueError("ROI membership must align with the native observations.")
        if any(existing.name == roi.name or existing.id == roi.id for existing in self.rois):
            raise ValueError("An ROI with that name already exists in this layer.")
        self.rois.append(roi)
        self._refresh_overlay()
        self._refresh_list()

    def save_roi(self):
        if self.draft is None or self.logger.running or not self.visible:
            return
        name = self.name_input.value.strip()
        if not name:
            while any(roi.name == f"ROI {self._counter}" for roi in self.rois):
                self._counter += 1
            name = f"ROI {self._counter}"
        roi = ROI(name, self.draft, self.draft_kind, palette_color(GROUP_PALETTE, self._counter - 1))
        if self._draft_seed_info:
            roi.seed_group, roi.seed_count = self._draft_seed_info
        try:
            self.add_roi(roi)
        except ValueError as exc:
            self._error(exc)
            return
        self._counter += 1
        self.cancel()
        self.name_input.value = ""
        adata = self.sample.adata
        obs_names = adata.obs_names.copy()

        def done(mask):
            if self.disposed or roi not in self.rois:
                if not self.disposed:
                    self.logger.set_message("ROI deleted.")
                    self._start_pending_preview()
                return
            if self.sample.adata is not adata or not obs_names.equals(self.sample.adata.obs_names):
                failed(ValueError("Layer data changed during ROI computation; delete the ROI and draw it again."))
                return
            roi.mask = mask
            self.logger.success(f"{escape(roi.name)}: {roi.count:,} observations")
            self._refresh_list()
            self._start_pending_preview()

        def failed(exc):
            if self.disposed or roi not in self.rois:
                if not self.disposed:
                    self.logger.set_message("ROI deleted.")
                    self._start_pending_preview()
                return
            roi.error = str(exc)
            self._error(exc)
            self._refresh_list()
            self._start_pending_preview()

        self.logger.run('ROI membership', lambda: self.compute_membership(roi.geometry), done, failed)

    def get_rois(self):
        return list(self.rois)

    def focus(self, roi):
        x0, y0, x1, y1 = roi.geometry.bounds
        px, py = .05 * (x1 - x0), .05 * (y1 - y0)
        self.component.figure.x_range.update(start=x0 - px, end=x1 + px)
        self.component.figure.y_range.update(start=y0 - py, end=y1 + py)

    def select(self, roi):
        if roi.mask is not None and roi.mask.any():
            self.component.set_selected_obs_names(self.sample.adata.obs_names[roi.mask])

    def create_group(self, roi):
        if roi.mask is None or not roi.mask.any():
            return
        name = self.component.groups_manager.add_mask_group(f"ROI: {roi.name}", roi.mask, self.doc, roi.id)
        self.logger.set_message(f"Created group: {escape(name)}")
        self._refresh_list()
        return name

    def remove_roi(self, roi):
        if isinstance(roi, str):
            roi = next((item for item in self.rois if item.id == roi), None)
        if roi is None or not any(item is roi for item in self.rois):
            return
        groups = self.component.groups_manager
        linked = groups.group_for_roi(roi.id)
        if linked:
            groups.delete_group(linked, self.doc)
        self.rois = [item for item in self.rois if item is not roi]
        self._refresh_overlay()
        self._refresh_list()

    def set_visible(self, visible):
        self.visible = visible
        self.renderer.visible = visible
        if not visible:
            self.cancel()

    def dispose(self):
        self.cancel()
        self.disposed = True
        self.component.groups_manager.remove_change_listener(self._groups_changed)
        if self.renderer in self.component.figure.renderers:
            self.component.figure.renderers.remove(self.renderer)

    def _error(self, exc):
        if self.logger is not None:
            self.logger.error(exc)

    def _refresh_overlay(self):
        shapes = [(roi.geometry, roi.color) for roi in self.rois if roi.visible]
        if self.draft is not None:
            shapes.append((self.draft, '#ffffff'))
        polygons = [[shape] if shape.geom_type == 'Polygon' else list(shape.geoms) for shape, _ in shapes]
        self.source.data = dict(xs=[[[list(ring.coords.xy[0]) for ring in [poly.exterior, *poly.interiors]]
                                    for poly in parts] for parts in polygons],
                                ys=[[[list(ring.coords.xy[1]) for ring in [poly.exterior, *poly.interiors]]
                                    for poly in parts] for parts in polygons],
                                color=[color for _, color in shapes])

    def _refresh_list(self):
        if self._layout is None or self.disposed:
            return
        children = []
        for roi in self.rois:
            count = 'Computing membership…' if roi.mask is None else f'{roi.count:,} observations'
            if roi.error:
                count = f'Error: {escape(roi.error)}'
            title = Div(text=f"<b style='color:{roi.color}'>{escape(roi.name)}</b> · {roi.kind}<br>{count}")
            if roi.seed_group is not None:
                title.text += f"<br>{roi.seed_count:,} seeds from {escape(roi.seed_group)}"
            focus = Button(label='Focus', width=65)
            select = Button(label='Select', width=65, disabled=not roi.count)
            group = Button(label='Create group', width=110,
                           disabled=not roi.count or bool(self.component.groups_manager.group_for_roi(roi.id)))
            delete = Button(label='Delete', width=65, button_type='danger')
            visible = Checkbox(label='Outline', active=roi.visible)
            focus.on_click(lambda item=roi: self.focus(item))
            select.on_click(lambda item=roi: self.select(item))
            group.on_click(lambda item=roi: self._group_click(item))
            delete.on_click(lambda item=roi: self.remove_roi(item))
            visible.on_change('active', lambda attr, old, new, item=roi: self._toggle_outline(item, new))
            children.append(column(title, row(focus, select, visible), row(group, delete), sizing_mode='stretch_width'))
        self.roi_list.children = children or [Div(text='No ROIs saved in this layer.')]

    def _toggle_outline(self, roi, visible):
        roi.visible = visible
        self._refresh_overlay()

    def _group_click(self, roi):
        try:
            self.create_group(roi)
        except Exception as exc:
            self._error(exc)

    def get_widgets(self, doc):
        self.doc = doc
        if self._layout is None:
            self.logger = LoadingLogger(doc)
            self.name_input = TextInput(title='ROI name (optional)', placeholder='ROI 1')
            self.method = Select(title='Geometry', value='box', options=[('box', 'Bounding box'), ('lasso', 'Lasso'), ('radius', 'Seed + radius')])
            self.radius = NumericInput(title=f'Radius ({self.radius_scale()[1]})', value=50., low=0, visible=False)
            self.seed_source = Select(title='Seed source', value='map',
                                      options=[('map', 'Click on map'), ('group', 'All bins in a group')], visible=False)
            self.seed_group = Select(title='Seed group', options=[('', 'Choose a group')], value='', visible=False)
            self.seed_hint = Div(text='', visible=False)
            def method_changed(attr, old, new):
                self.cancel()
                self._update_seed_controls()
                self.preview_group_radius()
            self.method.on_change('value', method_changed)
            self.seed_source.on_change('value', method_changed)
            self.seed_group.on_change('value', lambda attr, old, new: self.preview_group_radius())
            self.radius.on_change('value', lambda attr, old, new: self.preview_group_radius())
            self.draw_button = Button(label='Draw', button_type='primary', width=100)
            self.save_button = Button(label='Save ROI', button_type='success', width=90, disabled=True)
            cancel = Button(label='Cancel', width=70)
            self.draw_button.on_click(self.begin)
            self.save_button.on_click(self.save_roi)
            cancel.on_click(self.cancel)
            self.roi_list = column(sizing_mode='stretch_width')
            self._layout = column(self.name_input, self.method, self.radius, self.seed_source, self.seed_group,
                                  self.seed_hint, row(self.draw_button, self.save_button, cancel),
                                  self.logger.div, self.roi_list, sizing_mode='stretch_width')
            self._refresh_list()
            self._groups_changed()
            self.logger.set_message('Choose a method and draw an ROI.')
        return [self._layout]

    def _update_seed_controls(self):
        self.radius.visible = self.method.value == 'radius'
        self.seed_source.visible = self.radius.visible and hasattr(self.sample, 'AGGREGATION_FACTORS')
        self.seed_group.visible = self._group_mode()
        self.seed_hint.visible = self._group_mode()
        self.draw_button.label = 'Preview group' if self._group_mode() else 'Draw'
