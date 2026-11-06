
from dataclasses import dataclass, field

from html import escape
import uuid

import threading

import numpy as np
import pandas as pd

import shapely
from shapely.geometry import Point, Polygon, box
from shapely.affinity import scale

from bokeh.events import SelectionGeometry, Tap
from bokeh.layouts import column, row
from bokeh.models import AutocompleteInput, BoxSelectTool, Button, Checkbox, ColumnDataSource, CustomJS, Div, LassoSelectTool, NumericInput, Select, TapTool, TextInput

from .... import settings
from ...utils.loading import LoadingLogger
from ...utils.palette import GROUP_PALETTE, palette_color
from ....startup.utils.progress import report
from ....core.roi import ROI, ROIGeometry, _ROIOverlay  # noqa: F401  (re-exported)


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


class ROIManager(ROIGeometry):
    def __init__(self, component):
        self.component = component
        super().__init__()
        self.doc = None
        self.logger = None
        self.draft = None
        self.draft_kind = None
        self.visible = True
        self.disposed = False

        self._layout = None
        self._preview_generation = 0
        self._preview_cancel = None
        self._preview_timeout = None
        self._pending_preview = None
        self._draft_seed_info = None
        self._draft_excluded_seed_mask = None
        self._requested_seed_mask = None
        self._draft_seed_barcode = None

        figure = component.figure
        if not hasattr(figure, '_roi_drawing_controller'):
            figure._roi_drawing_controller = _DrawingController(figure)

        self.drawing = figure._roi_drawing_controller

        self.source = ColumnDataSource(dict(xs=[], ys=[], color=[]))

        self.renderer = figure.multi_polygons('xs', 'ys', source=self.source, fill_alpha=.06,
                                       line_color='color', fill_color='color', line_width=2,
                                       selection_fill_alpha=.06, nonselection_fill_alpha=.06,
                                       level='overlay', name=f'{component.key}_rois')

        self.raster_source = ColumnDataSource(dict(image=[], x=[], y=[], dw=[], dh=[]))
        self.raster_renderer = figure.image_rgba(
            image='image', x='x', y='y', dw='dw', dh='dh', source=self.raster_source,
            level='overlay', name=f'{component.key}_roi_rasters')

        component.groups_manager.add_change_listener(self._groups_changed)

    @property
    def sample(self):
        return self.component.sample









    def _group_mode(self):
        return (self._layout is not None and self.method.value == 'radius'
                and self.seed_source.value == 'group' and hasattr(self.sample, 'AGGREGATION_FACTORS'))

    def _barcode_mode(self):
        return (self._layout is not None and self.method.value == 'radius'
                and self.seed_source.value == 'barcode' and hasattr(self.sample, 'AGGREGATION_FACTORS'))




    def _barcode_input_changed(self, attr, old, new):
        if not self._barcode_mode() or self.disposed:
            return
        if self._draft_seed_barcode is not None and new.strip() != self._draft_seed_barcode:
            self.cancel()
            self.seed_hint.text = 'Choose an exact barcode to preview the radius ROI.'
        self.seed_barcode.completions = self._barcode_suggestions(new)

    def preview_barcode_radius(self, focus=False):
        if not self._barcode_mode() or self.disposed or not self.visible or self.logger.running:
            return
        self.cancel()
        barcode = self.seed_barcode.value.strip()
        if not barcode:
            self.seed_hint.text = 'Type a native barcode or bin index and choose a suggestion.'
            self.logger.set_message('Choose a barcode to preview its radius ROI.')
            return
        try:
            x, y = self.barcode_seed_center(barcode)
            factors, unit = self.radius_scale()
            geometry = self.seed_radius_geometry(x, y, self.radius.value, *factors)
            self.radius.title = f'Radius ({unit})'
            self.set_draft(geometry, 'radius')
        except (ValueError, TypeError, AttributeError) as exc:
            self.seed_hint.text = ''
            self._error(exc)
            return
        self._draft_seed_barcode = barcode
        self.seed_hint.text = f'Seed: {escape(barcode)} · native bin centre ({x:.2f}, {y:.2f}).'
        self.logger.set_message('Barcode preview ready. Save ROI to keep it.')
        if focus:
            self._focus_geometry(geometry)

    def _preview_radius(self):
        if self._barcode_mode():
            self.preview_barcode_radius()
        else:
            self.preview_group_radius()

    def _groups_changed(self):
        self._refresh_list()
        if self._layout is not None:
            names = list(self.component.groups_manager._groups)
            current = self.seed_group.value
            self.seed_group.options = [('', 'Choose a group')] + [(name, name) for name in names]
            if current and current not in names:
                self.seed_group.value = ''
            elif current and self._group_mode() and self._requested_seed_mask is not None:
                group, previous = self._requested_seed_mask
                values = self.sample.adata.obs.get(f'{settings.GROUP_MASK_PREFIX}{current}')
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
            key = f'{settings.GROUP_MASK_PREFIX}{group}'
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
        exclude_group = self.group_region.value == 'minus_group'
        self._pending_preview = (self._preview_generation, group, radius, factors, mask,
                                 self.sample.adata, self.sample.adata.obs_names.copy(), exclude_group)
        self.seed_hint.text = f"All {int(mask.sum()):,} native bins are seeds; their neighbourhoods form one ROI."
        if exclude_group:
            self.seed_hint.text += " Source-group bin shapes and members are excluded."
        if self.logger.running:
            return
        if immediate:
            self._start_pending_preview()
        else:
            self._preview_timeout = self.doc.add_timeout_callback(self._start_pending_preview, settings.ROI_PREVIEW_DEBOUNCE_MS)

    def _start_pending_preview(self):
        self._preview_timeout = None
        if self.logger.running or self._pending_preview is None or self.disposed or not self._group_mode():
            return
        generation, group, radius, factors, mask, adata, names, exclude_group = self._pending_preview
        self._pending_preview = None
        cancel = self._preview_cancel = threading.Event()

        def done(result):
            if generation != self._preview_generation or self.disposed:
                if not self.disposed:
                    self.logger.set_message("Preview cancelled.")
                    self._start_pending_preview()
                return
            key = f'{settings.GROUP_MASK_PREFIX}{group}'
            if (self.sample.adata is not adata or not names.equals(adata.obs_names)
                    or key not in adata.obs or not np.array_equal(mask, adata.obs[key].fillna(False).to_numpy(dtype=bool))):
                self.preview_group_radius(immediate=True)
                return
            geometry, count, overlay = result
            self.set_draft(geometry, 'radius', _overlay=overlay)
            self._draft_seed_info = (group, count)
            self._draft_excluded_seed_mask = mask.copy() if exclude_group else None
            detail = " Source group excluded." if exclude_group else ""
            if overlay.raster is not None:
                detail += " Complex outline shown as a bitmap preview; membership uses exact geometry."
            self.logger.success(f"Preview ready: one ROI from all {count:,} group-bin seeds.{detail} Save ROI to keep it.")

        def failed(exc):
            if self.disposed:
                return
            if generation == self._preview_generation:
                self._error(exc)
            else:
                self.logger.set_message("Preview cancelled.")
            self._start_pending_preview()

        def work():
            geometry, count = self.compute_group_radius(group, radius, factors, mask, cancel.is_set, exclude_group)
            report("Preparing ROI display")
            overlay = self._prepare_overlay(geometry, cancel.is_set)
            return geometry, count, overlay

        self.logger.run('Group-radius preview', work, done, failed)

    def set_draft(self, geometry, kind, *, _overlay=None):
        if _overlay is None:
            geometry = self.validate_geometry(geometry)
            _overlay = self._prepare_overlay(geometry)
        self.draft = geometry
        self._overlay_cache[id(geometry)] = _overlay
        self.draft_kind = kind
        self._draft_seed_info = None
        self._draft_excluded_seed_mask = None
        self._draft_seed_barcode = None
        self._refresh_overlay()
        if self._layout is not None:
            self.save_button.disabled = self.logger.running
            self.logger.set_message("Draft ready. Save ROI or cancel.")

    def begin(self):
        if self._barcode_mode():
            self.preview_barcode_radius(focus=True)
            return
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
        self._draft_excluded_seed_mask = None
        self._draft_seed_barcode = None
        self._refresh_overlay()
        if self._layout is not None:
            self.save_button.disabled = True
            if not self.logger.running:
                self.logger.set_message("Choose a method and draw an ROI.")

    def add_roi(self, roi):
        if id(roi.geometry) not in self._overlay_cache:
            self.validate_geometry(roi.geometry)
        if roi.mask is not None:
            roi.mask = np.asarray(roi.mask, dtype=bool).copy()
            if roi.mask.shape != (self.sample.adata.n_obs,):
                raise ValueError("ROI membership must align with the native observations.")
        if roi.excluded_seed_mask is not None:
            roi.excluded_seed_mask = np.asarray(roi.excluded_seed_mask, dtype=bool).copy()
            if roi.excluded_seed_mask.shape != (self.sample.adata.n_obs,):
                raise ValueError("Excluded group membership must align with the native observations.")
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
        roi.seed_barcode = self._draft_seed_barcode
        if self._draft_seed_info:
            roi.seed_group, roi.seed_count = self._draft_seed_info
            if self._draft_excluded_seed_mask is not None:
                roi.excluded_seed_mask = self._draft_excluded_seed_mask.copy()
        
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

        def membership():
            mask = self.compute_membership(roi.geometry)
            # Intersects includes boundary contact. Remove the captured source bins
            # explicitly even when they touch the edge of a subtracted polygon.
            return mask if roi.excluded_seed_mask is None else mask & ~roi.excluded_seed_mask

        self.logger.run('ROI membership', membership, done, failed)

    def get_rois(self):
        return list(self.rois)

    def focus(self, roi):
        self._focus_geometry(roi.geometry)

    def _focus_geometry(self, geometry):
        x0, y0, x1, y1 = geometry.bounds
        px, py = settings.ROI_FOCUS_PADDING * (x1 - x0), settings.ROI_FOCUS_PADDING * (y1 - y0)
        self.component.figure.x_range.update(start=x0 - px, end=x1 + px)
        self.component.figure.y_range.update(start=y0 - py, end=y1 + py)

    def select(self, roi):
        if roi.mask is not None and roi.mask.any():
            self.component.set_selected_obs_names(self.sample.adata.obs_names[roi.mask])

    def unselect(self, roi):
        if roi.mask is not None:
            members = set(self.sample.adata.obs_names[roi.mask])
            self.component.set_selected_obs_names(self.component.get_selected_obs_names() - members)

    def create_group(self, roi, prepared=None):
        if roi.mask is None or not roi.mask.any():
            return

        name = self.component.groups_manager.add_mask_group(f"ROI: {roi.name}", roi.mask, self.doc, roi.id, prepared=prepared)
        self.logger.set_message(f"Created group: {escape(name)}")
        self._refresh_list()
        return name

    def remove_roi(self, roi, prepared=None):
        if isinstance(roi, str):
            roi = next((item for item in self.rois if item.id == roi), None)
        if roi is None or not any(item is roi for item in self.rois):
            return
        groups = self.component.groups_manager
        linked = groups.group_for_roi(roi.id)
        if linked:
            groups.delete_group(linked, self.doc, prepared=prepared)
        self.unselect(roi)
        self.rois = [item for item in self.rois if item is not roi]
        self._refresh_overlay()
        self._refresh_list()

    def set_visible(self, visible):
        self.visible = visible
        self.renderer.visible = visible
        self.raster_renderer.visible = visible
        if not visible:
            self.cancel()

    def dispose(self):
        self.cancel()
        self.disposed = True
        self.component.groups_manager.remove_change_listener(self._groups_changed)
        for renderer in (self.renderer, self.raster_renderer):
            if renderer in self.component.figure.renderers:
                self.component.figure.renderers.remove(renderer)
        self._overlay_cache.clear()
        self._group_union_cache = None
        self._barcode_cache = None

    def _error(self, exc):
        if self.logger is not None:
            self.logger.error(exc)



    def _refresh_overlay(self):
        shapes = [(roi.geometry, roi.color) for roi in self.rois if roi.visible]
        if self.draft is not None:
            shapes.append((self.draft, '#ffffff'))
        keep = {id(roi.geometry) for roi in self.rois}
        if self.draft is not None:
            keep.add(id(self.draft))
        self._overlay_cache = {key: value for key, value in self._overlay_cache.items() if key in keep}
        vectors = dict(xs=[], ys=[], color=[])
        rasters = dict(image=[], x=[], y=[], dw=[], dh=[])
        for geometry, color in shapes:
            overlay = self._overlay_cache.get(id(geometry))
            if overlay is None:
                overlay = self._prepare_overlay(geometry)
                self._overlay_cache[id(geometry)] = overlay
            if overlay.raster is None:
                vectors['xs'].append(overlay.xs)
                vectors['ys'].append(overlay.ys)
                vectors['color'].append(color)
            else:
                x0, y0, x1, y1 = overlay.bounds
                rasters['image'].append(overlay.image(color))
                rasters['x'].append(x0)
                rasters['y'].append(y0)
                rasters['dw'].append(x1 - x0)
                rasters['dh'].append(y1 - y0)
        self.source.data = vectors
        self.raster_source.data = rasters

    def _refresh_list(self):
        if self._layout is None or self.disposed:
            return

        children = []
        for roi in self.rois:
            count = 'Computing membership…' if roi.mask is None else f'{roi.count:,} observations'
            if roi.error:
                count = f'Error: {escape(roi.error)}'
            
            title = Div(text=f"<b style='color:{roi.color}'>{escape(roi.name)}</b> · {roi.kind}<br>{count}")
            overlay = self._overlay_cache.get(id(roi.geometry))
            if overlay is not None and overlay.raster is not None:
                title.text += "<br>Overview image · exact membership"
            if roi.seed_group is not None:
                title.text += f"<br>{roi.seed_count:,} seeds from {escape(roi.seed_group)}"
                if roi.excluded_seed_mask is not None:
                    title.text += " · source group excluded"
            if roi.seed_barcode is not None:
                title.text += f'<br>Seed: {escape(roi.seed_barcode)}'
            
            focus = Button(label='Focus', width=65)
            select = Button(label='Select', width=65, disabled=not roi.count)
            unselect = Button(label='Unselect', width=75, disabled=not roi.count)
            group = Button(label='Create group', width=110,
                           disabled=not roi.count or bool(self.component.groups_manager.group_for_roi(roi.id)))
            delete = Button(label='Delete', width=65, button_type='danger')
            visible = Checkbox(label='Outline', active=roi.visible)
            
            focus.on_click(lambda item=roi: self.focus(item))
            select.on_click(lambda item=roi: self.select(item))
            unselect.on_click(lambda item=roi: self.unselect(item))
            group.on_click(lambda item=roi: self._group_click(item))
            delete.on_click(lambda item=roi: self._delete_click(item))
            visible.on_change('active', lambda attr, old, new, item=roi: self._toggle_outline(item, new))
            
            children.append(column(title, row(focus, select, unselect), row(group, delete, visible), sizing_mode='stretch_width'))
        
        self.roi_list.children = children or [Div(text='No ROIs saved in this layer.')]

    def _toggle_outline(self, roi, visible):
        roi.visible = visible
        self._refresh_overlay()

    def _group_click(self, roi):
        if self.logger.running or roi.mask is None or not roi.mask.any():
            return
        groups = self.component.groups_manager
        if groups.group_for_roi(roi.id):
            return
        # Reserve the same unique name used by add_mask_group.
        base = name = f"ROI: {roi.name}"
        suffix = 2
        while name in groups._groups:
            name = f"{base} ({suffix})"
            suffix += 1
        definitions = {**groups._groups, name: []}
        fixed = dict(groups._fixed_masks)
        fixed[name] = pd.Series(roi.mask.copy(), index=self.sample.adata.obs_names.copy())

        def done(prepared):
            if any(item is roi for item in self.rois):
                self.create_group(roi, prepared=prepared)

        groups.run_group_action(self.doc, f"Creating group '{name}'",
                                lambda: groups._prepare_groups(definitions, fixed), done,
                                f"Group '{escape(name)}' created.", logger=self.logger)

    def _delete_click(self, roi):
        if self.logger.running:
            return
        groups = self.component.groups_manager
        linked = groups.group_for_roi(roi.id)
        if linked is None:
            self.remove_roi(roi)
            self.logger.set_message('ROI deleted.')
            return
        remaining = {name: conditions for name, conditions in groups._groups.items() if name != linked}
        groups.run_group_action(self.doc, f"Removing ROI and group '{linked}'",
                                lambda: groups._prepare_groups(remaining),
                                lambda prepared: self.remove_roi(roi, prepared=prepared),
                                f"ROI and group '{escape(linked)}' deleted.", logger=self.logger)

    def get_widgets(self, doc):
        self.doc = doc
        if self._layout is None:
            self.logger = LoadingLogger(doc)
            self.name_input = TextInput(title='ROI name (optional)', placeholder='ROI 1')
            self.method = Select(title='Geometry', value='box', options=[('box', 'Bounding box'), ('lasso', 'Lasso'), ('radius', 'Seed + radius')])
            self.radius = NumericInput(title=f'Radius ({self.radius_scale()[1]})', value=settings.ROI_DEFAULT_RADIUS, low=0, visible=False)
            self.seed_source = Select(title='Seed source', value='map',
                                      options=[('map', 'Click on map'), ('barcode', 'Find barcode / bin index'),
                                               ('group', 'All bins in a group')], visible=False)
            self.seed_barcode = AutocompleteInput(title='Barcode / bin index', placeholder='Type a native barcode…',
                                                   completions=[], min_characters=1,
                                                   max_completions=settings.ROI_BARCODE_MAX_COMPLETIONS,
                                                   case_sensitive=True, search_strategy='starts_with',
                                                   restrict=False, visible=False)
            # Suggestions arrive from the server after the original keyup.
            # Refresh the focused dropdown without committing the typed value.
            self.seed_barcode.js_on_change('completions', CustomJS(args=dict(barcode=self.seed_barcode), code="""
                const view = Bokeh.index.find_one(barcode)
                if (view?.input_el?.matches(':focus')) {
                    view.input_el.dispatchEvent(new KeyboardEvent('keyup', {key: 'Unidentified', bubbles: true}))
                }
            """))
            self.seed_group = Select(title='Seed group', options=[('', 'Choose a group')], value='', visible=False)
            self.group_region = Select(title='ROI region', value='full',
                                       options=[('full', 'Full ROI'), ('minus_group', 'ROI minus group')], visible=False)
            self.seed_hint = Div(text='', visible=False)
            def method_changed(attr, old, new):
                self.cancel()
                self._update_seed_controls()
                self._preview_radius()
            self.method.on_change('value', method_changed)
            self.seed_source.on_change('value', method_changed)
            self.seed_group.on_change('value', lambda attr, old, new: self.preview_group_radius())
            self.seed_barcode.on_change('value_input', self._barcode_input_changed)
            self.seed_barcode.on_change('value', lambda attr, old, new: self.preview_barcode_radius(focus=True))
            self.radius.on_change('value', lambda attr, old, new: self._preview_radius())
            self.group_region.on_change('value', lambda attr, old, new: self.preview_group_radius())
            self.draw_button = Button(label='Draw', button_type='primary', width=100)
            self.save_button = Button(label='Save ROI', button_type='success', width=90, disabled=True)
            cancel = Button(label='Cancel', width=70)
            self.draw_button.on_click(self.begin)
            self.save_button.on_click(self.save_roi)
            cancel.on_click(self.cancel)
            self.roi_list = column(sizing_mode='stretch_width')
            self._layout = column(self.name_input, self.method, self.radius, self.seed_source, self.seed_barcode, self.seed_group, self.group_region,
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
        self.group_region.visible = self._group_mode()
        self.seed_barcode.visible = self._barcode_mode()
        self.seed_hint.visible = self._group_mode() or self._barcode_mode()
        self.draw_button.label = ('Preview barcode' if self._barcode_mode() else
                                 'Preview group' if self._group_mode() else 'Draw')
