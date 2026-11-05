
from pathlib import Path
import threading
import uuid
import numpy as np

from bokeh.models import (
    Div,
    HoverTool, 
    CheckboxGroup,
    AutocompleteInput,
    Select,
    Slider,
    ColorPicker,
    Button,
    CustomJS,
    Tabs,
    TabPanel,
    ColumnDataSource, 
    DataTable, 
    TableColumn,
    AllIndices,
    BooleanFilter,
    CDSView
)
from bokeh.layouts import column
from bokeh.transform import transform

from .... import settings

from ...source.manager import SourceManager, DynamicSourceManager
from ...source.parameter import StaticParameter, Scaler, CategoricalColor
from ...source.lazyloader import LazyLoader

from ....help import help_header

from ...layout.layer import (
    ClusteringManager,
    DeconvolutionManager,
    GroupsManager,
    ROIManager,
    DEManager,
    MetadataComparisonManager,
    ParamForm
)

from .utils.selection import native_selection, selection_indices
from .utils.hover_names import HoverNames

NOT_SEQUENCED_VAL = settings.NOT_SEQUENCED_VAL


DOWNLOADABLE_DATA_PATH = settings.DOWNLOADABLE_DATA_PATH

class Layer:

    AVAILABLE_ATTRS = ['color', 'alpha', 'fill_alpha', 'line_width', 'line_color', 'line_alpha']
    TABLE_ATTRS = ['attr', 'is_static', 'colname', 'value', 'minv', 'maxv', 'not_sequenced_val']

    CLUSTER_COLOR_LABEL = settings.CLUSTER_COLOR_LABEL

    adata_path = None

    def __init__(
            self,
            figure,
            sample,
            lazy_loading=False,
            hover_color=None,
            alias='Spot',
            key=None,
            name=None,
            ):
        
        self.figure = figure
        # alias is the unit ("Cell"); key/name tell apart several layers of one
        # unit, e.g. one cells layer per segmentation method
        self.key = key or f"{alias.lower()}s"
        self.name = name or f"{alias}s"

        self.lazy_loading = lazy_loading

        if self.lazy_loading:
            self.manager = DynamicSourceManager(sample)
        else:
            self.manager = SourceManager(sample)
        self.manager.component = self
        # Lazy layers keep row names on the server; hover fetches them on demand
        self.hover_names = (HoverNames(lambda: self.manager.loaded_obs_names, lambda: self.manager.row_version)
                            if lazy_loading else None)

        self.hover_color = hover_color
        self.alias = alias

        self.glyph_renderer = None
        self.hover = None

        self.checkbox = None
        self.slider = None
        self.colorpicker = None
        self.autocomplete_input = None
        self.expression_scale_select = None

        self.clustering_manager = ClusteringManager(self)
        self.deconvolution_manager = DeconvolutionManager(self)
        self.groups_manager = GroupsManager(self)
        self.roi_manager = ROIManager(self)
        self.de_manager = DEManager(self)
        self.metadata_comparison_manager = MetadataComparisonManager(self)

        self.selected_obs_names = set()
        self._selection_listeners = []
        self._columns_listeners = []      # called when obs columns change
        self._clustering_listeners = []   # called when a clustering finishes
        self._suspend_selection_sync = False
        self._selection_callback_bound = False
        
        self.sample_type = f"{self.alias.lower()}s"

        # Exports share a public static folder: unguessable per-export names keep
        # one user's file from being overwritten or fetched by another.
        self._export_name = None
        self._view = None

    @property
    def sample(self):
        return self.manager.sample

    @property
    def static_adata_path(self):
        if self._export_name is None:
            return None
        return (DOWNLOADABLE_DATA_PATH / self._export_name).resolve()

    @property
    def static_adata_url(self):
        return f"/hd_viewer/static/data/{self._export_name}"

    @property
    def widgets_name(self):
        return f"{self.key}_widgets"
    
    @property
    def param_table_name(self):
        return f"{self.key}_params"
    
    @property
    def view(self):
        if self._view is None or len(self._view) == 0:
            self._view = np.ones(self.manager.n_geometries, dtype=bool)
        return self._view

    def make_view(self):
        """CDSView for a new glyph: no filter at all while every row is shown."""
        booleans = np.asarray(self.view, dtype=bool)
        return CDSView(filter=AllIndices() if booleans.all() else BooleanFilter(booleans=booleans))

    def glyph_specs(self):
        """Glyph property → source column; colour-coded columns go through their
        parameter's palette."""
        return {param.attr: transform(param.attr, param.mapper) if getattr(param, 'mapper', None) is not None
                else param.attr for param in self.manager.parameters}

    def apply_glyph_constants(self, constants):
        """Per-level constants (bin width/height/angle) for glyphs that use them."""
    
    @staticmethod
    def _remove_export(path):
        try:
            Path(path).unlink(missing_ok=True)
        except OSError:
            pass

    def _cleanup_download_file(self):
        path = self.adata_path or self.static_adata_path
        if path is None:
            return
        self._remove_export(path)
        self.adata_path = None

    def save_adata(self):
        if self.lazy_loading:
            return None

        DOWNLOADABLE_DATA_PATH.mkdir(parents=True, exist_ok=True)
        self._export_name = f"{self.key}-{uuid.uuid4().hex}.h5ad"
        adata_path = self.static_adata_path
        try:
            with self.sample.fetch_lock:
                self.sample.adata.write_h5ad(adata_path)
        except Exception:
            if adata_path.exists():
                adata_path.unlink(missing_ok=True)
            self.adata_path = None
            raise

        self.adata_path = adata_path
        return adata_path

    def attach_loader(self, doc):
        if self.lazy_loading and self.glyph_renderer is not None:
            loader = LazyLoader(self.manager, self.figure)
            loader.add_range_callback(self, doc)
            self.loader = loader

            def _initial_load():
                loader.force_reload()
            doc.add_next_tick_callback(_initial_load)

    def refresh_layout(self, doc):
        widgets_layout = doc.get_model_by_name(self.widgets_name)
        widgets_layout.children = [self.get_widgets(doc)]

    def remove_existing_glyph(self):
        if self.glyph_renderer:
            self.figure.renderers.remove(self.glyph_renderer)
            self.glyph_renderer = None

    def refresh_glyph(self):
        self.remove_existing_glyph()
        self.add_glyph()

    def remove_existing_hover(self):
        if self.hover:
            self.figure.tools.remove(self.hover)
            self.hover = None

    def refresh_hover(self):
        # Reuse the tool: adding/removing tools makes BokehJS rebuild the toolbar,
        # which could leave the main figure without one.
        self.add_hover()

    def refresh_figure(self):
        self.refresh_glyph()
        self.refresh_hover()
        self._sync_color_picker()

    def _sync_color_picker(self):
        if self.colorpicker is None:
            return
        color = self.manager.get_param('color')
        hidden = self.glyph_renderer is not None and not self.glyph_renderer.visible
        self.colorpicker.disabled = hidden or color is None or not color.is_static

    @property
    def obs_name_values(self):
        return np.asarray(self.manager.loaded_obs_names, dtype=str)

    def get_selected_obs_names(self):
        return set(self.selected_obs_names)

    def selection_indices_for(self, obs_names):
        return selection_indices(self.sample, self.selected_obs_names, obs_names)

    def add_columns_listener(self, callback):
        if callback not in self._columns_listeners:
            self._columns_listeners.append(callback)

    def add_clustering_listener(self, callback):
        if callback not in self._clustering_listeners:
            self._clustering_listeners.append(callback)

    def add_selection_listener(self, callback):
        if callback not in self._selection_listeners:
            self._selection_listeners.append(callback)

    def notify_selection_changed(self):
        for callback in self._selection_listeners:
            callback(self)

    def apply_selected_obs_names(self):
        if self.glyph_renderer is None:
            return

        obs_names = self.obs_name_values
        if obs_names.size == 0 or len(self.selected_obs_names) == 0:
            indices = []
        else:
            indices = self.selection_indices_for(obs_names)
            
        self._suspend_selection_sync = True
        self.glyph_renderer.data_source.selected.indices = indices
        self._suspend_selection_sync = False

    def set_selected_obs_names(self, obs_names, notify=True):
        self.selected_obs_names = native_selection(self.sample, obs_names)
        self.apply_selected_obs_names()
        if notify:
            self.notify_selection_changed()

    def _on_source_selection_change(self, attr, old, new):
        if self._suspend_selection_sync:
            return
        obs_names = self.obs_name_values
        selected_indices = np.asarray(new, dtype=int)
        valid_indices = selected_indices[(selected_indices >= 0) & (selected_indices < obs_names.size)]
        self.selected_obs_names = native_selection(self.sample, obs_names[valid_indices])
        self.notify_selection_changed()

    def attach_selection_callbacks(self):
        if self.glyph_renderer is None:
            return
        selected = self.glyph_renderer.data_source.selected
        if not self._selection_callback_bound:
            selected.on_change("indices", self._on_source_selection_change)
            self._selection_callback_bound = True

    def add_glyph(self):
        pass

    def remove_parameter(self, attr):
        self.manager.remove_parameter(attr)
        self.refresh_figure()
        self.update_hover_tooltips()

    def add_parameter(
            self, 
            attr='alpha', 
            static=False,
            val=0.5,
            min_val=0.0,
            max_val=1.0,
            colname=None,
            static_not_sequenced=True,
            not_sequenced_val=0.0,
            colname_hovered=False
            ):
        
        if static:
            self.manager.add_parameter(
                StaticParameter(
                    val, 
                    attr=attr, 
                    not_sequenced_val=not_sequenced_val
                    )
            )
        else:
            assert attr in self.AVAILABLE_ATTRS, "Incorrect attribute"
            
            if "color" in attr:
                param = CategoricalColor(
                    colname,
                    attr=attr,
                    static_not_sequenced=static_not_sequenced,
                    not_sequenced_val=not_sequenced_val,
                    colname_hovered=colname_hovered
                )
                # Clusters and groups have assigned colours: use them before any
                # data is drawn, so the map, the DE table and the DE plots agree
                param.category_colors = self.known_category_colors(colname)
                self.manager.add_parameter(param)
            else:
                self.manager.add_parameter(
                    Scaler(
                        min_val,
                        max_val,
                        colname,
                        attr=attr,
                        static_not_sequenced=static_not_sequenced,
                        not_sequenced_val=not_sequenced_val,
                        colname_hovered=colname_hovered
                        )
                )
        self.refresh_figure()

    def known_category_colors(self, colname):
        """Assigned colours of the categories of `colname`, if it holds clusters or groups."""
        clustering = getattr(self, "clustering_manager", None)
        if clustering is not None and colname == clustering.cluster_col and clustering._color_map:
            return dict(clustering._color_map)

        groups = getattr(self, "groups_manager", None)
        if groups is not None and colname == settings.GROUPS_COL and groups._color_map:
            return dict(groups._color_map)
        
        return {}

    def get_hover_tooltips(self):
        items = []
        cluster_col = getattr(getattr(self, "clustering_manager", None), "cluster_col", None)
        for param in self.manager.parameters:
            if (not param.is_static and param.colname_hovered
                    and param.colname not in (settings.GROUPS_COL, cluster_col)
                    and param.attr != 'alpha'):
                label = param.colname
                if self.sample.expression_scale != "counts" and param.colname in self.sample.adata.var_names:
                    label += f" ({'log-normalized' if self.sample.expression_scale == 'log1p' else 'normalized'})"
                items.append(f"{label}: @{param.colname}")

        if settings.GROUPS_COL in self.manager.source.data:
            items.append(f"Group: @{{{settings.GROUPS_COL}}}")
        if cluster_col is not None and cluster_col in self.manager.source.data:
            items.append(f"Cluster: @{{{cluster_col}}}")

        if len(items) > 0:
            return ", ".join(items)
        return HoverNames.TEMPLATE if self.hover_names is not None else "@obs_name"

    def add_hover(self):
        if not self.glyph_renderer:
            return
        if self.hover is None:
            self.hover = HoverTool(
                tooltips=self.get_hover_tooltips(),
                renderers=[self.glyph_renderer],
                formatters=self.hover_names.formatters if self.hover_names is not None else {},
            )
            if self.hover_names is not None:
                self.hover_names.bind(self.hover)
            self.figure.add_tools(self.hover)
        else:
            self.hover.update(renderers=[self.glyph_renderer], tooltips=self.get_hover_tooltips())
        
    def sync_hover_names(self):
        """Row positions changed (rows dropped or reloaded): cached names are stale."""
        if self.hover_names is not None:
            self.hover_names.sync_version()

    def update_hover_tooltips(self):
        if self.hover:
            self.hover.tooltips = self.get_hover_tooltips()

    def set_visibility(self, visible:bool, apply_checkbox=False):
        self.glyph_renderer.visible = visible
        if self.lazy_loading and hasattr(self.manager, 'set_active'):
            # Hidden: drop the rows and stop loading; shown: load the viewport
            was_active = self.manager.active
            self.manager.set_active(visible)
            if not visible:
                self._view = self.manager.view_filter.copy()
                self.change_glyph_view()
            elif not was_active and getattr(self, 'loader', None) is not None:
                self.loader.force_reload()
        self.roi_manager.set_visible(visible)
        for widget in self._all_widgets:
            if widget is not self.checkbox or apply_checkbox:
                widget.disabled = not visible
        self._sync_color_picker()

    def get_visible_checkbox(self):
        self.checkbox = CheckboxGroup(labels=[f"Show {self.alias}s"], active=[0])

        def on_checkbox_change(attr, old, new):
            visible = 0 in new
            self.set_visibility(visible)

        self.checkbox.on_change('active', on_checkbox_change)

        return self.checkbox
    
    def change_glyph_view(self):
        """Show the rows of `_view`, sent as a binary mask and updated in place;
        with every row visible no mask is sent at all."""
        if self.glyph_renderer is None or self._view is None:
            return
        booleans = np.asarray(self._view, dtype=bool)
        view = self.glyph_renderer.view
        if booleans.all():
            if not isinstance(view.filter, AllIndices):
                view.filter = AllIndices()
        elif isinstance(view.filter, BooleanFilter):
            current = view.filter.booleans
            if not (isinstance(current, np.ndarray) and current.shape == booleans.shape
                    and np.array_equal(current, booleans)):
                view.filter.booleans = booleans
        else:
            view.filter = BooleanFilter(booleans=booleans)
    
    @property
    def show_not_sequenced(self):
        checkbox = getattr(self, 'notsequenced_checkbox', None)
        return checkbox is None or 0 in checkbox.active

    def get_not_sequenced_visible_checkbox(self):
        # Keep the user's choice when the widgets are rebuilt
        active = [0] if self.show_not_sequenced else []
        self.notsequenced_checkbox = CheckboxGroup(labels=[f"Show not sequenced {self.alias}s"], active=active)

        static_is_nan = self.manager.get_nan_vector()

        def on_checkbox_change(attr, old, new):
            visible = 0 in new

            if self.lazy_loading:
                is_nan = self.manager.get_nan_vector()
                # view_filter is what every reload rebuilds the view from
                current_filter = np.asarray(self.manager.view_filter, dtype=bool)
            else:
                is_nan = static_is_nan
                current_filter = np.asarray(
                    self.view,
                    dtype=bool,
                )
            
            if is_nan.shape[0] != current_filter.shape[0]:
                raise ValueError(
                    f"Unexpected mismatch between glyph view filter and sample data shape. Glyph view filter has {current_filter.shape[0]} rows, and sample data has {is_nan.shape[0]} rows"
                )

            booleans = current_filter.copy()
            booleans[is_nan] = visible

            if self.lazy_loading:
                self.manager.view_filter = booleans.copy()
            self._view = booleans
            self.change_glyph_view()

        self.notsequenced_checkbox.on_change("active", on_checkbox_change)
        return self.notsequenced_checkbox
    
    def get_download_button(self, host=settings.SERVER_HOST, port=settings.DOWNLOAD_SERVER_PORT):
        if self.lazy_loading:
            return None

        button = Button(label="Download")
        # Tags carry the completed export URL. Trigger the browser download only
        # after the server has finished writing, rather than on the initial click.
        button.js_on_change('tags', CustomJS(code="""
            const file = cb_obj.tags[0]
            if (file == null || file.url == null) return
            const link = document.createElement('a')
            link.href = file.url
            link.download = file.filename
            link.style.display = 'none'
            document.body.appendChild(link)
            link.click()
            link.remove()
        """))

        def download():
            button.disabled = True
            button.label = "Preparing download…"
            try:
                self.save_adata()
                # Each export removes only its own file, after the browser has had
                # time to fetch it; a fixed 3 s window truncated large downloads.
                timer = threading.Timer(settings.DOWNLOAD_RETENTION_SECONDS,
                                        self._remove_export, args=(self.adata_path,))
                timer.daemon = True
                timer.start()
                button.tags = [dict(url=self.static_adata_url, filename=f"{self.key}.h5ad")]
                button.label = "Download"
            except Exception:
                self._cleanup_download_file()
                button.label = "Download failed — retry"
                raise
            finally:
                button.disabled = False

        button.on_click(download)
        self.download_button = button
        return button

    def get_colnames(self):
        return self.sample.get_colnames()
    
    def get_autocomplete_input(self, doc):
        self.autocomplete_input = AutocompleteInput(
                                    title="Gene/Metadata", 
                                    value=self.get_colnames()[0], # TODO
                                    completions=self.get_colnames(), 
                                    min_characters=settings.AUTOCOMPLETE_MIN_CHARACTERS,
                                    case_sensitive=settings.AUTOCOMPLETE_CASE_SENSITIVE,
                                    search_strategy=settings.AUTOCOMPLETE_SEARCH_STRATEGY,
                                    restrict=True
                                    )
        
        def on_input_change(_, old, new):
            doc.add_next_tick_callback(lambda: self.show_colname(new))

        self.autocomplete_input.on_change('value', on_input_change)

        return self.autocomplete_input

    def get_expression_scale_select(self, doc):
        """Display genes as counts, normalized or log-normalized counts."""
        self.expression_scale_select = Select(
            title="Expression values",
            value=self.sample.expression_scale,
            options=list(self.sample.EXPRESSION_SCALES.items()),
        )

        def on_change(_, old, new):
            doc.add_next_tick_callback(lambda: self.set_expression_scale(new))

        self.expression_scale_select.on_change('value', on_change)
        return self.expression_scale_select

    def set_expression_scale(self, scale):
        if scale == self.sample.expression_scale:
            return
        with self.sample.fetch_lock:
            self.sample.expression_scale = scale
        self.manager.refresh_values()
        self.update_hover_tooltips()

    def show_colname(self, colname):
        """Display `colname`, deciding from the current parameters rather than the
        previous input value: clustering/groups change the colour behind the
        autocomplete's back, so `old` may not be what is on screen."""
        not_sequenced_color = NOT_SEQUENCED_VAL[f"{self.alias.lower()}s"]["color"]
        color_param = self.manager.get_param('color')
        alpha_param = self.manager.get_param('alpha')
        has_dynamic_alpha = alpha_param is not None and not alpha_param.is_static

        if self.sample.is_continous(colname):
            # Continuous: static colour, alpha scaled by the values of `colname`
            if color_param is not None and not color_param.is_static:
                self.add_parameter(
                    attr='color',
                    val=self.colorpicker.color if self.colorpicker else '#6baed6',
                    static=True,
                    not_sequenced_val=not_sequenced_color,
                )
            if has_dynamic_alpha:
                if alpha_param.colname != colname:
                    # Only the alpha follows the selector: extra parameters (Parameters
                    # tab) may read the same gene and must keep it
                    self.manager.update_colname(alpha_param.colname, colname, attrs=("alpha",))
                maxv = self.slider.value if self.slider else alpha_param.maxv
                self.manager.update_param('alpha', new_minv=0, new_maxv=maxv, avoid_scaling=False)
        else:
            # Categorical: one colour per category, uniform alpha
            if color_param is None or color_param.is_static or color_param.colname != colname:
                self.add_parameter(
                    attr='color',
                    colname=colname,
                    static=False,
                    not_sequenced_val=not_sequenced_color,
                    colname_hovered=True,
                )
            if has_dynamic_alpha:
                self.manager.update_param('alpha', avoid_scaling=True)

        self.update_hover_tooltips()

        self._sync_color_picker()

    def get_available_attrs(self):
        avail_attrs = []
        for param in self.AVAILABLE_ATTRS:
            if param not in [p.attr for p in self.manager.parameters]:
                avail_attrs.append(param)

        return avail_attrs
    
    
    def get_widgets(self, doc):
        self._all_widgets = []  # flat list for set_visibility

        try:
            general_widgets = [self.get_visible_checkbox()]
            self._all_widgets.append(self.checkbox)

            clear_selection = Button(label="Clear selection", width=130)
            clear_selection.on_click(lambda: self.set_selected_obs_names([]))
            general_widgets.append(clear_selection)
            self._all_widgets.append(clear_selection)

            general_widgets.append(self.get_not_sequenced_visible_checkbox())
            self._all_widgets.append(self.notsequenced_checkbox)

            if not self.lazy_loading:
                dl = self.get_download_button()
                if dl:
                    general_widgets.append(dl)
                    self._all_widgets.append(dl)

            for param in self.manager.parameters:
                if 'alpha' in param.attr and not param.is_static:
                    self.slider = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title=f"{self.alias}s alpha")

                    def update_alpha(attr, old, new, _param=param):
                        def set_new_max_alpha():
                            self.manager.update_param(_param.attr, new_maxv=new)
                        doc.add_next_tick_callback(set_new_max_alpha)

                    self.slider.on_change('value', update_alpha)
                    general_widgets.append(self.slider)
                    self._all_widgets.append(self.slider)

            # --- Color Tab ---
            color_widgets = []

            if any([not param.is_static for param in self.manager.parameters]):
                autocomplete = self.get_autocomplete_input(doc)
                color_widgets.append(autocomplete)
                self._all_widgets.append(autocomplete)
                scale_select = self.get_expression_scale_select(doc)
                color_widgets.append(scale_select)
                self._all_widgets.append(scale_select)

            for param in self.manager.parameters:
                if param.attr == 'color':
                    color = (param.value if param.is_static else self.colorpicker.color
                             if self.colorpicker is not None else settings.INITIAL_ATTRS["color"])
                    self.colorpicker = ColorPicker(title=f"Select {self.alias}s Color", color=color)

                    def update_color(attr, old, new):
                        def set_new_color():
                            current = self.manager.get_param('color')
                            if not self.colorpicker.disabled and current is not None and current.is_static:
                                self.manager.update_param('color', new_val=new)

                        doc.add_next_tick_callback(set_new_color)

                    self.colorpicker.on_change('color', update_color)
                    color_widgets.append(self.colorpicker)
                    self._all_widgets.append(self.colorpicker)

            self._sync_color_picker()

            # --- Clustering Tab ---
            def on_run_clustering(doc):
                
                self.refresh_completions()
                    
                if self.autocomplete_input is not None and self.clustering_manager.cluster_col in self.autocomplete_input.completions:
                    self.autocomplete_input.value = self.clustering_manager.cluster_col

                self.update_hover_tooltips()
                self.refresh_hover()

                self.de_manager.refresh_de_clusters()
                self.metadata_comparison_manager.refresh_clusters()
                for callback in self._clustering_listeners:
                    callback(self)
                on_markers_ready(doc)       # placeholder until the new markers are ranked

            def on_markers_ready(doc):
                if hasattr(self, "_metadata_panel"):
                    def _refresh(doc=doc):
                        self._metadata_panel.update_marker_table(self)
                    doc.add_next_tick_callback(_refresh)

            clustering_widgets = self.clustering_manager.get_widgets(doc, on_run=on_run_clustering,
                                                                     on_markers=on_markers_ready)
            self._all_widgets.extend(clustering_widgets)

            deconvolution_widgets = self.deconvolution_manager.get_widgets(doc)
            self._all_widgets.extend(deconvolution_widgets)

            groups_widgets = self.groups_manager.get_widgets(doc)
            self._all_widgets.extend(groups_widgets)

            roi_widgets = self.roi_manager.get_widgets(doc)
            self._all_widgets.extend(roi_widgets)

            de_widgets = self.de_manager.get_widgets(doc)
            self._all_widgets.extend(de_widgets)

            metadata_widgets = self.metadata_comparison_manager.get_widgets(doc)
            self._all_widgets.extend(metadata_widgets)

            # --- Tab panels ---
            # Each tab opens with its title and a "?" explaining it
            def tab(title, topic, widgets):
                return TabPanel(child=column(help_header(title, topic), *widgets), title=title)

            tab_panels = [
                tab("General", "general", general_widgets),
                tab("Color", "color", color_widgets or [Div(text="No color options.")]),
                tab("Clustering", "clustering", clustering_widgets),
                tab("Deconvolution", "deconvolution", deconvolution_widgets),
                tab("Groups", "groups", groups_widgets),
                tab("ROIs", "rois", roi_widgets),
                tab("DE", "de", de_widgets),
                tab("Metadata comparison", "metadata_comparison", metadata_widgets),
            ]

            self.add_param_form = ParamForm(self)
            form_layout = self.add_param_form.get_layout(doc)
            self._all_widgets.append(form_layout)

            tab_panels.append(TabPanel(child=form_layout, title="Parameters"))

            self.tabs_widget = Tabs(tabs=tab_panels)
            return self.tabs_widget
            
        except Exception as e:
            raise

    def get_widgets_layout(self, doc):
        return column(self.get_widgets(doc), name=self.widgets_name)
    
    def refresh_completions(self):
        # Display gene/metadata
        if self.autocomplete_input is not None:
            self.autocomplete_input.completions = self.get_colnames()
        for callback in self._columns_listeners:
            callback(self)

        # Group definition
        if self.groups_manager is not None:
            self.groups_manager.refresh_completions()

        # Metadata comparison
        if self.metadata_comparison_manager is not None:
            self.metadata_comparison_manager.refresh_columns()

        # Extra paramater
        if self.add_param_form is not None:
            self.add_param_form.refresh_completions()


