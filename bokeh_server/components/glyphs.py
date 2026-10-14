
import os
from pathlib import Path
from .utils.groups import GroupsManager
from .bin_info import BinInfoOverlay
import numpy as np

from bokeh.models import (
    Div,
    HoverTool, 
    CheckboxGroup,
    AutocompleteInput,
    Slider,
    ColorPicker,
    Button,
    CustomJS,
    Tabs,
    TabPanel,
    ColumnDataSource, 
    DataTable, 
    TableColumn,

    BooleanFilter,
    CDSView
)
from bokeh.layouts import column

from ..manager import SourceManager, DynamicSourceManager
from ..parameter import StaticParameter, Scaler, CategoricalColor

from ..loader import LazyLoader

from .utils.param_form import ParamForm
from .utils.clustering import ClusteringManager

from . import NOT_SEQUENCED_VAL


STATIC_PATH = Path("/bokeh_server/static/data")

BASE_DIR = Path(__file__).resolve().parent  
DOWNLOADABLE_DATA_PATH = (BASE_DIR.parent / "static" / "data").resolve()
# static/data está en .gitignore: en un clon nuevo no existe y save_adata fallaría.
DOWNLOADABLE_DATA_PATH.mkdir(parents=True, exist_ok=True)

class Component:

    AVAILABLE_ATTRS = ['color', 'alpha', 'fill_alpha', 'line_width', 'line_color', 'line_alpha']
    TABLE_ATTRS = ['attr', 'is_static', 'colname', 'value', 'minv', 'maxv', 'not_sequenced_val']

    CLUSTER_COLOR_LABEL = "cluster_color"

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

        self.hover_color = hover_color
        self.alias = alias

        self.glyph_renderer = None
        self.hover = None

        self.checkbox = None
        self.slider = None
        self.colorpicker = None
        self.autocomplete_input = None

        self.clustering_manager = ClusteringManager(self)
        self.groups_manager = GroupsManager(self)
        self.selected_obs_names = set()
        self._selection_listeners = []
        self._suspend_selection_sync = False
        self._selection_callback_bound = False
        self.sample_type = f"{self.alias.lower()}s"

        self._view = None

    @property
    def sample(self):
        return self.manager.sample

    @property
    def static_adata_path(self):
        if self.adata_path:
            return STATIC_PATH / self.adata_path.name

    @property
    def widgets_name(self):
        return f"{self.key}_widgets"
    
    @property
    def param_table_name(self):
        return f"{self.key}_params"
    
    @property
    def view(self):
        if not self._view:
            self._view = [True] * self.manager.n_geometries
        return self._view
    
    def save_adata(self):
        if not self.lazy_loading:
            i = 0
            while(True):
                adata_path = (DOWNLOADABLE_DATA_PATH / f"{i}.h5ad").resolve()
                if not os.path.exists(adata_path):
                    break
                i += 1
            self.adata_path = adata_path
            self.sample.adata.write_h5ad(self.adata_path)

    def attach_loader(self, doc):
        if self.lazy_loading and self.glyph_renderer is not None:
            loader = LazyLoader(self.manager, self.figure)
            loader.add_range_callback(self, doc)
            self.loader = loader

            def _initial_load():
                loader.force_reload()
            doc.add_next_tick_callback(_initial_load)

    def refresh_widgets_layout(self, doc):
        widgets_layout = doc.get_model_by_name(self.widgets_name)
        widgets_layout.children = [self.get_widgets(doc)]

    def refresh_params_layout(self, doc):
        params_layout = doc.get_model_by_name(self.param_table_name)
        params_layout.children = [self.get_params_table()]

    def refresh_layout(self, doc):
        self.refresh_widgets_layout(doc)
        # self.refresh_params_layout(doc)

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
        self.remove_existing_hover()
        self.add_hover()

    def refresh_figure(self):
        self.refresh_glyph()
        self.refresh_hover()

    @property
    def obs_name_values(self):
        return np.asarray(self.manager.source.data.get("obs_name", []), dtype=str)

    def get_selected_obs_names(self):
        return set(self.selected_obs_names)

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
            indices = np.flatnonzero(np.isin(obs_names, list(self.selected_obs_names))).tolist()
            
        self._suspend_selection_sync = True
        self.glyph_renderer.data_source.selected.indices = indices
        self._suspend_selection_sync = False

    def set_selected_obs_names(self, obs_names, notify=True):
        self.selected_obs_names = {str(obs_name) for obs_name in obs_names}
        self.apply_selected_obs_names()
        if notify:
            self.notify_selection_changed()

    def _on_source_selection_change(self, attr, old, new):
        if self._suspend_selection_sync:
            return
        obs_names = self.obs_name_values
        selected_indices = np.asarray(new, dtype=int)
        valid_indices = selected_indices[selected_indices < obs_names.size]
        self.selected_obs_names = set(obs_names[valid_indices].tolist())
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
        pass

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
                self.manager.add_parameter(
                    CategoricalColor(
                        colname,
                        attr=attr,
                        static_not_sequenced=static_not_sequenced,
                        not_sequenced_val=not_sequenced_val,
                        colname_hovered=colname_hovered
                        )
                )
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

    #def get_hover_tooltips(self):
    #    items = []
    #    for param in self.manager.parameters:
    #        if not param.is_static and param.colname_hovered:
    #            items.append(f"{param.colname}: @{param.colname}")
#
    #    if len(items) > 0:
    #        return ", ".join(items)
    #    else:
    #        return ""

    def get_hover_tooltips(self):
        items = []
        for param in self.manager.parameters:
            if (not param.is_static and param.colname_hovered
                    and param.colname not in ('custom_clusters', 'leiden')
                    and param.attr != 'alpha'):
                items.append(f"{param.colname}: @{param.colname}")

        # Añadir leiden solo si tiene valores reales (no todos None)
        if 'custom_clusters' in self.manager.source.data:
            items.append("Group: @custom_clusters")
        if 'leiden' in self.manager.source.data:
            items.append("Cluster: @leiden")

        if len(items) > 0:
            return ", ".join(items)
        else:
            return "@obs_name"

    def add_hover(self):
        if self.glyph_renderer:
            self.hover = HoverTool(
                tooltips=self.get_hover_tooltips(),
                renderers=[self.glyph_renderer]
            )
            self.figure.add_tools(self.hover)
        
    def update_hover_tooltips(self):
        if self.hover:
            self.hover.tooltips = self.get_hover_tooltips()

    def set_visibility(self, visible:bool, apply_checkbox=False):
       self.glyph_renderer.visible = visible
       for widget in self._all_widgets:
           if widget is not self.checkbox or apply_checkbox:
               widget.disabled = not visible

    def get_visible_checkbox(self):
        self.checkbox = CheckboxGroup(labels=[f"Show {self.alias}s"], active=[0])

        def on_checkbox_change(attr, old, new):
            visible = 0 in new
            self.set_visibility(visible)

        self.checkbox.on_change('active', on_checkbox_change)

        return self.checkbox
    
    def change_glyph_view(self):
        if self._view:
            self.glyph_renderer.view.filter = BooleanFilter(booleans=self._view)
    
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
            self._view = booleans.tolist()
            self.change_glyph_view()

        self.notsequenced_checkbox.on_change("active", on_checkbox_change)
        return self.notsequenced_checkbox
    
    def get_download_button(self, host='0.0.0.0', port=5002):
        self.save_adata()
        if self.static_adata_path:
            self.download_button = Button(label="Download")
            self.download_button.js_on_click(
                CustomJS(
                        code=f"""
                        window.open('{self.static_adata_path}');
                        """,
                )
            )

            return self.download_button

    def get_colnames(self):
        return self.sample.get_colnames()
    
    def get_autocomplete_input(self, doc):
        self.autocomplete_input = AutocompleteInput(
                                    title="Gene/Metadata", 
                                    value=self.get_colnames()[0], # TODO
                                    completions=self.get_colnames(), 
                                    min_characters=2, 
                                    case_sensitive=False, 
                                    restrict=True
                                    )
        
        def on_input_change(_, old, new):
            doc.add_next_tick_callback(lambda: self.show_colname(new))

        self.autocomplete_input.on_change('value', on_input_change)

        return self.autocomplete_input

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
                    self.manager.update_colname(alpha_param.colname, colname)
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

            for param in self.manager.parameters:
                if 'color' in param.attr:
                    self.colorpicker = ColorPicker(title=f"Select {self.alias}s Color", color='#6baed6')

                    def update_color(attr, old, new, _param=param):
                        def set_new_color():
                            self.manager.update_param('color', new_val=new)
                            # self.refresh_params_layout(doc)
                        doc.add_next_tick_callback(set_new_color)

                    self.colorpicker.on_change('color', update_color)
                    color_widgets.append(self.colorpicker)
                    self._all_widgets.append(self.colorpicker)

            # --- Clustering Tab ---
            def on_run_clustering(doc):
                if self.autocomplete_input:
                    self.autocomplete_input.completions = self.get_colnames()
                    # Reflect what is displayed, so the previous gene can be picked again
                    if self.clustering_manager.cluster_col in self.autocomplete_input.completions:
                        self.autocomplete_input.value = self.clustering_manager.cluster_col

                    self.add_param_form.colname_input.completions = self.get_colnames()

                self.update_hover_tooltips()
                self.refresh_hover()
                
                if hasattr(self, "_metadata_panel"):
                    def _refresh(doc=doc):
                        self._metadata_panel.update_marker_table(self)
                    doc.add_next_tick_callback(_refresh)

            groups_widgets = self.groups_manager.get_widgets(doc)
            self._all_widgets.extend(groups_widgets)

            clustering_widgets = self.clustering_manager.get_widgets(doc, on_run=on_run_clustering)
            self._all_widgets.extend(clustering_widgets)

            # --- Tab panels ---
            tab_panels = [
                TabPanel(child=column(*general_widgets), title="General"),
                TabPanel(child=column(*color_widgets) if color_widgets else column(Div(text="No color options.")), title="Color"),
                TabPanel(child=column(*clustering_widgets), title="Clustering"),
                TabPanel(child=column(*groups_widgets), title="Groups"),
            ]

            for param in self.manager.parameters:
                if not param.is_static and 'alpha' not in param.attr and 'color' not in param.attr:
                    tab_panels.append(
                        TabPanel(child=column(Div(text=f"Parameter: <b>{param.attr}</b>")), title=param.attr)
                    )

            if len(self.get_available_attrs()) > 0:
                self.add_param_form = ParamForm(self)
                form_layout = self.add_param_form.get_layout(doc)
                self._all_widgets.append(form_layout)
                tab_panels.append(
                    TabPanel(child=column(form_layout), title="＋ Add param")
                )

            self.tabs_widget = Tabs(tabs=tab_panels)
            return self.tabs_widget
            
        except Exception as e:
            raise

    def get_widgets_layout(self, doc):
        return column(self.get_widgets(doc), name=self.widgets_name)
    

    @property
    def params_source(self):
        source_dict = {
            column: [getattr(param, column, None) for param in self.manager.parameters]
            for column in self.TABLE_ATTRS
        }

        if hasattr(self, '_params_source'):
            self._params_source.data = source_dict
        else:
            self._params_source = ColumnDataSource(data=source_dict)

        return self._params_source

    def get_params_table(self):
        columns = [TableColumn(field=column, title=column) for column in self.TABLE_ATTRS]
        data_table = DataTable(source=self.params_source, columns=columns, width=400, height=300)
        return data_table
    
    def get_table_layout(self):
        return column(self.get_params_table(), name=self.param_table_name)

    

class Circles(Component):

    def __init__(self, radius, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius

    def add_glyph(self):
        self.glyph_renderer = self.figure.circle(
                                    x='x', 
                                    y='y', 
                                    radius=self.radius, 
                                    hover_color=self.hover_color,
                                    **{param.attr: param.attr for param in self.manager.parameters},
                                    source=self.manager.source,
                                    view=CDSView(filter=BooleanFilter(booleans=self.view)),
                                )
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()
        
class Squares(Component):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bin_info = BinInfoOverlay(self)

    def add_glyph(self):
        # Both native footprints and axis-aligned aggregate cells supply their
        # corners; x/y remain scalar centers for lazy loading and selection.
        self.glyph_renderer = self.figure.patches(
            xs='xs',
            ys='ys',
            hover_color=self.hover_color,
            **{param.attr: param.attr for param in self.manager.parameters},
            source=self.manager.source,
            view=CDSView(filter=BooleanFilter(booleans=self.view)),
        )
        self.bin_info.bind_renderer(self.glyph_renderer)
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()
        
class Patches(Component):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_glyph(self):
        self.glyph_renderer = self.figure.patches(
                                    xs='x', 
                                    ys='y',
                                    hover_color=self.hover_color,
                                    **{param.attr: param.attr for param in self.manager.parameters},
                                    source=self.manager.source,
                                    view=CDSView(filter=BooleanFilter(booleans=self.view)),
                                )
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()
