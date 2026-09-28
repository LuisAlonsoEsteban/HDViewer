
import os
from pathlib import Path

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
            ):
        
        self.figure = figure

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
        return f"{self.alias.lower()}s_widgets"
    
    @property
    def param_table_name(self):
        return f"{self.alias.lower()}s_params"
    
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

    def get_hover_tooltips(self):
        items = []
        for param in self.manager.parameters:
            if not param.is_static and param.colname_hovered:
                items.append(f"{param.colname}: @{param.colname}")

        if len(items) > 0:
            return ", ".join(items)
        else:
            return ""
        
    def update_hover_tooltips(self):
        if self.hover:
            self.hover.tooltips = self.get_hover_tooltips()

    def add_hover(self):
        if self.glyph_renderer:
            self.hover = HoverTool(tooltips=[(self.alias, self.get_hover_tooltips())], renderers=[self.glyph_renderer])
            self.figure.add_tools(self.hover)

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
    
    def get_not_sequenced_visible_checkbox(self):
        self.notsequenced_checkbox = CheckboxGroup(labels=[f"Show not sequenced {self.alias}s"], active=[0])

        static_is_nan = self.manager.get_nan_vector()

        def on_checkbox_change(attr, old, new):
            visible = 0 in new

            if self.lazy_loading:
                is_nan = self.manager.get_nan_vector()
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
            def update_colname():
                if self.sample.is_continous(old) == self.sample.is_continous(new):
                    self.manager.update_colname(old, new)
                else:
                    if self.sample.is_continous(new):
                        if self.slider:
                            self.manager.update_param('alpha', new_minv=0, new_maxv=self.slider.value)
                        if self.colorpicker:
                            self.add_parameter(
                                attr='color',
                                val=self.colorpicker.color,
                                static=True,
                                not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.alias.lower()}s"]["color"],
                            )
                    else:
                        self.manager.update_param('alpha', avoid_scaling=True)
                        if self.autocomplete_input:
                            self.add_parameter(
                                attr='color',
                                colname=self.autocomplete_input.value,
                                static=False,
                                not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.alias.lower()}s"]["color"],
                                colname_hovered=True,
                            )

                self.update_hover_tooltips()
                # self.refresh_params_layout(doc)
                
            doc.add_next_tick_callback(update_colname)

        self.autocomplete_input.on_change('value', on_input_change)

        return self.autocomplete_input
    
    def get_available_attrs(self):
        avail_attrs = []
        for param in self.AVAILABLE_ATTRS:
            if param not in [p.attr for p in self.manager.parameters]:
                avail_attrs.append(param)

        return avail_attrs
    
    
    def get_widgets(self, doc):
        self._all_widgets = []  # flat list for set_visibility

        # --- General Tab ---
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
                self.autocomplete_input.value = self.clustering_manager.CLUSTER_COL

                self.add_param_form.colname_input.completions = self.get_colnames()

            if hasattr(self, "_metadata_panel"):
                def _refresh(doc=doc):
                    self._metadata_panel.update_marker_table(self)
                doc.add_next_tick_callback(_refresh)

        clustering_widgets = self.clustering_manager.get_widgets(doc, on_run_clustering)
        self._all_widgets.extend(clustering_widgets)

        # --- Tab panels ---
        tab_panels = [
            TabPanel(child=column(*general_widgets), title="General"),
            TabPanel(child=column(*color_widgets) if color_widgets else column(Div(text="No color options.")), title="Color"),
            TabPanel(child=column(*clustering_widgets), title="Clustering"),
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

    def __init__(self, angle, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.width = width
        # self.height = height
        self.angle=angle

    def add_glyph(self):
        # width = 'width' if 'width' in self.manager.source.data else self.width
        # height = 'height' if 'height' in self.manager.source.data else self.height
        self.glyph_renderer = self.figure.rect(
                                x='x', 
                                y='y', 
                                # width=width, 
                                # height=height,
                                width='width', 
                                height='height',
                                angle=self.angle,
                                hover_color=self.hover_color,
                                **{param.attr: param.attr for param in self.manager.parameters},
                                source=self.manager.source,
                                view=CDSView(filter=BooleanFilter(booleans=self.view)),
                            )
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
