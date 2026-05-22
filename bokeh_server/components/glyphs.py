
import os

from bokeh.models import (
    HoverTool, 
    CheckboxGroup,
    AutocompleteInput,
    Slider,
    ColorPicker,
    Button,
    CustomJS,
    Div,
    ColumnDataSource, 
    DataTable, 
    TableColumn
)
from bokeh.layouts import column

from ..manager import SourceManager, DynamicSourceManager
from ..parameter import StaticParameter, Scaler, CategoricalColor

from ..loader import LazyLoader

from .utils.param_form import ParamForm

from pathlib import Path

STATIC_PATH = Path("/bokeh_server/static/data")

BASE_DIR = Path(__file__).resolve().parent  
DOWNLOADABLE_DATA_PATH = (BASE_DIR.parent / "static" / "data").resolve()

class Component:

    AVAILABLE_ATTRS = ['color', 'alpha', 'fill_alpha', 'line_width', 'line_color', 'line_alpha']
    TABLE_ATTRS = ['attr', 'is_static', 'colname', 'value', 'minv', 'maxv', 'not_sequenced_val']

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

        self.sample = sample
        self.lazy_loading = lazy_loading

        self.set_manager()

        self.hover_color = hover_color
        self.alias = alias

        self.glyph_renderer = None
        self.hover = None
        
        if not self.lazy_loading:
            self.save_adata()

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

    def set_manager(self):
        if self.lazy_loading:
            self.manager = DynamicSourceManager(self.sample)
        else:
            self.manager = SourceManager(self.sample)

    def attach_loader(self, doc):
        if self.lazy_loading and self.glyph_renderer is not None:
            loader = LazyLoader(self.manager, self.figure)
            loader.add_range_callback(self, doc)

    def refresh_widgets_layout(self, doc):
        widgets_layout = doc.get_model_by_name(self.widgets_name)
        widgets_layout.children = self.get_widgets(doc)

    def refresh_params_layout(self, doc):
        params_layout = doc.get_model_by_name(self.param_table_name)
        params_layout.children = [self.get_params_table()]

    def refresh_layout(self, doc):
        self.refresh_widgets_layout(doc)
        self.refresh_params_layout(doc)

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

    def add_glyph(self):
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
        if self.hover_color:
            self.hover = HoverTool(tooltips=[(self.alias, self.get_hover_tooltips())], renderers=[self.glyph_renderer])
            self.figure.add_tools(self.hover)

    #def set_visibility(self, visible:bool, apply_checkbox=False):
    #    self.glyph_renderer.visible = visible
    #    for widget in self.widgets:
    #        if (widget is not self.checkbox and widget is not self.add_param_button) or apply_checkbox:
    #            widget.disabled = not visible

    def set_visibility(self, visible:bool, apply_checkbox=False):
        self.glyph_renderer.visible = visible
        for widget in self.widgets:
            is_clustering_widget = (
                hasattr(self, 'resolution_slider') and widget is self.resolution_slider
            ) or (
                hasattr(widget, 'label') and widget.label in ["Run clustering", "Running..."]
            )
            if is_clustering_widget:
                continue
            if (widget is not self.checkbox and widget is not self.add_param_button) or apply_checkbox:
                widget.disabled = not visible

    def get_visible_checkbox(self):
        self.checkbox = CheckboxGroup(labels=[f"Show {self.alias}s"], active=[0])

        def on_checkbox_change(attr, old, new):
            visible = 0 in new
            self.set_visibility(visible)

        self.checkbox.on_change('active', on_checkbox_change)

        return self.checkbox
    
    def get_download_button(self, host='0.0.0.0', port=5002):
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
        autocomplete_input = AutocompleteInput(
                                    title="Gene/Metadata", 
                                    value=self.get_colnames()[0], # TODO
                                    completions=self.get_colnames(), 
                                    min_characters=2, 
                                    case_sensitive=False, 
                                    restrict=True
                                    )
        
        def on_input_change(_, old, new):
            def update_colname():
                self.manager.update_colname(old, new)
                self.update_hover_tooltips()
                self.refresh_params_layout(doc)

            doc.add_next_tick_callback(update_colname)

        autocomplete_input.on_change('value', on_input_change)

        return autocomplete_input
    
    def get_available_attrs(self):
        avail_attrs = []
        for param in self.AVAILABLE_ATTRS:
            if param not in [p.attr for p in self.manager.parameters]:
                avail_attrs.append(param)

        return avail_attrs
    
    def get_add_param_button(self, doc):
        self.add_param_button = Button(label='Add parameter')

        def handler():
            if not hasattr(self, 'add_param_form'):
                self.add_param_form = ParamForm(self)
            form_layout = self.add_param_form.get_layout(doc)
            parent_layout = doc.get_model_by_name('spots_widgets')
            
            def append_form():
                self.add_param_button.visible = False
                parent_layout.children.append(form_layout)
                
            doc.add_next_tick_callback(append_form)
    
        self.add_param_button.on_click(handler)
        return self.add_param_button
    
    def get_widgets(self, doc):
        self.widgets = [self.get_visible_checkbox()]

        if not self.lazy_loading:
            self.widgets.append(self.get_download_button())

        if any([not param.is_static for param in self.manager.parameters]):
            self.widgets.append(self.get_autocomplete_input(doc))

        for param in self.manager.parameters:

            if 'alpha' in param.attr and not param.is_static:
                widget = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title=f"{self.alias}s alpha")

                def update_alpha(attr, old, new, _param=param):
                    def set_new_max_alpha():
                        self.manager.update_param(_param.attr, new_maxv=new)
                    doc.add_next_tick_callback(set_new_max_alpha)
                
                widget.on_change('value', update_alpha)

            elif 'color' in param.attr:
                widget = ColorPicker(title=f"Select {self.alias}s Color", color='#6baed6')

                def update_color(attr, old, new, _param=param):
                    def set_new_color():
                        self.manager.update_param('color', new_val=new)
                        self.refresh_params_layout(doc)

                    doc.add_next_tick_callback(set_new_color)

                widget.on_change('color', update_color)

            self.widgets.append(widget)

        if len(self.get_available_attrs()) > 0:
            self.widgets.append(self.get_add_param_button(doc))

        # Clustering
        self.resolution_slider = Slider(
            start=0.1, end=2.0, value=0.9, step=0.1, title="Clustering resolution"
        )

        clustering_button = Button(label="Run clustering", button_type="success")
        clustering_button = Button(label="Run clustering", button_type="success")
        clustering_info = Div(text="<p style='color:gray;font-size:12px'>Pulsa Run clustering para ver resultados</p>")

        def run_clustering():
            resolution = self.resolution_slider.value
            def do_clustering():
                clustering_button.label = "Running..."
                clustering_button.disabled = True
                try:
                    self.sample.run_clustering(resolution=resolution)
                    autocomplete = next(
                        (w for w in self.widgets if isinstance(w, AutocompleteInput)), None
                    )
                    if autocomplete:
                        autocomplete.completions = self.get_colnames()
                        autocomplete.value = 'leiden'
                    counts = self.sample.adata.obs['leiden'].value_counts().sort_index()
                    info_html = "<b>Clusters encontrados:</b><br>"
                    for cluster, count in counts.items():
                        info_html += f"Cluster {cluster}: {count} spots<br>"
                    # Añadir genes marcadores al panel
                    markers = self.sample.get_marker_genes(n_genes=5)
                    if markers:
                        info_html += "<br><b>Top genes por cluster:</b><br>"
                        for cluster, genes in markers.items():
                            gene_names = ", ".join([g[0] for g in genes])
                            info_html += f"Cluster {cluster}: {gene_names}<br>"
                    clustering_info.text = info_html
                    clustering_button.label = "Run clustering"
                except Exception as e:
                    clustering_button.label = f"Error: {str(e)[:30]}"
                finally:
                    clustering_button.disabled = False
            doc.add_next_tick_callback(do_clustering)

        clustering_button.on_click(run_clustering)

        self.widgets.append(self.resolution_slider)
        self.widgets.append(clustering_button)
        self.widgets.append(clustering_info)
        return self.widgets
    
    def get_widgets_layout(self, doc):
        return column(*self.get_widgets(doc), name=self.widgets_name)
    

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
                                )
        
class Squares(Component):

    def __init__(self, width, height, angle, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.angle=angle

    def add_glyph(self):
        self.glyph_renderer = self.figure.rect(
                                x='x', 
                                y='y', 
                                width=self.width, 
                                height=self.height,
                                angle=self.angle,
                                hover_color=self.hover_color,
                                **{param.attr: param.attr for param in self.manager.parameters},
                                source=self.manager.source,
                            )
        
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
                                )



