from bokeh.layouts import column, row
from bokeh.models import Button, ColumnDataSource, Div, Select, TabPanel, Tabs
from bokeh.plotting import figure

from . import Panel


class EmbeddingPanel(Panel):

    METHOD_LABELS = {
        "pca": "PCA",
        "umap": "UMAP",
        "tsne": "t-SNE",
    }

    def __init__(self):
        super().__init__(name="embedding_panel")

        self.components = {}
        self._syncing = False
        self._component_views = {}
        self.tabs = Tabs(tabs=[])

        self.set_children([self.tabs])

    def set_components(self, components):
        component_map = {}
        for component in components:
            component_map[component.sample_type] = component
            component.add_selection_listener(self._on_component_selection_change)
        self.components = component_map

        existing_keys = set(self._component_views.keys())
        current_keys = set(self.components.keys())

        for sample_type, component in self.components.items():
            if sample_type not in self._component_views:
                self._component_views[sample_type] = self._create_component_view(component)

        for sample_type in existing_keys - current_keys:
            del self._component_views[sample_type]

        self.tabs.tabs = [
            self._component_views[sample_type]["panel"]
            for sample_type in self.components.keys()
        ]

        if len(self.tabs.tabs) == 0:
            self.tabs.active = 0
        elif self.tabs.active >= len(self.tabs.tabs):
            self.tabs.active = 0

    def _create_component_view(self, component):
        source = ColumnDataSource(data={
            "obs_name": [],
            "x": [],
            "y": [],
            "color": [],
        })
        source.selected.on_change(
            "indices",
            lambda attr, old, new, sample_type=component.sample_type: self._on_embedding_selection_change(sample_type, new),
        )

        method_select = Select(
            title="Method",
            options=[(key, label) for key, label in self.METHOD_LABELS.items()],
            value="pca",
        )
        run_button = Button(label="Run embedding", button_type="success")
        status_div = Div(
            text=f"Select a method and run an embedding for {component.alias.lower()}s."
        )

        plot = figure(
            title=f"{component.alias} embedding",
            height=350,
            visible=False,
            tools="pan,wheel_zoom,box_select,lasso_select,tap,save,reset",
        )
        plot.scatter(
            x="x",
            y="y",
            source=source,
            size=6,
            color="color",
            alpha=0.8,
        )
        plot.toolbar.logo = None

        def _on_run():
            self._run_embedding(component.sample_type)

        run_button.on_click(_on_run)

        content = column(
            row(method_select, run_button, sizing_mode="stretch_width"),
            status_div,
            plot,
        )
        panel = TabPanel(child=content, title=component.alias)

        return {
            "component": component,
            "source": source,
            "method_select": method_select,
            "run_button": run_button,
            "status_div": status_div,
            "plot": plot,
            "panel": panel,
        }

    def _run_embedding(self, sample_type):
        view = self._component_views.get(sample_type)
        if view is None:
            return

        component = view["component"]
        method = view["method_select"].value

        try:
            df = component.sample.get_embedding_dataframe(method)
            df["color"] = self._get_colors_for_component(component, df["obs_name"].tolist())
            view["source"].data = {
                "obs_name": df["obs_name"].astype(str).tolist(),
                "x": df["x"].tolist(),
                "y": df["y"].tolist(),
                "color": df["color"].tolist(),
            }
            view["plot"].xaxis.axis_label = f"{self.METHOD_LABELS[method]} 1"
            view["plot"].yaxis.axis_label = f"{self.METHOD_LABELS[method]} 2"
            view["plot"].visible = True
            view["status_div"].text = (
                f"<p style='color:green'>Loaded {len(df)} {component.alias.lower()} embeddings.</p>"
            )
            self._update_embedding_selection_from_component(sample_type)
        except Exception as exc:
            view["status_div"].text = f"<p style='color:red'>Error: {exc}</p>"
            view["plot"].visible = False

    def _get_colors_for_component(self, component, obs_names):
        source_data = component.manager.source.data
        obs_name_values = [str(value) for value in source_data.get("obs_name", [])]
        color_values = source_data.get("color", [])
        color_map = {
            obs_name: color
            for obs_name, color in zip(obs_name_values, color_values)
        }

        default_color = "#6baed6"
        if len(color_values) > 0:
            default_color = color_values[0]

        return [color_map.get(str(obs_name), default_color) for obs_name in obs_names]

    def _update_embedding_selection_from_component(self, sample_type):
        view = self._component_views.get(sample_type)
        if view is None:
            return

        component = view["component"]
        selected = component.get_selected_obs_names()
        obs_names = [str(value) for value in view["source"].data.get("obs_name", [])]
        indices = [i for i, obs_name in enumerate(obs_names) if obs_name in selected]

        self._syncing = True
        view["source"].selected.indices = indices
        self._syncing = False

    def _on_embedding_selection_change(self, sample_type, new):
        if self._syncing:
            return

        view = self._component_views.get(sample_type)
        if view is None:
            return

        component = view["component"]
        obs_names = [str(value) for value in view["source"].data.get("obs_name", [])]
        selected_obs_names = {obs_names[i] for i in new if i < len(obs_names)}
        component.set_selected_obs_names(selected_obs_names, notify=False)

    def _on_component_selection_change(self, component):
        self._update_embedding_selection_from_component(component.sample_type)
