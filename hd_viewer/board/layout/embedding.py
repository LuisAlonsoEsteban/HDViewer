
import pandas as pd
import numpy as np

from bokeh.layouts import column, row
from bokeh.models import AutocompleteInput, Button, ColumnDataSource, HoverTool, Select, TabPanel, Tabs
from bokeh.plotting import figure

from ... import settings
from ...core import embedding
from ...settings.config import embedding_kwargs
from . import Panel
from ..utils.loading import LoadingLogger
from ...help import help_header


class EmbeddingPanel(Panel):

    METHOD_LABELS = {
        "pca": "PCA",
        "umap": "UMAP",
        "tsne": "t-SNE",
    }

    DEFAULT_COLOR = "#6baed6"
    MISSING_COLOR = "#555555"

    def __init__(self, doc):
        super().__init__(name="embedding_panel")
        self.doc = doc

        self.components = {}
        self._syncing = False
        self._component_views = {}
        self._embedding_keys = {
            "pca": "X_pca_2d",
            "umap": "X_umap",
            "tsne": "X_tsne",
        }
        self.tabs = Tabs(tabs=[])

        self.set_children([help_header("Embeddings", "embedding"), self.tabs])

    def _get_component_context(self, component):
        sample = getattr(component.manager, "sample", None)
        clustering_manager = getattr(component, "clustering_manager", None)
        if sample is None or clustering_manager is None:
            raise RuntimeError(f"Embedding requires a component with both a sample and a clustering manager: {component!r}")
        return sample, clustering_manager

    def set_components(self, components):
        component_map = {}
        for component in components:
            component_map[component.key] = component
            component.add_selection_listener(self._on_component_selection_change)
            component.add_columns_listener(self._on_component_columns_change)
            component.add_clustering_listener(self._on_component_clustered)
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
            "value": [],
        })
        source.selected.on_change(
            "indices",
            lambda attr, old, new, sample_type=component.key: self._on_embedding_selection_change(sample_type, new),
        )

        method_select = Select(
            title="Method",
            options=[(key, label) for key, label in self.METHOD_LABELS.items()],
            value="pca",
        )
        run_button = Button(label="Run embedding", button_type="success")
        _, clustering_manager = self._get_component_context(component)
        color_input = AutocompleteInput(
            title="Color by",
            completions=component.get_colnames(),
            value=clustering_manager.cluster_col if clustering_manager.is_clustered else "",
            placeholder="Gene or metadata column",
            min_characters=settings.AUTOCOMPLETE_MIN_CHARACTERS,
            case_sensitive=settings.AUTOCOMPLETE_CASE_SENSITIVE,
            search_strategy=settings.AUTOCOMPLETE_SEARCH_STRATEGY,
            restrict=True,
        )
        color_input.on_change(
            "value", lambda attr, old, new, sample_type=component.key: self._recolor(sample_type))
        status = LoadingLogger(
            self.doc,
            text=f"Select a method and run an embedding for {component.name}.",
        )

        plot = figure(
            title=f"{component.name} embedding",
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
        plot.add_tools(HoverTool(tooltips=[("", "@obs_name"), ("value", "@value")]))
        plot.toolbar.logo = None

        def _on_run():
            self._run_embedding(component.key)

        run_button.on_click(_on_run)

        content = column(
            row(method_select, run_button, sizing_mode="stretch_width"),
            color_input,
            status.div,
            plot,
        )
        panel = TabPanel(child=content, title=component.name)

        return {
            "component": component,
            "source": source,
            "method_select": method_select,
            "run_button": run_button,
            "color_input": color_input,
            "status": status,
            "plot": plot,
            "panel": panel,
        }

    def compute_embedding(self, component, method: str, n_pcs=30, n_neighbors=15):
        sample, clustering_manager = self._get_component_context(component)
        return embedding.compute(sample, clustering_manager, method, n_pcs=n_pcs, n_neighbors=n_neighbors)

    def get_embedding_dataframe(self, component, method: str, resolution_level=None, *, n_pcs=30, n_neighbors=15) -> pd.DataFrame:
        sample, clustering_manager = self._get_component_context(component)
        return embedding.frame(sample, clustering_manager, method, n_pcs=n_pcs, n_neighbors=n_neighbors)

    def _run_embedding(self, sample_type):
        view = self._component_views.get(sample_type)
        if view is None:
            return

        component = view["component"]
        method = view["method_select"].value
        status = view["status"]
        if status.running:
            return

        # Embed at the level currently displayed (Bins may be aggregated)
        resolution_level = getattr(component.manager, "effective_resolution_level", None)

        def _work():
            return self.get_embedding_dataframe(component, method, resolution_level=resolution_level, **embedding_kwargs())

        def _on_error(exc):
            view["run_button"].disabled = False
            status.error(exc)
            view["plot"].visible = False

        view["run_button"].disabled = True
        status.run(
            f"{self.METHOD_LABELS[method]} embedding",
            _work,
            on_done=lambda df: self._on_embedding_done(sample_type, method, df),
            on_error=_on_error,
        )

    def _on_embedding_done(self, sample_type, method, df):
        view = self._component_views.get(sample_type)
        if view is None:
            return

        component = view["component"]
        view["run_button"].disabled = False
        try:
            obs_names = df["obs_name"].astype(str).tolist()
            colors, values = self.colors_for(component, view["color_input"].value, obs_names)
            view["source"].data = {
                "obs_name": obs_names,
                "x": df["x"].tolist(),
                "y": df["y"].tolist(),
                "color": colors,
                "value": values,
            }
            view["plot"].xaxis.axis_label = f"{self.METHOD_LABELS[method]} 1"
            view["plot"].yaxis.axis_label = f"{self.METHOD_LABELS[method]} 2"
            view["plot"].visible = True
            view["status"].success(f"Loaded {len(df)} {component.alias.lower()} embeddings.")
            self._update_embedding_selection_from_component(sample_type)
        except Exception as exc:
            view["status"].error(exc)
            view["plot"].visible = False

    def colors_for(self, component, colname, obs_names):
        """Colours (and hover values) of `obs_names` by `colname`, over every
        embedded observation rather than only the map's loaded viewport.
        Clusters and groups keep their assigned colours; genes follow the
        layer's expression scale."""
        palette = {}
        if colname and colname in component.get_colnames():
            palette = dict(component.known_category_colors(colname))
            color_param = component.manager.get_param("color")
            if color_param is not None and not color_param.is_static and color_param.colname == colname:
                palette = {**getattr(color_param, "category_colors", {}), **palette}
        return embedding.colours(component.manager.sample, colname, obs_names, palette)

    def _recolor(self, sample_type):
        view = self._component_views.get(sample_type)
        if view is None:
            return
        obs_names = list(view["source"].data.get("obs_name", []))
        if not obs_names:
            return
        try:
            colors, values = self.colors_for(view["component"], view["color_input"].value, obs_names)
        except Exception as exc:
            view["status"].error(exc)
            return
        view["source"].data.update(color=colors, value=values)

    def _on_component_columns_change(self, component):
        view = self._component_views.get(component.key)
        if view is not None:
            view["color_input"].completions = component.get_colnames()

    def _on_component_clustered(self, component):
        """New clusters have their colours assigned: show them on the embedding."""
        view = self._component_views.get(component.key)
        if view is None:
            return
        cluster_col = component.clustering_manager.cluster_col
        view["color_input"].completions = component.get_colnames()
        if view["color_input"].value == cluster_col:
            self._recolor(component.key)
        else:
            view["color_input"].value = cluster_col     # recolours through on_change

    def _update_embedding_selection_from_component(self, sample_type):
        view = self._component_views.get(sample_type)
        if view is None:
            return

        component = view["component"]
        obs_names = [str(value) for value in view["source"].data.get("obs_name", [])]
        indices = component.selection_indices_for(obs_names)

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
        self._update_embedding_selection_from_component(component.key)
