
import pandas as pd
import numpy as np

from bokeh.models import ColumnDataSource

from ... import settings
from ...core.layer_state import (  # noqa: F401  (re-exported)
    LayerRows, ParameterSet, _compact, parameter_column, source_columns,
)


class SourceManager(ParameterSet):
    """Static layers (spots): every row is loaded once into a ColumnDataSource."""

    def __init__(self, sample):
        self.sample = sample
        self.component = None
        self.parameters = []
        self.is_empty = True

        self.is_sequenced = []

        self.source = ColumnDataSource(data={'x': [], 'y': []})
        self.cache = pd.DataFrame()

    @property
    def n_geometries(self):
        return self.sample.n_geometries

    @property
    def _data(self):
        return self.source.data

    def _set(self, name, values):
        self.source.data[name] = values

    def reset(self):
        self.source = ColumnDataSource(data={'x': [], 'y': []})
        for param in self.parameters:
            self.source.data[param.attr] = []
        self.is_empty = True

    def get_nan_vector(self):
        return self.sample.get_obs_names().isna().to_numpy()

    @property
    def loaded_obs_names(self):
        """Observation names of the source rows, in row order."""
        return np.asarray(self.source.data.get('obs_name', []), dtype=object)

    row_version = 0
    
    def add_parameter(self, param):
        existing_param = self.get_param(param.attr)
        if existing_param is not None:
            self.parameters.remove(existing_param)
            if not existing_param.is_static and existing_param.colname in self.cache.columns:
                self.cache.drop(columns=existing_param.colname, inplace=True)
            if existing_param.attr in self.source.data.keys():
                del self.source.data[existing_param.attr]

        self.parameters.append(param)

        if param.is_static:
            param.set_is_sequenced(self.is_sequenced)
        else:
            if self.sample.has_adata():
                if self.sample.is_continous(param.colname):
                    param.min_counts, param.max_counts = self.sample.get_range(param.colname)
            else:
                raise ValueError("Non-static parameters cannot be added without an AnnData associated")
        
        if not self.is_empty:
            if param.is_static:
                self.source.data[param.attr] = param.get_values()
            else:
                if param.colname in self.cache.columns:
                    self.cache.drop(columns=param.colname, inplace=True)

                if param.colname not in self.cache.columns:
                    if param.colname in self.sample.adata.obs.columns:
                        # Leer de adata.obs y alinear con obs_names del source
                        obs_series = self.sample.adata.obs[param.colname]
                        obs_names = pd.Series(self.loaded_obs_names)

                        aligned = obs_names.map(obs_series)
                        
                        self.cache[param.colname] = aligned.values
                    else:
                        colname_data = self.sample.fetch_data(
                            colnames=param.colname, return_coords=False
                        )
                        self.cache[param.colname] = colname_data[param.colname]

                valores = param.get_values(
                    pd.Series(self.cache[param.colname]).copy()
                )
                self.source.data[param.attr] = parameter_column(valores, len(self.source.data['x']))
        else:
            self.reset()

    def remove_parameter(self, attr):
        """Drop the parameter bound to `attr`; the glyph falls back to its default."""
        param = self.get_param(attr)
        if param is None:
            return
        self.parameters.remove(param)
        data = dict(self.source.data)
        data.pop(attr, None)
        # Its column is only needed while another parameter still uses it
        if not param.is_static and param.colname not in self.get_colnames():
            data.pop(param.colname, None)
            if param.colname in self.cache.columns:
                self.cache.drop(columns=param.colname, inplace=True)
        self.source.data = data

    def update_colname(self, old_colname, new_colname, attrs=None):
        """Point the parameters reading `old_colname` at `new_colname`; with `attrs`,
        only those attributes (other parameters may read the same column)."""
        for param in self.get_params_by_colname(old_colname):
            if attrs is not None and param.attr not in attrs:
                continue
            param.colname = new_colname

            if hasattr(param, "min_counts"):
                min_val, max_val = self.sample.get_range(param.colname)
                if isinstance(min_val, str):
                    param.min_counts = 0
                    param.max_counts = 1
                else:
                    param.min_counts = min_val
                    param.max_counts = max_val

            if param.colname not in self.cache.columns:
                self.cache[param.colname] = self.sample.get_value(colname=param.colname)

            self.source.data[param.attr] = param.get_values(self.cache[param.colname].copy())
            self.source.data[param.colname] = self.cache[param.colname]

        if not self.get_params_by_colname(old_colname):
            self.source.data.pop(old_colname, None)
        
    def refresh_values(self):
        """Re-read gene values after the sample's expression scale changed."""
        for param in self._gene_parameters():
            if param.colname in self.cache.columns:
                self.cache.drop(columns=param.colname, inplace=True)
            self.update_colname(param.colname, param.colname, attrs=(param.attr,))

    def fetch_data(self):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames())
        colnames = self.get_colnames()
        if len(colnames) > 0:
            self.is_sequenced.clear()
            self.is_sequenced += (~fetched_data[self._sequenced_column(colnames)].isna()).to_numpy().tolist()

        for param in self.parameters:
            if not param.is_static:
                if param.colname not in self.cache.columns:
                    self.cache[param.colname] = fetched_data[param.colname]
                col_data = (
                    self.cache[param.colname]
                    if param.colname in self.cache.columns
                    else fetched_data[param.colname]
                )
                fetched_data[param.attr] = param.get_values(col_data.copy())
            else:
                fetched_data[param.attr] = param.get_values()

        cols_to_drop = [
            param.colname for param in self.parameters
            if not param.is_static and not param.colname_hovered
        ]
        fetched_data.drop(cols_to_drop, axis=1, inplace=True)
        self.source.data = source_columns(fetched_data)
        if self.component is not None:
            self.component.apply_selected_obs_names()
        self.is_empty = False

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self._set_column(param.attr, parameter_column(
                param.update(*args, **kwargs), len(self.source.data['x'])))
        else:
            if param.colname not in self.cache.columns:
                if self.is_empty:
                    param.update(pd.Series(dtype=object), *args, **kwargs)
                    return
                self.cache[param.colname] = self.sample.get_value(colname=param.colname)
            self.source.data[param.attr] = param.update(self.cache[param.colname].copy(), *args, **kwargs)
        if self.component is not None:
            self.component.apply_selected_obs_names()


from concurrent.futures import ThreadPoolExecutor
from functools import partial
from time import perf_counter

# Viewport loads of every layer and session (DynamicSourceManager.fetch_data_async)
FETCH_EXECUTOR = ThreadPoolExecutor(max_workers=settings.VIEWPORT_FETCH_WORKERS, thread_name_prefix="viewport-fetch")


class BokehSink:
    """LayerRows sink writing to a Bokeh ColumnDataSource."""

    def __init__(self, source):
        self.source = source

    @property
    def data(self):
        return self.source.data

    def stream(self, columns):
        self.source.stream(columns)

    def replace(self, columns):
        self.source.data = columns

    def set(self, name, values):
        self.source.data[name] = values


class DynamicSourceManager(LayerRows):
    """Lazily loaded layers (bins, cells) of the Bokeh UI: core row logic plus a
    ColumnDataSource, background loads on the document's loop and the layer
    component's hooks (selection, tooltips, visibility masks, hover names)."""

    def __init__(self, sample):
        self.component = None
        self._queued_doc = None
        self._context_callbacks = []
        self.source = ColumnDataSource(data={'x': [], 'y': []})
        super().__init__(sample, sink=BokehSink(self.source))

    def on_context_change(self, callback):
        self._context_callbacks.append(callback)

    def _notify_context_change(self):
        for callback in self._context_callbacks:
            callback()

    # -- LayerRows hooks -------------------------------------------------------

    def _on_context_change(self):
        self._notify_context_change()

    def _on_rows_changed(self):
        sync = getattr(self.component, 'sync_hover_names', None)
        if sync is not None:
            sync()

    def _show_not_sequenced(self):
        return self.component is None or self.component.show_not_sequenced

    def _after_style_change(self):
        if self.component is not None:
            self.component.apply_selected_obs_names()

    def _after_update(self):
        self._update_component()

    def _update_component(self):
        if self.component is not None:
            apply_constants = getattr(self.component, "apply_glyph_constants", None)
            if apply_constants is not None:
                apply_constants(self.glyph_constants)
            self.component.apply_selected_obs_names()

            self.component.update_hover_tooltips()
            cm = getattr(self.component, "clustering_manager", None)
            gm = getattr(self.component, "groups_manager", None)
            color_param = self.get_param('color')
            # El color puede ser un StaticParameter (sin colname).
            showing_groups = getattr(color_param, "colname", None) == settings.GROUPS_COL
            if showing_groups and gm is not None and getattr(gm, "_groups", None):
                gm._apply_visibility_mask()
            elif (not showing_groups and cm is not None and getattr(cm, "is_clustered", False)
                    and len(cm._active_clusters) < len(cm._cluster_labels)):
                cm._apply_visibility_mask()
            else:
                self.component._view = self.view_filter.copy()
                self.component.change_glyph_view()

    # -- background loads ---------------------------------------------------------

    def fetch_data_async(self, doc):
        """Like fetch_data, but Sample.fetch_data runs in a worker thread and the
        rows are added on a later tick of `doc`, so the map stays responsive. A
        result is dropped if the rows were cleared meanwhile (new level, colours...);
        calls made while a load runs are merged into one load that follows it."""
        if self.bbox is None or not self.active:
            return
        self._queued_doc = doc
        if not self._in_flight:
            self._start_async()

    def _start_async(self):
        doc, self._queued_doc = self._queued_doc, None
        self._in_flight = True
        request = self._request(list(self.bbox))
        generation = self._generation
        submitted = perf_counter()

        def work():
            started = perf_counter()
            data = self._fetch(request)
            return data, (started - submitted) * 1e3, (perf_counter() - started) * 1e3

        def done(future):                                  # in the worker thread
            try:
                doc.add_next_tick_callback(partial(self._finish_async, future, request, generation))
            except Exception as exc:                        # session closed meanwhile
                print(f"[fetch] result dropped: {exc!r}")

        FETCH_EXECUTOR.submit(work).add_done_callback(done)

    def _finish_async(self, future, request, generation):
        self._in_flight = False
        try:
            data, wait_ms, fetch_ms = future.result()
        except Exception as exc:
            print(f"[fetch] {type(self.sample).__name__} failed: {exc!r}")
            data = None
        if (data is not None and generation == self._generation
                and request["colnames"] == self.get_colnames()):
            self._add_fetched(data, request, wait_ms=wait_ms, fetch_ms=fetch_ms)
        if self._queued_doc is not None and self.bbox is not None:
            self._start_async()
