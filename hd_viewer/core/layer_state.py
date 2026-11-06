"""Rows of a lazily loaded map layer, without any plotting library.

`LayerRows` decides what a viewer holds for one layer: which rows are loaded for
the current view and display level, their encoded style columns (alpha, colour
codes), the raw values kept back for restyling, row names, visibility and the
per-level glyph constants. It writes every change through a *sink*:

- the HTTP API uses `RecordingSink`, whose recorded operations are sent to the
  browser and applied there to the BokehJS data source.

Subclasses customise the hooks (`_show_not_sequenced`, `_after_update`,
`_after_style_change`, `_on_rows_changed`, `_on_context_change`).
"""
import numpy as np
import pandas as pd

from .. import settings


def _compact(array):
    # Screen coordinates and styles need no double precision: half the bytes
    return array.astype(np.float32) if array.dtype == np.float64 else array


def source_columns(frame):
    """Keep numeric columns binary-serializable; polygon/string columns stay lists.

    Own the arrays so source updates cannot mutate a sample's cached DataFrame.
    Object columns deliberately stay one-dimensional lists, including equal-size
    polygon rings, which NumPy would otherwise expand into a rectangular array.
    """
    return {name: _compact(column.to_numpy(copy=True)) if column.dtype.kind in 'biuf'
            else column.tolist() for name, column in frame.items()}


def parameter_column(values, size):
    if np.isscalar(values):
        values = [values] * size
    array = np.asarray(values)
    return _compact(array) if array.dtype.kind in 'biuf' else list(values)


class RecordingSink:
    """In-memory columns plus the operations applied since the last `take()`:
    ('append', columns), ('replace', columns) and ('set', name, values)."""

    def __init__(self):
        self.data = {}
        self.ops = []

    def stream(self, columns):
        for name, values in columns.items():
            current = self.data.get(name)
            if current is None or len(current) == 0:
                self.data[name] = values
            elif isinstance(current, np.ndarray) and isinstance(values, np.ndarray):
                self.data[name] = np.concatenate((current, values))
            else:
                self.data[name] = list(current) + list(values)
        self.ops.append(('append', columns))

    def replace(self, columns):
        self.data = dict(columns)
        self.ops.append(('replace', dict(columns)))

    def set(self, name, values):
        self.data[name] = values
        self.ops.append(('set', name, values))

    def take(self):
        ops, self.ops = self.ops, []
        return ops


class ParameterSet:
    """Style parameters of a layer (shared by every source manager)."""

    def get_param(self, attr):
        for param in self.parameters:
            if param.attr == attr:
                return param
        return None

    def get_params_by_colname(self, colname):
        return [param for param in self.parameters if not param.is_static and param.colname == colname]

    def get_colnames(self):
        """Columns read by the parameters, in parameter order (deterministic)."""
        return list(dict.fromkeys(param.colname for param in self.parameters if not param.is_static))

    def _sequenced_column(self, colnames):
        """Column telling sequenced rows apart: the opacity's (a gene value,
        missing only where nothing was sequenced), not e.g. group labels, which
        are missing for every row outside the groups."""
        alpha = self.get_param('alpha')
        if alpha is not None and not alpha.is_static and alpha.colname in colnames:
            return alpha.colname
        return colnames[0]

    def _gene_parameters(self):
        adata = self.sample.adata
        return [param for param in self.parameters
                if not param.is_static and adata is not None and param.colname in adata.var_names]

    def _set_column(self, name, values):
        """Assign a column unless it is unchanged (e.g. a colour picker edit only
        changes the palette of a colour-coded column)."""
        current = self._data.get(name)
        if (isinstance(current, np.ndarray) and isinstance(values, np.ndarray)
                and current.shape == values.shape and np.array_equal(current, values)):
            return
        self._set(name, values)

    def colors(self, attr='color'):
        """Displayed hex colours of a colour-coded attribute, row by row."""
        param = self.get_param(attr)
        values = self._data.get(attr, [])
        if getattr(param, 'palette', None) is None:
            return list(values)
        return param.decode(values)


class LayerRows(ParameterSet):
    """Loaded rows of one layer view; see the module docstring."""

    GLYPH_CONSTANTS = ('width', 'height', 'angle')

    def __init__(self, sample, sink=None):
        self.sample = sample
        self.sink = RecordingSink() if sink is None else sink
        self.parameters = []
        self.is_sequenced = []
        self.cache = pd.DataFrame()
        # Bumped whenever the loaded rows are dropped (new level, colours...): a
        # background load started before that is stale and must not be added
        self._generation = 0
        self._in_flight = False
        # Raw, non-hovered values remain on the server for slider/style changes.
        # Only columns used by active parameters are retained.
        self._value_cache = {}
        self.indexes = []
        self.bbox = None
        self.fetch_kwargs = {}
        self.resolution_level = None
        self.fixed_aggregation_factor = None
        # Constant per display level (bins): sent once as glyph values, not per row
        self.glyph_constants = {}
        # Per-row geometry columns (x, y; spots also radius), kept when the rows
        # are cleared so the source always matches the glyph
        self.coordinate_columns = ('x', 'y')
        self._last_fetch_rows = 0
        # A hidden layer keeps nothing loaded and fetches nothing until shown
        self.active = True
        # Row names stay on the server (hover asks for them): positions of
        # `_obs_names` match the source rows; `row_version` changes with them
        self._obs_names = np.empty(0, dtype=object)
        self.row_version = 0
        self.clear_sources()

    # -- sink and hooks --------------------------------------------------

    @property
    def _data(self):
        return self.sink.data

    def _set(self, name, values):
        self.sink.set(name, values)

    def _show_not_sequenced(self):
        return True

    def _after_update(self):
        """Rows or styles changed (e.g. reapply visibility masks)."""

    def _after_style_change(self):
        """A style parameter changed without touching rows."""

    def _on_rows_changed(self):
        """Row positions changed: cached row names are stale."""

    def _on_context_change(self):
        """View, level or loaded rows were reset."""

    # -- state -------------------------------------------------------------

    @property
    def loaded_obs_names(self):
        return self._obs_names

    @property
    def n_geometries(self):
        return len(self.indexes)

    def _rows_changed(self):
        self.row_version += 1
        self._on_rows_changed()

    def set_active(self, active):
        """Hidden layers drop their rows (freeing the browser) and skip loads;
        the caller reloads the viewport when the layer is shown again."""
        if active == self.active:
            return
        self.active = active
        if not active:
            self.clear_sources()

    def needs_pruning(self):
        """Drop far rows only once the loaded rows exceed a budget: every prune
        resends all retained rows to the browser."""
        return len(self.indexes) > settings.PRUNE_ROW_BUDGET * max(self._last_fetch_rows, 1)

    def clear_sources(self):
        self._generation += 1
        self.indexes.clear()
        self.is_sequenced.clear()

        self.cache = self.cache[0:0]
        self._value_cache.clear()

        attr_dict = {param.attr: [] for param in self.parameters}
        hovered_dict = {param.colname: [] for param in self.parameters
                        if not param.is_static and param.colname_hovered}
        self.sink.replace({**{name: [] for name in self.coordinate_columns}, **attr_dict, **hovered_dict})
        self._obs_names = np.empty(0, dtype=object)
        self._rows_changed()
        self.view_filter = np.empty(0, dtype=bool)
        self._on_context_change()

    # -- parameters ----------------------------------------------------------

    def update_colname(self, old_colname, new_colname, attrs=None):
        """Point the parameters reading `old_colname` at `new_colname`; with `attrs`,
        only those attributes (other parameters may read the same column)."""
        previous_level = self.effective_resolution_level
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

        if not self._refresh_loaded_parameters(previous_level, refetch=True):
            self.clear_sources()
            self.fetch_data()

    def refresh_values(self):
        """Re-read gene values after the sample's expression scale changed,
        keeping the loaded geometry when possible."""
        params = self._gene_parameters()
        if not params:
            return
        for param in params:
            if hasattr(param, "min_counts"):
                param.min_counts, param.max_counts = self.sample.get_range(param.colname)
        if self._refresh_loaded_parameters(self.effective_resolution_level, refetch=True):
            return
        self.clear_sources()
        if self.bbox is not None:
            self.fetch_data()

    def add_parameter(self, param):
        # Replace any parameter already bound to this attribute (e.g. clustering
        # swaps the static color for the cluster column)
        previous_level = self.effective_resolution_level
        existing_param = self.get_param(param.attr)
        if existing_param is not None:
            self.parameters.remove(existing_param)
        self.parameters.append(param)

        if param.is_static:
            param.set_is_sequenced(self.is_sequenced)
        elif self.sample.is_continous(param.colname):
            param.min_counts, param.max_counts = self.sample.get_range(param.colname)

        if self._refresh_loaded_parameters(previous_level, refetch=not param.is_static):
            return

        # Changed resolution or incomplete row coverage requires a full reload.
        self.clear_sources()
        if self.bbox is not None:
            self.fetch_data()

    def remove_parameter(self, attr):
        param = self.get_param(attr)
        if param is None:
            return
        previous_level = self.effective_resolution_level
        self.parameters.remove(param)
        if self._refresh_loaded_parameters(previous_level):
            return
        self.clear_sources()
        if self.bbox is not None:
            self.fetch_data()

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self._set_column(param.attr, parameter_column(param.update(*args, **kwargs), len(self.indexes)))
        elif param.colname in self._data or param.colname in self._value_cache:
            values = self._data.get(param.colname, self._value_cache.get(param.colname))
            series = pd.Series(values)
            self._set(param.attr, parameter_column(param.update(series, *args, **kwargs), len(self.indexes)))
        else:
            param.update(pd.Series(dtype=object), *args, **kwargs)
            self.clear_sources()
            self.fetch_data()
            return
        self._after_style_change()

    def _refresh_loaded_parameters(self, previous_level, refetch=False):
        """Update values/style without retransmitting the loaded geometry.

        A pan may have retained rows outside the latest fetch box. If a new
        column cannot cover every loaded row, use the ordinary reload path.
        In-flight requests and resolution changes also use that path.
        """
        if not self.active:
            return True                 # nothing loaded; shown again → full reload
        if (not self.indexes or self._in_flight
                or previous_level != self.effective_resolution_level):
            return False
        colnames = self.get_colnames()
        data = self._data
        if refetch:
            request = self._request(self.bbox)
            fetched = self._fetch(dict(request, return_coords=False))
            if not fetched.index.is_unique:
                return False
            positions = fetched.index.get_indexer(self.indexes)
            if (positions < 0).any():
                return False
            values = {name: fetched[name].to_numpy()[positions] for name in colnames}
        else:
            values = {name: data.get(name, self._value_cache.get(name)) for name in colnames}
            if any(value is None for value in values.values()):
                return False

        if colnames:
            self.is_sequenced[:] = (~pd.isna(values[self._sequenced_column(colnames)])).tolist()
        updates = {}
        for param in self.parameters:
            if not refetch and not param.is_static and param.attr in data:
                continue
            if param.is_static:
                result = param.get_values(is_sequenced_values=self.is_sequenced)
            else:
                series = pd.Series(values[param.colname])
                if param.colname.startswith('Group'):
                    series = series.fillna('No asignado').astype(str)
                result = param.get_values(series)
                if param.colname_hovered:
                    updates[param.colname] = parameter_column(values[param.colname], len(self.indexes))
            # Static edits only need to send their own column. Other style
            # columns have not changed if their raw values are unchanged.
            if refetch or param.attr not in data or param.is_static:
                updates[param.attr] = parameter_column(result, len(self.indexes))
        hovered = {param.colname for param in self.parameters
                   if not param.is_static and param.colname_hovered}
        self._value_cache = {name: np.asarray(value) for name, value in values.items() if name not in hovered}
        keep = {*self.coordinate_columns, 'obs_name', *updates, *(param.attr for param in self.parameters), *hovered}
        if set(data) - keep:
            self.sink.replace({**{key: value for key, value in data.items() if key in keep}, **updates})
        else:
            for name, value in updates.items():
                self._set_column(name, value)
        self.view_filter = np.ones(len(self.indexes), dtype=bool)
        if not self._show_not_sequenced() and len(self.is_sequenced) == len(self.indexes):
            self.view_filter = np.asarray(self.is_sequenced, dtype=bool)
        self._generation += 1
        self._on_context_change()
        self._after_update()
        return True

    def get_nan_vector(self):
        names = pd.Series(self.loaded_obs_names, dtype=object)
        sequenced = names.isin(self.sample.adata.obs_names)
        # Aggregated bins use display names rather than native observation names.
        prefix = getattr(self.sample, 'aggregate_obs_prefix', None)
        if prefix is not None:
            summaries = names.dropna().astype(str)
            summaries = summaries[summaries.str.startswith(prefix)]
            suffixes = summaries.str.rsplit(':', n=1).str[-1]
            factors = pd.to_numeric(suffixes[suffixes.str.isdigit()].unique()).astype(int)
            for factor in factors:
                if factor in self.sample.AGGREGATION_FACTORS:
                    sequenced |= names.isin(self.sample.aggregation_names(factor).dropna())
        return ~sequenced.to_numpy(dtype=bool)

    # -- view context --------------------------------------------------------

    @property
    def minimum_aggregation_factor(self):
        """Cluster columns set a floor; larger display bins remain available."""
        get_forced_factor = getattr(self.sample, "get_forced_aggregation_factor", None)
        return get_forced_factor(self.get_colnames()) if get_forced_factor is not None else None

    @property
    def effective_resolution_level(self):
        requested = (self.sample.get_resolution_level_for_factor(self.fixed_aggregation_factor)
                     if self.fixed_aggregation_factor is not None else self.resolution_level)
        minimum = self.minimum_aggregation_factor
        if minimum is None:
            return requested
        get_factor = getattr(self.sample, 'get_aggregation_factor', lambda level: 1)
        return self.sample.get_resolution_level_for_factor(max(minimum, get_factor(requested)))

    def set_aggregation_factor(self, factor):
        """None follows zoom; 1 means native; other supported factors are fixed."""
        if factor is not None and factor not in self.sample.AGGREGATION_FACTORS:
            raise ValueError(f"Unsupported bin aggregation factor: {factor}")
        minimum = self.minimum_aggregation_factor
        if factor is not None and minimum is not None and factor < minimum:
            raise ValueError(f"Displayed cluster columns require an aggregation factor of at least {minimum}.")
        previous_level = self.effective_resolution_level
        self.fixed_aggregation_factor = factor
        if self.effective_resolution_level != previous_level:
            self.clear_sources()
        self._on_context_change()

    def set_fetch_context(self, bbox=None, resolution_level=None, **fetch_kwargs):
        previous_level = self.effective_resolution_level
        extra_changed = fetch_kwargs != self.fetch_kwargs

        self.bbox = bbox
        self.resolution_level = resolution_level
        self.fetch_kwargs = fetch_kwargs
        resolution_changed = self.effective_resolution_level != previous_level

        if resolution_changed or extra_changed:
            self.clear_sources()

        self._on_context_change()
        return resolution_changed or extra_changed

    # -- loading ---------------------------------------------------------------

    def _request(self, bbox):
        """Arguments of Sample.fetch_data for the current context."""
        fetch_kwargs = dict(self.fetch_kwargs)
        if self.effective_resolution_level is not None:
            fetch_kwargs["resolution_level"] = self.effective_resolution_level
        return dict(colnames=self.get_colnames(), bbox=bbox, **fetch_kwargs)

    def _fetch(self, request):
        with self.sample.fetch_lock:
            return self.sample.fetch_data(**request)

    def add_data(self, bbox=None):
        request = self._request(bbox)
        self._add_fetched(self._fetch(request), request)

    def fetch_data(self):
        if self.bbox is None or not self.active:
            return
        self.add_data(bbox=self.bbox)

    def _add_fetched(self, fetched_data, request, wait_ms=None, fetch_ms=None):
        """Add the rows of `fetched_data` that are not loaded yet."""
        if not self.active:
            return
        colnames = request["colnames"]
        self._last_fetch_rows = fetched_data.shape[0]
        new_data = fetched_data[~fetched_data.index.isin(self.indexes)].copy()

        if new_data.shape[0] == 0:
            return

        constants = [name for name in self.GLYPH_CONSTANTS if name in new_data.columns]
        for name in constants:
            self.glyph_constants[name] = float(new_data[name].iloc[0])
        new_data = new_data.drop(columns=constants)

        geometry = [name for name in new_data.columns if name not in colnames and name != 'obs_name'
                    and name not in {param.attr for param in self.parameters}]
        if set(geometry) - set(self.coordinate_columns):
            self.coordinate_columns = tuple(dict.fromkeys((*self.coordinate_columns, *geometry)))

        new_is_sequenced = []
        if len(colnames) > 0:
            new_is_sequenced = (~new_data[self._sequenced_column(colnames)].isna()).to_numpy().tolist()

        for param in self.parameters:
            if not param.is_static:
                vector_datos = new_data[param.colname].copy()
                if param.colname.startswith("Group"):
                    vector_datos = vector_datos.fillna("No asignado").astype(str)

                new_data[param.attr] = param.get_values(vector_datos)
            else:
                new_data[param.attr] = param.get_values(is_sequenced_values=new_is_sequenced)

        cols_to_drop = [param.colname for param in self.parameters
                        if not param.is_static and not param.colname_hovered
                        and not any(not other.is_static and other.colname == param.colname
                                    and other.colname_hovered for other in self.parameters)]
        for name in set(cols_to_drop):
            values = new_data[name].to_numpy(copy=True)
            previous = self._value_cache.get(name)
            self._value_cache[name] = values if previous is None else np.concatenate((previous, values))
        new_data = new_data.drop(cols_to_drop, axis=1, inplace=False)

        names = new_data['obs_name'].to_numpy(dtype=object)
        new_data = new_data.drop(columns=['obs_name'])
        existing_cols = set(self._data.keys())
        present_cols = set(new_data.columns)

        for missing in existing_cols - present_cols:
            new_data[missing] = np.nan

        extra = present_cols - existing_cols - set(self.coordinate_columns if existing_cols else ())
        if extra:
            new_data = new_data.drop(columns=list(extra), inplace=False)

        self.sink.stream(source_columns(new_data))
        self._obs_names = np.concatenate((self._obs_names, names))

        self.indexes.extend(new_data.index.tolist())
        new_visible = np.ones(new_data.shape[0], dtype=bool)
        if not self._show_not_sequenced() and len(new_is_sequenced) == new_data.shape[0]:
            new_visible = np.asarray(new_is_sequenced, dtype=bool)
        self.view_filter = np.concatenate((self.view_filter, new_visible))
        self.is_sequenced += new_is_sequenced

        self._after_update()

    # -- pruning -----------------------------------------------------------------

    def remove_indcs(self, rem_indcs):
        if len(rem_indcs) == 0 or len(self.indexes) == 0:
            return

        keep_mask = ~pd.Index(self.indexes).isin(set(rem_indcs))

        self.sink.replace({key: (values[keep_mask] if isinstance(values, np.ndarray)
                                 else [value for value, keep in zip(values, keep_mask) if keep])
                           for key, values in self._data.items()})

        self._value_cache = {name: values[keep_mask] for name, values in self._value_cache.items()}
        if len(self.cache) == len(keep_mask):
            self.cache = self.cache.iloc[keep_mask].reset_index(drop=True)
        if len(self.is_sequenced) == len(keep_mask):
            self.is_sequenced[:] = np.asarray(self.is_sequenced, dtype=bool)[keep_mask].tolist()

        self.view_filter = self.view_filter[keep_mask]
        self._obs_names = self._obs_names[keep_mask]
        self._rows_changed()

        self.indexes = np.asarray(self.indexes, dtype=object)[keep_mask].tolist()
        self._after_update()

    def remove_data(self, bbox):
        fetch_kwargs = dict(self.fetch_kwargs)
        if self.effective_resolution_level is not None:
            fetch_kwargs["resolution_level"] = self.effective_resolution_level
        fetched_data = self._fetch(dict(bbox=bbox, **fetch_kwargs))
        self.remove_indcs(fetched_data['obs_name'].to_numpy())

    @staticmethod
    def _centers(values):
        if len(values) > 0 and np.ndim(values[0]) > 0:
            return np.array([np.mean(v) if len(v) else np.nan for v in values], dtype=float)
        return np.asarray(values, dtype=float)

    def remove_further_data(self, x, y, distance=1000):
        distances = np.hypot(self._centers(self._data['x']) - x, self._centers(self._data['y']) - y)
        removed_indexes = np.asarray(self.indexes, dtype=object)[distances > distance].tolist()
        self.remove_indcs(removed_indexes)
