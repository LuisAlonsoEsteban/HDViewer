

import pandas as pd
import numpy as np

from bokeh.models import ColumnDataSource

class SourceManager:

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

    def reset(self):
        self.source = ColumnDataSource(data={'x': [], 'y': []})
        for param in self.parameters:
            self.source.data[param.attr] = []
        self.is_empty = True

    def get_param(self, attr):
        for param in self.parameters:
            if param.attr == attr:
                return param
        return None
    
    def get_params_by_colname(self, colname):
        params = []
        for param in self.parameters:
            if not param.is_static and param.colname == colname:
                params.append(param)
        return params

    def get_nan_vector(self):
        dummy_vector = self.sample.get_value(colname=self.get_colnames()[0])
        return dummy_vector.isna().to_numpy()
    
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
                        obs_names = pd.Series(self.source.data['obs_name'])
                        
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
                current_data = dict(self.source.data)
                current_data[param.attr] = valores
                self.source.data = current_data
        else:
            self.reset()

        
    def update_colname(self, old_colname, new_colname):
        for param in self.get_params_by_colname(old_colname):
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
        
    def get_colnames(self):
        return list({param.colname for param in self.parameters if not param.is_static})

    def fetch_data(self):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames())
        colnames = self.get_colnames()
        if len(colnames) > 0:
            self.is_sequenced.clear()
            self.is_sequenced += (~fetched_data[colnames[0]].isna()).to_numpy().tolist()

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
        self.source.data = fetched_data.to_dict(orient='list')
        if self.component is not None:
            self.component.apply_selected_obs_names()
        self.is_empty = False

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self.source.data[param.attr] = param.update(*args, **kwargs)
        else:
            if param.colname not in self.cache.columns:
                if self.is_empty:
                    param.update(pd.Series(dtype=object), *args, **kwargs)
                    return
                self.cache[param.colname] = self.sample.get_value(colname=param.colname)
            self.source.data[param.attr] = param.update(self.cache[param.colname].copy(), *args, **kwargs)
        if self.component is not None:
            self.component.apply_selected_obs_names()


import math

class DynamicSourceManager(SourceManager):


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.indexes = []
        self.bbox = None
        self.fetch_kwargs = {}
        self.resolution_level = None
        self.fixed_aggregation_factor = None
        self._context_callbacks = []
        self.clear_sources()

    def on_context_change(self, callback):
        self._context_callbacks.append(callback)

    def _notify_context_change(self):
        for callback in self._context_callbacks:
            callback()

    @property
    def n_geometries(self):
        return len(self.indexes)

    def clear_sources(self):
        self.indexes.clear()
        self.is_sequenced.clear()

        self.cache = self.cache[0:0]

        if len(self.parameters) > 0:
            attr_dict = {param.attr: [] for param in self.parameters}
            hovered_dict = {param.colname: [] for param in self.parameters 
                        if not param.is_static and param.colname_hovered}
        else:
            attr_dict = {}
            hovered_dict = {}
        
        self.source.data = {
            'x': [], 
            'y': [], 
            'width': [],
            'height': [],
            **({'angle': [], 'xs': [], 'ys': []} if hasattr(self.sample, 'native_angle') else {}),
            'obs_name': [],
            **attr_dict, 
            **hovered_dict
        }
        self.view_filter = np.empty(0, dtype=bool)
        self._notify_context_change()

    def update_colname(self, old_colname, new_colname):
        for param in self.get_params_by_colname(old_colname):
            param.colname = new_colname

            if hasattr(param, "min_counts"):
                min_val, max_val = self.sample.get_range(param.colname)
                if isinstance(min_val, str):
                    param.min_counts = 0
                    param.max_counts = 1
                else:
                    param.min_counts = min_val
                    param.max_counts = max_val

        self.clear_sources()
        self.fetch_data()

    def add_parameter(self, param):
        # Replace any parameter already bound to this attribute (e.g. clustering
        # swaps the static color for the cluster column)
        existing_param = self.get_param(param.attr)
        if existing_param is not None:
            self.parameters.remove(existing_param)
        self.parameters.append(param)

        if param.is_static:
            param.set_is_sequenced(self.is_sequenced)
        elif self.sample.is_continous(param.colname):
            param.min_counts, param.max_counts = self.sample.get_range(param.colname)

        if not param.is_static or self.is_empty:
            self.clear_sources()
            if existing_param is not None and self.bbox is not None:
                self.fetch_data()
        else:
            self.source.data[param.attr] = param.get_values()

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self.source.data[param.attr] = param.update(*args, **kwargs)
        elif param.colname in self.source.data:
            # The cache is not filled with lazy loading: rescale the loaded values
            series = pd.Series(self.source.data[param.colname])
            self.source.data[param.attr] = param.update(series, *args, **kwargs)
        else:
            param.update(pd.Series(dtype=object), *args, **kwargs)
            self.clear_sources()
            self.fetch_data()
            return

        if self.component is not None:
            self.component.apply_selected_obs_names()

    def get_nan_vector(self):
        colnames = self.get_colnames()
        if len(colnames) > 0 and colnames[0] in self.source.data:
            return pd.Series(self.source.data[colnames[0]]).isna().to_numpy()
        obs_names = pd.Series(self.source.data['obs_name'])
        return obs_names.isna().to_numpy()
        
    @property
    def effective_resolution_level(self):
        """Cluster locking takes precedence over a fixed level, then zoom."""
        get_forced_factor = getattr(self.sample, "get_forced_aggregation_factor", None)
        if get_forced_factor is not None:
            factor = get_forced_factor(self.get_colnames())
            if factor is not None:
                return self.sample.get_resolution_level_for_factor(factor)
        if self.fixed_aggregation_factor is not None:
            return self.sample.get_resolution_level_for_factor(self.fixed_aggregation_factor)
        return self.resolution_level

    def set_aggregation_factor(self, factor):
        """None follows zoom; 1 means native; other supported factors are fixed."""
        if factor is not None and factor not in self.sample.AGGREGATION_FACTORS:
            raise ValueError(f"Unsupported bin aggregation factor: {factor}")
        previous_level = self.effective_resolution_level
        self.fixed_aggregation_factor = factor
        if self.effective_resolution_level != previous_level:
            self.clear_sources()
        self._notify_context_change()

    def set_fetch_context(self, bbox=None, resolution_level=None, **fetch_kwargs):
        previous_level = self.effective_resolution_level
        extra_changed = fetch_kwargs != self.fetch_kwargs

        self.bbox = bbox
        self.resolution_level = resolution_level
        self.fetch_kwargs = fetch_kwargs
        resolution_changed = self.effective_resolution_level != previous_level

        if resolution_changed or extra_changed:
            self.clear_sources()

        self._notify_context_change()
        return resolution_changed or extra_changed

    def add_data(self, bbox=None):
        colnames = self.get_colnames()
        fetch_kwargs = dict(self.fetch_kwargs)
        if self.effective_resolution_level is not None:
            fetch_kwargs["resolution_level"] = self.effective_resolution_level
            
        fetched_data = self.sample.fetch_data(colnames=colnames, bbox=bbox, **fetch_kwargs)
        new_data = fetched_data[~fetched_data.index.isin(self.indexes)]

        # Sin filas nuevas no hay nada que añadir (y get_values([]) de los
        # parámetros estáticos devolvería el vector completo, no uno vacío).
        if new_data.shape[0] == 0:
            return

        new_is_sequenced = []
        if len(colnames) > 0:
            new_is_sequenced = (~new_data[colnames[0]].isna()).to_numpy().tolist()

        for param in self.parameters:
            if not param.is_static:
                
                # ─── PARCHE DE RESCATE DE CLUSTERS ───────────────────────────
                # Si el parámetro que se va a pintar coincide con tus grupos personalizados,
                # nos aseguramos de mapear los posibles NaNs o valores vacíos 
                # al estado de 'No asignado' para que el mapeador de color los reconozca.
                vector_datos = new_data[param.colname].copy()
                if param.colname.startswith("Group"):
                    vector_datos = vector_datos.fillna("No asignado").astype(str)
                # ─────────────────────────────────────────────────────────────

                new_data[param.attr] = param.get_values(vector_datos)
            else:
                new_data[param.attr] = param.get_values(is_sequenced_values=new_is_sequenced)

        #new_data = new_data.drop(self.get_colnames(), axis=1, inplace=False)
        cols_to_drop = [param.colname for param in self.parameters 
                if not param.is_static and not param.colname_hovered]
        new_data = new_data.drop(cols_to_drop, axis=1, inplace=False)
        
        # Reconciliar columnas con el esquema del source (p.ej. cells no trae width/height)
        existing_cols = set(self.source.data.keys())
        present_cols  = set(new_data.columns)

        # Columnas que el source espera pero new_data no trae → rellenar con NaN
        for missing in existing_cols - present_cols:
            new_data[missing] = np.nan

        # Columnas extra que el source no conoce → descartar (evita el error inverso)
        extra = present_cols - existing_cols
        if extra:
            new_data = new_data.drop(columns=list(extra), inplace=False)

        self.source.stream(new_data.to_dict(orient='list'))

        # Actualizar índices solo tras un stream correcto: si algo falla antes,
        # self.indexes y source.data no se desincronizan (rompía remove_further_data).
        self.indexes.extend(new_data.index.tolist())
        new_visible = np.ones(new_data.shape[0], dtype=bool)
        if (self.component is not None and not self.component.show_not_sequenced
                and len(new_is_sequenced) == new_data.shape[0]):
            new_visible = np.asarray(new_is_sequenced, dtype=bool)
        self.view_filter = np.concatenate((self.view_filter, new_visible))
        self.is_sequenced += new_is_sequenced

        if self.component is not None:
            self.component.apply_selected_obs_names()

            self.component.update_hover_tooltips()
            cm = getattr(self.component, "clustering_manager", None)
            gm = getattr(self.component, "groups_manager", None)
            color_param = self.get_param('color')
            # El color puede ser un StaticParameter (sin colname).
            showing_groups = getattr(color_param, "colname", None) == "custom_clusters"
            if showing_groups and gm is not None and getattr(gm, "_clusters", None):
                gm._apply_visibility_mask()
            elif (not showing_groups and cm is not None and getattr(cm, "is_clustered", False)
                    and len(cm._active_clusters) < len(cm._cluster_labels)):
                cm._apply_visibility_mask()
            else:
                self.component._view = self.view_filter.tolist()
                self.component.change_glyph_view()


    def fetch_data(self):
        self.add_data(bbox=self.bbox)

    def remove_indcs(self, rem_indcs):
        if len(rem_indcs) == 0 or len(self.indexes) == 0:
            return

        rem_idx_set = set(rem_indcs)
        keep_mask = ~pd.Index(self.indexes).isin(rem_idx_set)

        new_source_data = {}
        for key, values in self.source.data.items():
            values_array = np.asarray(values, dtype=object)
            new_source_data[key] = values_array[keep_mask].tolist()
        self.source.data = new_source_data

        self.view_filter = self.view_filter[keep_mask]
        if self.component is not None:
            self.component._view = self.view_filter.tolist()

        indexes_array = np.asarray(self.indexes, dtype=object)
        self.indexes = indexes_array[keep_mask].tolist()

    def remove_data(self, bbox):
        fetch_kwargs = dict(self.fetch_kwargs)
        if self.effective_resolution_level is not None:
            fetch_kwargs["resolution_level"] = self.effective_resolution_level
        fetched_data = self.sample.fetch_data(bbox=bbox, **fetch_kwargs)
        removed_indexes = fetched_data['obs_name'].to_numpy()
        self.remove_indcs(removed_indexes)

    @staticmethod
    def distance(p1, p2):
        if isinstance(p2[0], tuple) or isinstance(p2[1], tuple):
            p2 = (
                sum(p2[0]) / len(p2[0]),
                sum(p2[1]) / len(p2[1])
            )
            
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    @staticmethod
    def _centers(values):
        """Una coordenada por fila. En cells cada fila es la lista de vértices
        del polígono (patches): se usa su media para que el array sea 1D."""
        if len(values) > 0 and np.ndim(values[0]) > 0:
            return np.array([np.mean(v) if len(v) else np.nan for v in values], dtype=float)
        return np.asarray(values, dtype=float)

    def remove_further_data(self, x, y, distance=1000):
        x_values = self._centers(self.source.data['x'])
        y_values = self._centers(self.source.data['y'])

        distances = np.hypot(x_values - x, y_values - y)
        indexes_array = np.asarray(self.indexes, dtype=object)
        
        removed_indexes = indexes_array[distances > distance].tolist()
        self.remove_indcs(removed_indexes)
