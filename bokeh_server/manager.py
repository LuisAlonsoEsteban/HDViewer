

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
        print(f"[parameters lista] {[(p.colname, id(p)) for p in self.parameters if not p.is_static]}")
        existing_param = self.get_param(param.attr)
        if existing_param is not None:
            self.parameters.remove(self.get_param(param.attr))
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
                        print(f"[obs_names] primeros 3: {obs_names.dropna().head(3).tolist()}")
                        print(f"[obs_series index] primeros 3: {obs_series.index[:3].tolist()}")
                        
                        aligned = obs_names.map(obs_series)
                        print(f"[aligned] no-NaN: {aligned.dropna().unique()[:5]}")
                        
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

        del self.source.data[old_colname]
        
    def get_colnames(self):
        return list({param.colname for param in self.parameters if not param.is_static})

    def fetch_data(self):
        print(f"[fetch_data] EJECUTÁNDOSE")
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
            'obs_name': [],
            **attr_dict, 
            **hovered_dict
        }
        self.view_filter = np.empty(0, dtype=bool)

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
        if self.get_param(param.attr) is None:
            self.parameters.append(param)

            if param.is_static:
                param.set_is_sequenced(self.is_sequenced)
            else:
                param.min_counts, param.max_counts = self.sample.get_range(param.colname)

            if not param.is_static or self.is_empty:
                self.clear_sources()
            else:
                self.source.data[param.attr] = param.get_values()
        else:
            raise ValueError(f"A parameter for attribute '{param.attr}' is already present.")
        
    def get_nan_vector(self):
        colnames = self.get_colnames()
        if len(colnames) > 0 and colnames[0] in self.source.data:
            return pd.Series(self.source.data[colnames[0]]).isna().to_numpy()
        obs_names = pd.Series(self.source.data['obs_name'])
        return obs_names.isna().to_numpy()
        
    def set_fetch_context(self, bbox=None, resolution_level=None, **fetch_kwargs):
        resolution_changed = resolution_level != self.resolution_level
        extra_changed = fetch_kwargs != self.fetch_kwargs

        self.bbox = bbox
        self.resolution_level = resolution_level
        self.fetch_kwargs = fetch_kwargs

        if resolution_changed or extra_changed:
            self.clear_sources()

        return resolution_changed or extra_changed

    def add_data(self, bbox=None):
        colnames = self.get_colnames()
        fetch_kwargs = dict(self.fetch_kwargs)
        if self.resolution_level is not None:
            fetch_kwargs["resolution_level"] = self.resolution_level
            
        fetched_data = self.sample.fetch_data(colnames=colnames, bbox=bbox, **fetch_kwargs)
        new_data = fetched_data[~fetched_data.index.isin(self.indexes)]

        self.indexes.extend(new_data.index.tolist())
        self.view_filter = np.concatenate((
            self.view_filter,
            np.ones(new_data.shape[0], dtype=bool),
        ))

        new_is_sequenced = []
        if len(colnames) > 0:
            new_is_sequenced = (~new_data[colnames[0]].isna()).to_numpy().tolist()
            self.is_sequenced += new_is_sequenced

        for param in self.parameters:
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
        
        self.source.stream(new_data.to_dict(orient='list'))

        if self.component is not None:
            self.component.apply_selected_obs_names()

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
        if self.resolution_level is not None:
            fetch_kwargs["resolution_level"] = self.resolution_level
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

    def remove_further_data(self, x, y, distance=1000):
        x_values = np.asarray(self.source.data['x'], dtype=float)
        y_values = np.asarray(self.source.data['y'], dtype=float)

        distances = np.hypot(x_values - x, y_values - y)
        indexes_array = np.asarray(self.indexes, dtype=object)
        
        removed_indexes = indexes_array[distances > distance].tolist()
        self.remove_indcs(removed_indexes)
