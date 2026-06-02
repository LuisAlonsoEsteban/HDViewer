

import pandas as pd
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
                if param.colname not in self.cache.columns:
                    colname_data = self.sample.fetch_data(colnames=param.colname, return_coords=False)
                    self.cache[param.colname] = colname_data[param.colname]
                self.source.data[param.attr] = param.get_values(self.cache[param.colname].copy())
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
        return list(set([param.colname for param in self.parameters if not param.is_static]))

    def fetch_data(self):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames())

        if len(self.get_colnames()) > 0:
            self.is_sequenced.clear()
            self.is_sequenced += (~fetched_data[self.get_colnames()[0]].isna()).values.tolist()

        for param in self.parameters:
            if not param.is_static:
                if param.colname not in self.cache.columns:
                    self.cache[param.colname] = fetched_data[param.colname]
                fetched_data[param.attr] = param.get_values(fetched_data[param.colname].copy())
            else:
                fetched_data[param.attr] = param.get_values()

        fetched_data.drop([param.colname for param in self.parameters if not param.is_static and not param.colname_hovered], axis=1, inplace=True)
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
            'obs_name': [],
            **attr_dict, 
            **hovered_dict
        }
        self.view_filter = []

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
        obs_names = pd.Series(self.source.data['obs_name'])
        return obs_names.isna().to_numpy()
        
    def add_data(self, bbox=None):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames(), bbox=bbox)
        new_data = fetched_data[~fetched_data.index.isin(self.indexes)]

        self.indexes += new_data.index.tolist()
        self.view_filter += [True] * new_data.shape[0]

        if len(self.get_colnames()) > 0:
            new_is_sequenced = (~new_data[self.get_colnames()[0]].isna()).values.tolist()
            self.is_sequenced += new_is_sequenced

        for param in self.parameters:
            if not param.is_static:

                # if param.colname in self.cache.columns:
                #     self.cache[param.colname] = pd.concat((self.cache[param.colname], new_data[param.colname]), axis=0) # ValueError: cannot reindex on an axis with duplicate labels
                # else:
                #     self.cache[param.colname] = new_data[param.colname]

                new_data[param.attr] = param.get_values(new_data[param.colname].copy())
            else:
                new_data[param.attr] = param.get_values(is_sequenced_values=new_is_sequenced)

        #new_data = new_data.drop(self.get_colnames(), axis=1, inplace=False)
        cols_to_drop = [param.colname for param in self.parameters 
                if not param.is_static and not param.colname_hovered]
        new_data = new_data.drop(cols_to_drop, axis=1, inplace=False)
        
        self.source.stream(new_data.to_dict(orient='list'))

        if self.component is not None:
            self.component.apply_selected_obs_names()

            self.component._view = self.view_filter
            self.component.change_glyph_view()


    def fetch_data(self):
        self.add_data(bbox=self.bbox)

    def remove_indcs(self, rem_indcs):
        for key in self.source.data.keys():
            self.source.data[key] = [value for idx, value in zip(self.indexes, self.source.data[key]) if idx not in rem_indcs]
        
        self.view_filter = [value for idx, value in zip(self.indexes, self.view_filter) if idx not in rem_indcs]
        self.component._view = self.view_filter

        self.indexes = [ind for ind in self.indexes if ind not in rem_indcs]

    def remove_data(self, bbox):
        fetched_data = self.sample.fetch_data(bbox=bbox)
        removed_indexes = fetched_data['obs_name'].values.tolist()
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
        distances = [DynamicSourceManager.distance((x,y), (_x,_y)) for _x, _y in zip(self.source.data['x'], self.source.data['y'])]
        removed_indexes = [i for d, i in zip(distances, self.indexes) if d > distance]
        print(f"Removing {len(removed_indexes)} bins")
        self.remove_indcs(removed_indexes)
