

import pandas as pd
from bokeh.models import ColumnDataSource

class SourceManager:

    def __init__(self, sample):
        self.sample = sample
        self.parameters = []
        self.is_empty = True

        self.is_sequenced = []

        self.source = ColumnDataSource(data={'x': [], 'y': []})
        self.cache = pd.DataFrame()

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
    
    def add_parameter(self, param):
        if self.get_param(param.attr) is None:
            self.parameters.append(param)
            if param.is_static:
                param.set_is_sequenced(self.is_sequenced)
            else:
                if self.sample.has_adata():
                    param.min_counts, param.max_counts = self.sample.get_range(param.colname)
                else:
                    raise ValueError("Non-static parameters cannot be added without an AnnData associated")
            
            if not self.is_empty:
                if param.is_static:
                    self.source.data[param.attr] = param.get_values()
                elif param.colname in self.cache.columns:
                    self.source.data[param.attr] = param.get_values(self.cache[param.colname].copy())
            else:
                self.reset()

        else:
            raise ValueError(f"A parameter for attribute '{param.attr}' is already present.")
        
    def update_colname(self, old_colname, new_colname):
        for param in self.get_params_by_colname(old_colname):
            param.colname = new_colname
            param.min_counts, param.max_counts = self.sample.get_range(param.colname)

            if param.colname not in self.cache.columns:
                colname_data = self.sample.fetch_data(colnames=param.colname, return_coords=False)
                self.cache[param.colname] = colname_data[param.colname]

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
        self.is_empty = False

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self.source.data[param.attr] = param.update(*args, **kwargs)
        else:
            self.source.data[param.attr] = param.update(self.cache[param.colname].copy(), *args, **kwargs)


import math

class DynamicSourceManager(SourceManager):


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.indexes = []
        self.bbox = None

    def clear_sources(self):
        self.indexes.clear()
        self.is_sequenced.clear()

        self.cache = self.cache[0:0]

        if len(self.parameters) > 0:
            attr_dict = {param.attr: [] for param in self.parameters}
        else:
            attr_dict = dict()
        self.source.data = {'x': [], 'y': [], **attr_dict}

    def update_colname(self, old_colname, new_colname):
        for param in self.get_params_by_colname(old_colname):
            param.colname = new_colname
            param.min_counts, param.max_counts = self.sample.get_range(param.colname)
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
        
    def add_data(self, bbox=None):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames(), bbox=bbox)
        new_data = fetched_data[~fetched_data.index.isin(self.indexes)]

        self.indexes += new_data.index.tolist()
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

        new_data = new_data.drop(self.get_colnames(), axis=1, inplace=False)
        self.source.stream(new_data.to_dict(orient='list'))

    def fetch_data(self):
        self.add_data(bbox=self.bbox)

    def remove_indcs(self, rem_indcs):
        for key in self.source.data.keys():
            self.source.data[key] = [value for idx, value in zip(self.indexes, self.source.data[key]) if idx not in rem_indcs]
        self.indexes = [ind for ind in self.indexes if ind not in rem_indcs]

    def remove_data(self, bbox):
        fetched_data = self.sample.fetch_data(bbox=bbox)
        removed_indexes = fetched_data.index.tolist()
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

