
import pandas as pd
from bokeh.models import ColumnDataSource

class StaticSourceManager:

    is_sequenced = []

    source = ColumnDataSource(data={'x': [], 'y': []})
    cache = pd.DataFrame()

    def __init__(self, sample):
        self.sample = sample
        self.parameters = []
        self.is_empty = True

    def reset(self):
        self.source.data['x'], self.source.data['y'] = [], []
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
            
            if not self.is_empty:
                if param.is_static:
                    self.source.data[param.attr] = param.get_values()
                elif param.colname in self.cache.columns:
                    self.source.data[param.attr] = param.get_values(self.cache[param.colname])
            else:
                self.reset()

        else:
            raise ValueError(f"A parameter for attribute '{param.attr}' is already present.")
        
    def update_colname(self, old_colname, new_colname):
        for param in self.get_params_by_colname(old_colname):
            param.colname = new_colname
            if param.colname not in self.cache.columns:
                colname_data = self.sample.fetch_data(colnames=param.colname, return_coords=False)
                self.cache[param.colname] = colname_data[param.colname]
            self.source.data[param.attr] = param.get_values(self.cache[param.colname])
        
    def get_colnames(self):
        return list(set([param.colname for param in self.parameters if not param.is_static]))

    def fetch_data(self):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames())

        if len(self.get_colnames()) > 0:
            self.is_sequenced += (~fetched_data[self.get_colnames()[0]].isna()).values.tolist()

        for param in self.parameters:
            if not param.is_static:
                if param.colname not in self.cache.columns:
                    self.cache[param.colname] = fetched_data[param.colname]
                fetched_data[param.attr] = param.get_values(fetched_data[param.colname])
            else:
                fetched_data[param.attr] = param.get_values()

        fetched_data.drop(self.get_colnames(), axis=1, inplace=True)
        self.source.data = fetched_data.to_dict(orient='list')
        self.is_empty = False

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self.source.data[param.attr] = param.update(*args, **kwargs)
        else:
            self.source.data[param.attr] = param.update(self.cache[param.colname].copy(), *args, **kwargs)
    



########################################################################################

class SourceManager:


    def __init__(self, sample):
        self.sample = sample
        self.parameters = []
        self.init_sources()

        self.indexes = []
        self.is_sequenced = []

    def init_sources(self):
        if len(self.parameters) > 0:
            attr_dict = {param.attr: [] for param in self.parameters}
        else:
            attr_dict = dict()
        self.colnames_data = pd.DataFrame()
        self.source = ColumnDataSource(data={'x': [], 'y': [], **attr_dict})

    def clear_sources(self):
        self.indexes.clear()
        self.is_sequenced.clear()

        self.colnames_data = self.colnames_data[0:0]

        if len(self.parameters) > 0:
            attr_dict = {param.attr: [] for param in self.parameters}
        else:
            attr_dict = dict()
        self.source.data = {'x': [], 'y': [], **attr_dict}

    def reset(self):
        self.indexes.clear()
        self.is_sequenced.clear()
        self.init_sources()

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
    
    def update_colname(self, old_colname, new_colname):
        print(old_colname, new_colname)
        for param in self.get_params_by_colname(old_colname):
            param.colname = new_colname
            # if param.colname not in self.colnames_data.columns:
            #     colname_data = self.sample.fetch_data(colnames=param.colname, bbox=self.bbox, return_coords=False)
            #     self.colnames_data[param.colname] = colname_data[self.indexes]
            # self.source.data[param.attr] = param.get_values(self.colnames_data[param.colname])
        self.clear_sources()
    
    def is_empty(self):
        return len(self.source.data['x']) == 0

    def add_parameter(self, param):
        if self.get_param(param.attr) is None:
            self.parameters.append(param)

            if param.is_static:
                param.set_is_sequenced(self.is_sequenced)

            if not param.is_static or self.is_empty():
                self.reset()
            else:
                self.source.data[param.attr] = param.get_values()
        else:
            raise ValueError(f"A parameter for attribute '{param.attr}' is already present.")
        
    def get_colnames(self):
        return list(set([param.colname for param in self.parameters if not param.is_static]))

    def add_data(self, bbox=None):
        fetched_data = self.sample.fetch_data(colnames=self.get_colnames(), bbox=bbox)
        new_data = fetched_data[~fetched_data.index.isin(self.indexes)]

        self.indexes += new_data.index.tolist()
        if len(self.get_colnames()) > 0:
            self.is_sequenced += (~new_data[self.get_colnames()[0]].isna()).values.tolist()

        for param in self.parameters:
            if not param.is_static:
                if param.colname not in self.colnames_data:
                    self.colnames_data[param.colname] = new_data[param.colname]
                new_data[param.attr] = param.get_values(new_data[param.colname])
            else:
                new_data[param.attr] = param.get_values()

        new_data = new_data.drop(self.get_colnames(), axis=1, inplace=False)
        self.source.stream(new_data.to_dict(orient='list'))

    def update_param(self, attr, *args, **kwargs):
        param = self.get_param(attr)
        if param.is_static:
            self.source.data[param.attr] = param.update(*args, **kwargs)
        else:
            print(self.colnames_data[param.colname])
            self.source.data[param.attr] = param.update(self.colnames_data[param.colname].copy(), *args, **kwargs)

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
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    def remove_further_data(self, x, y, distance=1000):
        distances = [SourceManager.distance((x,y), (_x,_y)) for _x, _y in zip(self.source.data['x'], self.source.data['y'])]
        removed_indexes = [i for d, i in zip(distances, self.indexes) if d > distance]
        self.remove_indcs(removed_indexes)
        

import math
from time import sleep

class Loader:

    STATIC_RANGE = 0.5
    OFFSET = 0.5

    VISIBLE_MAX_SIZE = 1000

    def __init__(self, sample, colnames=None):
        self.manager = SourceManager(sample, colnames=colnames) 

        self.window_size = None
        self.x = None
        self.y = None

        self.is_visible = False

    def get_bbox(self, x, y):
        return [
            x - self.window_size / 2,
            x + self.window_size / 2,
            y - self.window_size / 2,
            y + self.window_size / 2
        ]
    
    def offset_bbox(self, bbox, offset=0.0):
        return [
            bbox[0] - offset,
            bbox[1] + offset,
            bbox[2] - offset,
            bbox[3] + offset
        ]
    
    def get_offset_bbox(self, x, y):
        return self.offset_bbox(self.get_bbox(x, y), offset=self.window_size*self.OFFSET)
    
    @property
    def static_window(self):
        return self.window_size[0] * self.STATIC_RANGE, self.window_size[1] * self.STATIC_RANGE

    def fetch_data(self, x, y):
        self.manager.add_data(bbox=self.get_offset_bbox(x, y))
        return self.manager.source
    
    def on_fetch_data(self):
        self.manager.remove_further_data(self.x, self.y, distance=self.VISIBLE_MAX_SIZE*10)
    
    def on_change_centre(self, x, y):
        if (x - self.x) > self.static_window[0] and (y - self.y) > self.static_window[1]:
            new_data = self.fetch_data(x, y)
            self.x, self.y = x, y
        sleep(0.2)
        return new_data

    def on_change_zoom(self, rx, ry):
        self.window_size = (rx[1] - rx[0], ry[1] - ry[0])
        if self.window_size[0] < self.VISIBLE_MAX_SIZE and \
           self.window_size[1] < self.VISIBLE_MAX_SIZE:
            self.is_visible = True
        sleep(0.2)