
import math
from time import sleep

from manager import DynamicSourceManager

class Loader:

    OPTIMAL_BINS_AT_WINDOW = 50

    def __init__(self, samples:dict):
        self.samples = samples
        self.optimal_resolution_sample = self.samples[max(self.samples.keys())]
        self.manager = DynamicSourceManager(self.optimal_resolution_sample) 

        self.x_range = None
        self.y_range = None

        self.is_visible = False 

    @property
    def pixels_per_bin(self):
        return self.optimal_resolution_sample.pixels_per_bin

    @property
    def window_size(self):
        return self.x_range.end - self.x_range.start, self.y_range.end - self.y_range.start
    
    @property
    def x(self):
        return (self.x_range.start + self.x_range.end) / 2
    
    @property
    def y(self):
        return (self.y_range.start + self.y_range.end) / 2
    
    @property
    def static_window(self):
        return self.window_size[0] * self.STATIC_RANGE, self.window_size[1] * self.STATIC_RANGE

    def set_optimal_sample(self):
        visible_bins = {}
        for bins_size, sample in self.samples.items():
            visible_bins[bins_size] = self.window_size / sample.pixels_per_bin
        distances = {bins_size: abs(vb - self.OPTIMAL_BINS_AT_WINDOW) for bins_size, vb in visible_bins.items()}
        if self.optimal_resolution_sample.pixels_per_bin != self.samples[min(distances.keys())].pixels_per_bin:
            print("Resolution switched")
            self.optimal_resolution_sample = self.samples[min(distances, key=distances.get)]
            self.manager.reset()

    def add_data(self):
        self.manager.add_data(bbox=[
            self.x_range.start, 
            self.x_range.end,
            self.y_range.start, 
            self.y_range.end
            ])
        
    def on_change_ranges(self):
        self.set_optimal_sample()
        self.add_data()
    


# class Loader:

#     STATIC_RANGE = 0.6
#     OFFSET = 0.5

#     OPTIMAL_BINS_AT_WINDOW = 50

#     def __init__(self, samples:dict):
#         self.samples = samples
#         self.optimal_resolution_sample = self.samples[max(self.samples.keys())]
#         self.manager = DynamicSourceManager(self.optimal_resolution_sample) 

#         self.x_range = None
#         self.y_range = None

#         self.is_visible = False 

#     @property
#     def window_size(self):
#         return self.x_range.end - self.x_range.start, self.y_range.end - self.y_range.start
    
#     @property
#     def x(self):
#         return (self.x_range.start + self.x_range.end) / 2
    
#     @property
#     def y(self):
#         return (self.y_range.start + self.y_range.end) / 2
    
#     @property
#     def static_window(self):
#         return self.window_size[0] * self.STATIC_RANGE, self.window_size[1] * self.STATIC_RANGE

#     def set_optimal_sample(self):
#         visible_bins = []
#         for sample in self.samples:
#             visible_bins.append(self.window_size / sample.pixels_per_bin)
#         distances = [abs(vb - self.OPTIMAL_BINS_AT_WINDOW) for vb in visible_bins]
#         self.optimal_resolution_sample = self.samples[distances.index(min(distances))]
#         self.manager.reset()
    

#     def get_bbox(self, x, y):
#         return [
#             x - (self.window_size / 2),
#             x + (self.window_size / 2),
#             y - (self.window_size / 2),
#             y + (self.window_size / 2)
#         ]
    
#     def offset_bbox(self, bbox, offset=0.0):
#         return [
#             bbox[0] - offset,
#             bbox[1] + offset,
#             bbox[2] - offset,
#             bbox[3] + offset
#         ]
    
#     def get_offset_bbox(self, x, y):
#         return self.offset_bbox(self.get_bbox(x, y), offset=self.window_size*self.OFFSET)

#     def fetch_data(self, x, y):
#         self.manager.add_data(bbox=self.get_offset_bbox(x, y))
#         return self.manager.source
    
#     def on_change_centre(self, x, y):
#         if (x - self.x) > self.static_window[0] and (y - self.y) > self.static_window[1]:
#             new_data = self.fetch_data(x, y)
#             self.x, self.y = x, y
#         # sleep(0.2)
#         return new_data

#     def on_change_zoom(self, rx, ry):
#         self.window_size = (rx[1] - rx[0], ry[1] - ry[0])
#         if self.window_size[0] < self.VISIBLE_MAX_SIZE and \
#            self.window_size[1] < self.VISIBLE_MAX_SIZE:
#             self.is_visible = True
            
#     def on_fetch_data(self):
#         return
#         self.manager.remove_further_data(self.x, self.y, distance=self.VISIBLE_MAX_SIZE*10)