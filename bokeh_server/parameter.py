
class Parameter:

    def __init__(self, is_static):
        self.is_static = is_static


class StaticParameter(Parameter):

    is_sequenced = None

    def __init__(self, value, attr='color', not_sequenced_val="#b8b8b8"):
        super().__init__(is_static=True)

        self.attr = attr

        self.value = value
        self.not_sequenced_val = not_sequenced_val

    def set_is_sequenced(self, is_sequenced):
        self.is_sequenced = is_sequenced

    def get_values(self, is_sequenced_values:list=None):
        if is_sequenced_values:
            return [self.value if is_sequenced
                    else self.not_sequenced_val
                    for is_sequenced in is_sequenced_values]
        
        elif self.is_sequenced is not None and len(self.is_sequenced) > 0:
            return [self.value if is_sequenced
                    else self.not_sequenced_val
                    for is_sequenced in self.is_sequenced]
        
        else:
            return self.value
    
    def update(self, new_val=None, new_not_sequenced_val=None):
        if new_val is not None:
            self.value = new_val
        if new_not_sequenced_val is not None:
            self.new_not_sequenced_val = new_not_sequenced_val
            
        return self.get_values()


class DynamicParameter(Parameter):

    def __init__(
            self, 
            colname, 
            attr='color', 
            min_counts=0.0, 
            max_counts=None, 
            static_not_sequenced=True, 
            not_sequenced_val="#f5f5f5",
            colname_hovered=False
            ):
        
        super().__init__(is_static=False)

        self.colname = colname
        self.attr = attr

        self.min_counts = min_counts
        self.max_counts = max_counts

        self.static_not_sequenced = static_not_sequenced
        self.not_sequenced_val = not_sequenced_val

        self.colname_hovered = colname_hovered

    def handle_nans(self, series):
        not_sequenced = series.isna()
        series.loc[not_sequenced] = self.not_sequenced_val
        return series, ~not_sequenced
    
    def apply_filter(self, series, sequenced):
        pass

    def get_values(self, series):
        return self.apply_filter(*self.handle_nans(series))
    

import matplotlib.colors as mcolors

class Color(DynamicParameter):

    def __init__(self, c1, c2, *args, **kwargs):
        super().__init__(*args, **kwargs)
    #     self.cmap = Color.get_cmap(c1, c2)

    @staticmethod
    def get_cmap(color1, color2, name='custom_diverging'):
        return mcolors.LinearSegmentedColormap.from_list(
                name, 
                [color1, color2]
            )
    
    def apply_filter(self, values):
        norm = mcolors.Normalize(vmin=min(values), vmax=max(values))
        return self.cmap(norm(values))


class Scaler(DynamicParameter):

    def __init__(self, minv, maxv, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.minv = minv
        self.maxv = maxv

    def scale_values(self, values):
        range_val = self.max_counts - self.min_counts if self.max_counts != self.min_counts else 1
        scaled_values = [self.minv + (value - self.min_counts) * (self.maxv - self.minv) / range_val for value in values]
        return scaled_values

    def apply_filter(self, series, sequenced):
        if self.static_not_sequenced:
            sequenced_values = series[sequenced]

            if len(sequenced_values) > 0:
                scaled_values = self.scale_values(sequenced_values)
            else:
                scaled_values = []

            result_series = series.copy()
            result_series[sequenced] = scaled_values
            
            return result_series.tolist()
        else:
            return self.scale_values(series)
    
    def update(self, series, new_minv=None, new_maxv=None):
        if new_minv is not None:
            self.minv = new_minv if new_minv < self.maxv else self.maxv
        if new_maxv is not None:
            self.maxv = new_maxv if new_maxv > self.minv else self.minv
        return self.get_values(series)
