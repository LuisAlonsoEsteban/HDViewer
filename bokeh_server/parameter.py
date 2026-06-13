
import numpy as np
import pandas as pd


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

    def _build_values(self, is_sequenced_values):
        sequenced_array = np.asarray(is_sequenced_values, dtype=bool)
        if sequenced_array.size == 0:
            return []
        return np.where(
            sequenced_array,
            self.value,
            self.not_sequenced_val,
        ).tolist()

    def get_values(self, is_sequenced_values:list=None):
        if is_sequenced_values is not None and len(is_sequenced_values) > 0:
            return self._build_values(is_sequenced_values)
        
        elif self.is_sequenced is not None and len(self.is_sequenced) > 0:
            return self._build_values(self.is_sequenced)
        
        else:
            return self.value
    
    def update(self, new_val=None, new_not_sequenced_val=None):
        if new_val is not None:
            self.value = new_val
        if new_not_sequenced_val is not None:
            self.not_sequenced_val = new_not_sequenced_val
            
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
        series = series.copy().astype(object)
        not_sequenced = series.isna()
        print(f"[handle_nans] NaN count: {not_sequenced.sum()} / {len(series)}")
        series.loc[not_sequenced] = self.not_sequenced_val
        return series, ~not_sequenced
    
    def apply_filter(self, series, sequenced):
        pass

    def get_values(self, series):
        print(f"[get_values] series unique: {series.dropna().unique()[:5]}")
        result = self.apply_filter(*self.handle_nans(series))
        print(
            "[get_values] result counts:",
            pd.Series(result).value_counts()
        )
        print("GET_VALUES")
        print("self id =", id(self))
        print("CLASS =", type(self).__name__)
        print("SELF ID =", id(self))

        if hasattr(self, "category_colors"):
            print("category_colors =", self.category_colors)
        print("series unique =", pd.Series(series).dropna().unique())
        return result
    

class Scaler(DynamicParameter):

    def __init__(self, minv, maxv, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.minv = minv
        self.maxv = maxv
        self.avoid_scaling = False

    def scale_values(self, values):
        if self.avoid_scaling:
            return [1.0] * len(values)
        try:
            values_array = pd.to_numeric(values, errors="raise").to_numpy(dtype=float, copy=False)
        except (ValueError, TypeError):
            values_array = np.arange(len(values), dtype=float)
        range_val = self.max_counts - self.min_counts if self.max_counts != self.min_counts else 1
        scaled_values = self.minv + (values_array - self.min_counts) * (self.maxv - self.minv) / range_val
        return scaled_values.tolist()

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
    
    def update(self, series, new_minv=None, new_maxv=None, avoid_scaling=False):
        self.avoid_scaling = avoid_scaling
        if avoid_scaling:
            self.minv = self.maxv
        else:
            if new_minv is not None:
                self.minv = new_minv
                if new_minv > self.maxv:
                    self.maxv = new_minv
            if new_maxv is not None:
                self.maxv = new_maxv
                if new_maxv < self.minv:
                    self.minv = new_maxv
        return self.get_values(series)


class CategoricalColor(DynamicParameter):

    PALETTE = [
        "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
        "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
        "#dcbeff", "#9A6324", "#fffac8", "#800000", "#aaffc3",
        "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#ffffff"
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.category_colors = {}

    def ensure_colors(self, categories):
        missing_categories = pd.Index(categories).difference(self.category_colors.keys())
        for category in missing_categories:
            idx = len(self.category_colors) % len(self.PALETTE)
            self.category_colors[category] = self.PALETTE[idx]

    def apply_filter(self, series, sequenced):
        print(f"[CategoricalColor] series unique: {series.unique()[:5]}")
        print(f"[CategoricalColor] sequenced sum: {sequenced.sum()} / {len(sequenced)}")
        result_series = series.copy()
        sequenced_series = series[sequenced].astype(str)
        unique_categories = pd.unique(sequenced_series)
        self.ensure_colors(unique_categories)
        print(f"[CategoricalColor] category_colors: {self.category_colors}")
        result_series.loc[sequenced] = sequenced_series.map(self.category_colors)
        print(f"[CategoricalColor] result unique: {pd.unique(result_series)[:5]}")
        return result_series.tolist()
