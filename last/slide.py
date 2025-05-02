
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import numpy as np
import pandas as pd

from .image import Image
from .converter import Converter

from . import config

from shapely.affinity import affine_transform
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show

class PlotConfig:
    def __init__(self, **kwargs):
        self.cmap = kwargs.get('cmap', plt.cm.get_cmap('viridis'))

        self.color = kwargs.get('color', 'grey')
        self.color_by = kwargs.get('color_by', None)
        self.edgecolor = kwargs.get('edgecolor', 'black')
        self.edgecolor_by = kwargs.get('edgecolor_by', None)
        self.linewidth = kwargs.get('linewidth', 0.5)
        self.linewidth_by = kwargs.get('linewidth_by', None)
        self.alpha = kwargs.get('alpha', 0.5)
        self.alpha_by = kwargs.get('alpha_by', None)

    def set_color(self, merged_gdf, attr='color'):
        attr_by = getattr(self, f"{attr}_by")
        if attr_by in merged_gdf.columns:
            norm = mcolors.Normalize(vmin=merged_gdf[attr_by].min(), vmax=merged_gdf[attr_by].max())
            setattr(self, attr, self.cmap(norm(merged_gdf[attr_by])).tolist())

    def scale_values(self, merged_gdf, attr):
        return config[attr]['min'] + (merged_gdf[getattr(self, f"{attr}_by")] - merged_gdf[getattr(self, f"{attr}_by")].min()) * (config[attr]['max'] - config[attr]['min']) / (merged_gdf[getattr(self, f"{attr}_by")].max() - merged_gdf[getattr(self, f"{attr}_by")].min())

    def alpha_in_color(self):
        return (isinstance(self.color[0], list) and len(self.color[0]) == 4) or (isinstance(self.color[0], float) and len(self.color) == 4)

    def parse_attrs(self, merged_gdf):
        self.set_color(merged_gdf, 'color')
        self.set_color(merged_gdf, 'edgecolor')
        if self.linewidth_by in merged_gdf.columns:
            self.linewidth = self.scale_values(merged_gdf, 'linewidth')
        if not self.alpha_in_color(): 
            if self.alpha_by in merged_gdf.columns:
                self.alpha = self.scale_values(merged_gdf, 'alpha')
        else:
            self.alpha = None

    def get_nan_representation(self, attr):
        return config[attr]['nan']
    
    def set_nan_representation(self, attr, not_sequenced_indcs, n_polygons):
        if getattr(self, f"{attr}_by") is not None:
            if isinstance(getattr(self, attr), list):
                setattr(self, attr, [c if i not in not_sequenced_indcs else self.get_nan_representation(attr) for i,c in enumerate(getattr(self, attr))])
            elif isinstance(getattr(self, attr), pd.Series):
                getattr(self, attr).iloc[not_sequenced_indcs] = self.get_nan_representation(attr)
        else:
            setattr(self, attr, [getattr(self, attr) if i not in not_sequenced_indcs else self.get_nan_representation(attr) for i in range(n_polygons)])

    def set_not_sequenced_params(self, merged_gdf, view_not_sequenced=False):
        if view_not_sequenced and 'not_sequenced' in merged_gdf.columns:

            not_sequenced_indcs = np.where(merged_gdf['not_sequenced'])[0]
            n_polygons = len(merged_gdf)

            self.set_nan_representation('color', not_sequenced_indcs, n_polygons)
            self.set_nan_representation('edgecolor', not_sequenced_indcs, n_polygons)
            self.set_nan_representation('linewidth', not_sequenced_indcs, n_polygons)
            if not self.alpha_in_color():
                self.set_nan_representation('alpha', not_sequenced_indcs, n_polygons)

    def get_attributes(self):
        return self.color, self.edgecolor, self.linewidth, self.alpha
    
    def set_attributes(self, merged_gdf, view_not_sequenced=True):
        self.parse_attrs(merged_gdf)
        self.set_not_sequenced_params(merged_gdf, view_not_sequenced=view_not_sequenced)
        return self.get_attributes()


class Sample:

    merged_gdf = None

    def __init__(self, geometries, adata=None, transform=None, **kwargs):
        self.geometries = geometries
        self.adata = adata
        self.transform = transform
        self.plot_config = PlotConfig(**kwargs)
        self.factor = 1

    def has_data(self):
        return self.adata is not None

    def add_value(self, colname, on='id'):
        if colname in self.geometries.columns:
            self.merged_gdf = self.geometries
        elif self.has_data():
            if colname in self.adata.var_names:
                value_df = self.adata[:, colname].to_df()
            elif colname in self.adata.obs.columns:
                value_df = self.adata.obs[colname].to_frame()
            else:
                raise ValueError(f"'{colname}' not found in adata")

            value_df[on] = value_df.index.astype(self.geometries[on].dtype)
            
            self.merged_gdf = self.geometries.merge(value_df, on=on, how='left').copy()
            self.merged_gdf.loc[:,'not_sequenced'] = self.merged_gdf[colname].isna()

        else:
            raise ValueError(f"'{colname}' not found in geometries, and AnnData not available")

        return self.merged_gdf
    
    def add_single_transform(self, transform):
        if self.transform is None:
            self.transform = [transform]
        elif isinstance(self.transform, list):
            self.transform.append(transform)
        elif callable(self.transform):
            self.transform = [self.transform, transform]
        else:
            raise ValueError("Transform should be a callable or a list of callables.")
    
    def crop(self, gdf, bbox):
        if (gdf["geometry"].geom_type == "Point").any() and self.transform is not None and all(col in self.geometries.columns for col in ["pxl_row_in_fullres", "pxl_col_in_fullres"]): # Assume this is enough for knowing if a point is going to be expanded to bins
            half_bin_size = Converter(df_tissue_positions=self.geometries).visiumHD.pixels_per_bin / 2
            half_bin_size = half_bin_size if half_bin_size is not None else 0
            cropped_gdf = gdf.cx[bbox[0]-half_bin_size:bbox[1]+half_bin_size, bbox[2]-half_bin_size:bbox[3]+half_bin_size]
        else:
            cropped_gdf = gdf.cx[bbox[0]:bbox[1], bbox[2]:bbox[3]]
        return cropped_gdf
    
    def apply_transform(self, gdf):
        if self.transform is not None:
            if isinstance(self.transform, (list, tuple)):
                for t in self.transform:
                    gdf.loc[:,'geometry'] = gdf['geometry'].apply(t)
            else:
                gdf.loc[:,'geometry'] = gdf['geometry'].apply(self.transform)
        return gdf
    
    def plot(self, p, bbox=None, view_not_sequenced=False):
        if self.merged_gdf is not None:
            if bbox is not None:
                self.merged_gdf = self.crop(self.merged_gdf, bbox)
            self.merged_gdf = self.apply_transform(self.merged_gdf)
            color, edgecolor, linewidth, alpha = self.plot_config.set_attributes(self.merged_gdf, view_not_sequenced=view_not_sequenced)

            xs = self.merged_gdf.geometry.apply(lambda geom: np.array(geom.exterior.xy[0]))
            ys = self.merged_gdf.geometry.apply(lambda geom: np.array(geom.exterior.xy[1]))

            p.patches(xs, ys, fill_color=color, line_color=edgecolor, line_width=linewidth, fill_alpha=alpha)
            show(p)
        else:
            if bbox is not None:
                self.geometries = self.crop(self.geometries, bbox)
            self.geometries = self.apply_transform(self.geometries)
            color, edgecolor, linewidth, alpha = self.plot_config.set_attributes(self.geometries, view_not_sequenced=False)
            
            xs = self.geometries.geometry.apply(lambda geom: np.array(geom.exterior.xy[0]))
            ys = self.geometries.geometry.apply(lambda geom: np.array(geom.exterior.xy[1]))

            p.patches(xs, ys, fill_color=color, line_color=edgecolor, line_width=linewidth, fill_alpha=alpha)


class Slide:

    def __init__(self, img=None, samples=None):
        self.image = Image(img) if img is not None else None
        self.samples = samples
        if not self.samples:
            self.samples = []
        self.factor = 1
        

    def add_sample(self, sample):
        self.samples.append(sample)

    def clean_samples(self):
        self.samples = []
        self.factor = 1

    def set_img(self, img):
        self.image = Image(img)

    def get_image(self, bbox=None):
        downsampled_img, factor =  self.image.get_image_within_memory_limit(bbox) if self.image is not None else None
        if factor is not None: 
            self.factor = factor
            for sample in self.samples:
                sample.add_single_transform(lambda x: affine_transform(x, [1/factor, 0, 0, 1/factor, 0, 0]))
        return downsampled_img

    def plot(self, bbox:list=None, img=True, gene=None, view_not_sequenced=False, figsize=(10, 10)):
        p = figure(plot_width=figsize[0], plot_height=figsize[1], tools='pan,wheel_zoom,box_zoom,reset')
        
        if self.image is not None and img:
            downsampled_img = self.get_image(bbox)
            scaled_bbox = [int(e / self.factor) for e in bbox]
            p.image(image=[downsampled_img], x=scaled_bbox[0], y=scaled_bbox[2], dw=scaled_bbox[1]-scaled_bbox[0], dh=scaled_bbox[3]-scaled_bbox[2], alpha=0.7)
        
        for sample in self.samples:
            if gene is not None and sample.has_data():
                # TODO: Support subplots if several genes are passed
                sample.add_value(gene)
            sample.plot(p, bbox=bbox, view_not_sequenced=view_not_sequenced)

        if gene is not None:
            p.title.text = gene

        p.x_range.start = scaled_bbox[0]
        p.x_range.end = scaled_bbox[1]
        p.y_range.start = scaled_bbox[2]
        p.y_range.end = scaled_bbox[3]

        show(p)
