
import geopandas as gpd
import pandas as pd

class Sample:

    _max_y = None

    def __init__(
            self, 
            gdf, 
            adata=None, 
            pixels_per_bin=None, 
            on='index'):
        
        self.gdf = gdf
        self.adata = adata

        self.pixels_per_bin = pixels_per_bin

        self.on = on
        if self.adata is not None:
            assert self.on in self.gdf.columns, "On parameter is incorrect. It must be a column of 'gdf'"

    @property
    def n_geometries(self):
        return self.gdf.shape[0]
    
    @property
    def n_not_sequenced(self):
        colname = self.adata.var_names.tolist()[0]
        return int(self.get_value(self.adata, self.gdf, colname=colname).isna().sum())
                
    def has_adata(self):
        return self.adata is not None
    
    def get_range(self, colname):
        if colname in self.adata.var_names:
            return self.adata[:, colname].X.min(), self.adata[:, colname].X.max()
        elif colname in self.adata.obs.columns:
            return self.adata.obs[colname].min(), self.adata.obs[colname].max()
    
    def crop_geom(self, left, right, top, bottom):
        return self.gdf.cx[left:right, bottom:top].copy()

    def crop(self, left, right, top, bottom, crop_adata=False):
        cropped_gdf =  self.crop_geom(left, right, top, bottom)
        if crop_adata:
            cropped_mask = self.adata.obs_names.isin(cropped_gdf.index)
            cropped_adata = self.adata[cropped_mask, :]
            return cropped_gdf, cropped_adata
        else:
            return cropped_gdf
    
    def get_x(self, geom):
        pass

    def get_y(self, geom):
        pass
    
    def get_coordinates(self, gdf) -> pd.DataFrame:
        gdf['x'] = gdf.geometry.apply(lambda geom: self.get_x(geom))
        gdf['y'] = gdf.geometry.apply(lambda geom: self.get_y(geom))
        return gdf[['x', 'y']]

    def get_adata_vector(self, adata, colname):
        if colname in adata.var_names:
            return adata[:, colname].to_df()
        elif colname in adata.obs.columns:
            return adata.obs[colname].to_frame()
        else:
            raise ValueError(f"'{colname}' not found in adata")
         
    def get_value(self, adata, gdf, colname=None):
        if colname in gdf.columns:
            return gdf[colname].values.tolist()

        elif self.has_adata():
            value_df = self.get_adata_vector(adata, colname)
            value_df[self.on] = value_df.index.astype(gdf[self.on].dtype)
            
            merged_gdf = gdf.reset_index().merge(value_df, on=self.on, how='left').set_index('index')

            return merged_gdf[colname]
        else:
            raise ValueError(f"'{colname}' not found in geometries, and AnnData not available")

    
    def fetch_data(self, bbox=None, colnames=None, return_coords=True) -> pd.DataFrame:
        if bbox is not None:
            if self.has_adata():
                cropped_gdf, cropped_adata = self.crop(*bbox, crop_adata=True)
            else:
                cropped_gdf = self.crop(*bbox, crop_adata=False)
        else:
            cropped_gdf, cropped_adata = self.gdf, self.adata

        if return_coords:
            data = self.get_coordinates(cropped_gdf)
        else:
            data = pd.DataFrame()

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
            for colname in colnames:
                data[colname] = self.get_value(cropped_adata, cropped_gdf, colname)

        return data
    
    def get_colnames(self):
        return self.adata.var_names.tolist() + self.adata.obs.columns.tolist()

        
class Bins(Sample):
    def __init__(self, *args, half_bin_size=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.half_bin_size = half_bin_size or 0

    def crop_geom(self, left, right, top, bottom):
        if left is None: left = self.half_bin_size
        if right is None: right = self.half_bin_size
        if bottom is None: bottom = self.half_bin_size
        if top is None: top = self.half_bin_size

        return self.gdf.cx[left - self.half_bin_size : right + self.half_bin_size,
                           bottom - self.half_bin_size : top + self.half_bin_size].copy()
            
    def get_x(self, geom):
        return geom.centroid.x
    
    def get_y(self, geom):
        return geom.centroid.y
    
class Spots(Sample):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_x(self, geom):
        return geom.centroid.x
    
    def get_y(self, geom):
        return geom.centroid.y

class Cells(Sample):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_x(self, geom):
        return list(zip(*geom.exterior.coords))[0]
    
    def get_y(self, geom):
        return list(zip(*geom.exterior.coords))[1]

