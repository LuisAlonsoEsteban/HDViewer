
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

        self.adata_processed = None

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
            col = self.adata.obs[colname]
            if hasattr(col, 'cat'):
                col = col.cat.as_ordered()
            return col.min(), col.max()
    
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
         
    def get_value(self, adata=None, gdf=None, colname=None):
        if gdf is None:
            gdf = self.gdf

        if adata is None:
            adata = self.adata

        if colname in gdf.columns:
            return gdf[colname].values.tolist()

        elif self.has_adata():            
            value_df = self.get_adata_vector(adata, colname)
            value_df[self.on] = value_df.index.astype(gdf[self.on].dtype)
            
            merged_gdf = gdf.reset_index().merge(value_df, on=self.on, how='left').set_index('index')
            print(f"Colname: {colname}, NaNs after merge: {merged_gdf[colname].isna().sum()}/{len(merged_gdf)}")
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
    
    def get_marker_genes(self, n_genes=5):
        import scanpy as sc
        import pandas as pd
        if not hasattr(self, 'adata_processed') or 'rank_genes_groups' not in self.adata_processed.uns:
            return {}
        
        result = self.adata_processed.uns['rank_genes_groups']
        clusters = result['names'].dtype.names
        markers = {}
        for cluster in clusters:
            genes = [result['names'][cluster][i] for i in range(n_genes)]
            scores = [result['scores'][cluster][i] for i in range(n_genes)]
            markers[cluster] = list(zip(genes, [round(s, 2) for s in scores]))
        return markers
    
    def get_colnames(self):
        return self.adata.var_names.tolist() + self.adata.obs.columns.tolist()
    
    def is_continous(self, colname):
        if colname in self.get_colnames():
            if colname in self.adata.var_names:
                return True
            else:
                return pd.api.types.is_numeric_dtype(self.adata.obs[colname])
        else:
            raise ValueError(f"{colname} not found in 'adata'.")

    def run_clustering(self, resolution=0.9, n_pcs=30, n_neighbors=15):
        import scanpy as sc

        if not self.adata_processed:
            adata = self.adata.copy()
            adata.var_names_make_unique()

            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata)
            adata = adata[:, adata.var.highly_variable].copy()
            sc.pp.scale(adata, max_value=10)
            sc.pp.pca(adata, n_comps=n_pcs)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

            self.adata_processed = adata

        sc.tl.leiden(self.adata_processed, resolution=resolution)
#
        common_idx = self.adata.obs_names.intersection(self.adata_processed.obs_names)

        self.adata.obs["leiden"] = None
        self.adata.obs.loc[common_idx, "leiden"] = (
            self.adata_processed.obs.loc[common_idx, "leiden"]
            .astype(str)
            .values
        )

        sc.tl.rank_genes_groups(self.adata_processed, groupby='leiden', method='wilcoxon')
        

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

