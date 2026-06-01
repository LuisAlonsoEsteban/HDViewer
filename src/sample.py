
import geopandas as gpd
import numpy as np
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
        self._embedding_keys = {
            "pca": "X_pca_2d",
            "umap": "X_umap",
            "tsne": "X_tsne",
        }

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
        
    def get_dummy_colname(self, dtype="numpy"):
        colname = self.get_value(colname=self.get_colnames()[0]).isna()
        if dtype == "numpy":
            return colname.to_numpy()
        elif dtype == "list":
            return colname.values.tolist()

    def get_obs_names(self, adata=None, gdf=None):
        if gdf is None:
            gdf = self.gdf

        if adata is None:
            adata = self.adata

        if adata is None:
            return gdf[self.on].astype(str)

        obs_df = pd.DataFrame({
            "obs_name": adata.obs_names.astype(str),
            self.on: adata.obs_names.astype(gdf[self.on].dtype),
        })
        merged_gdf = gdf.reset_index().merge(obs_df, on=self.on, how='left').set_index('index')
        # return merged_gdf["obs_name"].astype(str)
        return merged_gdf["obs_name"]

    
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

        data['obs_name'] = self.get_obs_names(cropped_adata, cropped_gdf).values

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
            for colname in colnames:
                data[colname] = self.get_value(cropped_adata, cropped_gdf, colname)

        return data
    
    def get_marker_genes(self, n_genes=5):
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

    def prepare_adata_processed(self, n_pcs=30, n_neighbors=15):
        import scanpy as sc

        if self.adata_processed is not None:
            return self.adata_processed

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
        return self.adata_processed

    def _store_embedding(self, method, values):
        key = self._embedding_keys[method]
        full = np.full((self.adata.n_obs, 2), np.nan, dtype=float)
        obs_index = pd.Index(self.adata.obs_names.astype(str))
        processed_index = pd.Index(self.adata_processed.obs_names.astype(str))
        common_idx = obs_index.get_indexer(processed_index)
        valid = common_idx >= 0
        full[common_idx[valid], :] = values[valid, :]
        self.adata.obsm[key] = full

    def compute_embedding(self, method: str, **params):
        import scanpy as sc

        method = method.lower()
        if method not in self._embedding_keys:
            raise ValueError(f"Unsupported embedding method: {method}")

        key = self._embedding_keys[method]
        if key in self.adata.obsm:
            return self.adata.obsm[key]

        adata = self.prepare_adata_processed(**params)

        if method == "pca":
            if "X_pca" not in adata.obsm:
                sc.pp.pca(adata, n_comps=params.get("n_pcs", 30))
            values = np.asarray(adata.obsm["X_pca"][:, :2], dtype=float)
        elif method == "umap":
            if "neighbors" not in adata.uns:
                sc.pp.neighbors(adata, n_neighbors=params.get("n_neighbors", 15), n_pcs=params.get("n_pcs", 30))
            if "X_umap" not in adata.obsm:
                sc.tl.umap(adata)
            values = np.asarray(adata.obsm["X_umap"][:, :2], dtype=float)
        else:
            if "X_tsne" not in adata.obsm:
                sc.tl.tsne(adata)
            values = np.asarray(adata.obsm["X_tsne"][:, :2], dtype=float)

        self._store_embedding(method, values)
        return self.adata.obsm[key]

    def get_embedding_dataframe(self, method: str) -> pd.DataFrame:
        method = method.lower()
        key = self._embedding_keys[method]
        if key not in self.adata.obsm:
            self.compute_embedding(method)

        values = np.asarray(self.adata.obsm[key], dtype=float)
        df = pd.DataFrame(values, columns=["x", "y"], index=self.adata.obs_names.astype(str))
        df.index.name = "obs_name"
        df.reset_index(inplace=True)
        df = df.dropna(subset=["x", "y"]).copy()
        return df

    def run_clustering(self, resolution=0.9, n_pcs=30, n_neighbors=15):
        import scanpy as sc

        self.prepare_adata_processed(n_pcs=n_pcs, n_neighbors=n_neighbors)

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
