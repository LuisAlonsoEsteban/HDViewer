
import geopandas as gpd
import numpy as np
import pandas as pd
from collections import OrderedDict
import math

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
         
    #def get_value(self, adata=None, gdf=None, colname=None):
    #    if gdf is None:
    #        gdf = self.gdf
#
    #    if adata is None:
    #        adata = self.adata
#
    #    if colname in gdf.columns:
    #        return gdf[colname].values.tolist()
#
    #    elif self.has_adata():            
    #        value_df = self.get_adata_vector(adata, colname)
    #        value_df[self.on] = value_df.index.astype(gdf[self.on].dtype)
    #        
    #        # ← AÑADE ESTO TEMPORALMENTE
    #        print(f"[get_value] colname={colname}")
    #        print(f"  value_df index dtype: {value_df.index.dtype}")
    #        print(f"  gdf[self.on] dtype: {gdf[self.on].dtype}")
    #        print(f"  value_df[self.on] sample: {value_df[self.on].head(3).tolist()}")
    #        print(f"  gdf[self.on] sample: {gdf[self.on].head(3).tolist()}")
    #        print(f"  value_df unique non-null: {value_df[colname].dropna().unique()[:5]}")
    #        
    #        merged_gdf = gdf.reset_index().merge(value_df, on=self.on, how='left').set_index('index')
    #        print(f"  merged result unique: {merged_gdf[colname].dropna().unique()[:5]}")
    #        return merged_gdf[colname]
    #    else:
    #        raise ValueError(f"'{colname}' not found in geometries, and AnnData not available")

    def get_value(self, adata=None, gdf=None, colname=None):
        if gdf is None:
            gdf = self.gdf
        if adata is None:
            adata = self.adata

        if colname in gdf.columns:
            return gdf[colname].values.tolist()

        elif self.has_adata():
            value_df = self.get_adata_vector(adata, colname)
            
            # Columnas de obs (leiden, Group 1, etc.) → map directo por obs_names
            if colname in adata.obs.columns:
                obs_series = self.get_obs_names(adata, gdf)
                result = obs_series.map(value_df[colname])
                result.index = gdf.index
                return result
            
            # Genes (adata.var_names) → merge original por self.on
            value_df[self.on] = value_df.index.astype(gdf[self.on].dtype)
            merged_gdf = gdf.reset_index().merge(value_df, on=self.on, how='left').set_index('index')
            return merged_gdf[colname]
        
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
        
        cluster_labels = self.adata_processed.obs['leiden']

        for cluster in clusters:
            genes = [result['names'][cluster][i] for i in range(n_genes)]
            scores = [result['scores'][cluster][i] for i in range(n_genes)]
            pvals = [result['pvals_adj'][cluster][i] for i in range(n_genes)]
            
            # Máscara booleana para separar las células que están dentro y fuera del cluster actual
            cluster_mask = cluster_labels == cluster
            logfcs = []
            
            for gene in genes:
                if gene in self.adata_processed.var_names:
                    # Extraemos la expresión del gen de forma segura (sea matriz densa o dispersa)
                    expr = self.adata_processed[:, gene].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    
                    # Calculamos las medias logarítmicas
                    mean_in = expr[cluster_mask].mean()
                    mean_out = expr[~cluster_mask].mean()
                    
                    # Como los datos de Scanpy ya están en log1p, aplicamos la resta directa
                    if mean_in == 0 and mean_out == 0:
                        l2fc = 0.0
                    elif mean_out == 0:
                        # Si fuera del cluster es 0, usamos pseudocuenta de seguridad para evitar división por cero
                        l2fc = np.log2((mean_in + 1e-6) / 1e-6)
                    else:
                        l2fc = mean_in - mean_out
                        
                    logfcs.append(round(float(l2fc), 2))
                else:
                    logfcs.append(0.0)
            

            # Empaquetamos los 4 elementos asegurando que todos mantengan tipos numéricos puros (float)
            markers[cluster] = list(zip(
                genes, 
                [round(float(s), 2) for s in scores], 
                [float(p) for p in pvals], 
                [float(l) for l in logfcs]
            ))
            
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
    AGGREGATED_PREFIX = "agg:"
    AGGREGATION_FACTORS = (1, 2, 4, 8, 12)
    MAX_CACHE_ENTRIES = 24

    def __init__(self, *args, half_bin_size=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.half_bin_size = half_bin_size or 0
        self._multiscale_cache = OrderedDict()

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

    @property
    def native_width(self):
        return float(self.pixels_per_bin or 0)

    @property
    def native_height(self):
        return float(self.pixels_per_bin or 0)

    def get_aggregation_factor(self, resolution_level):
        if resolution_level is None:
            return 1
        if isinstance(resolution_level, str) and resolution_level.startswith(self.AGGREGATED_PREFIX):
            return max(1, int(resolution_level.split(":", 1)[1]))
        return 1

    def get_resolution_level_for_factor(self, factor):
        if factor <= 1:
            return None
        return f"{self.AGGREGATED_PREFIX}{int(factor)}"

    def _cache_result(self, cache_key, value):
        self._multiscale_cache[cache_key] = value
        self._multiscale_cache.move_to_end(cache_key)
        while len(self._multiscale_cache) > self.MAX_CACHE_ENTRIES:
            self._multiscale_cache.popitem(last=False)

    def _align_bbox(self, bbox, cell_width, cell_height):
        if bbox is None:
            min_x, min_y, max_x, max_y = self.gdf.total_bounds
            bbox = [min_x, max_x, min_y, max_y]

        left, right, bottom, top = bbox

        aligned_left = np.floor(left / cell_width) * cell_width
        aligned_right = np.ceil(right / cell_width) * cell_width
        aligned_bottom = np.floor(bottom / cell_height) * cell_height
        aligned_top = np.ceil(top / cell_height) * cell_height

        return [aligned_left, aligned_right, aligned_bottom, aligned_top]

    def _get_numeric_value(self, adata, gdf, colname):
        values = self.get_value(adata=adata, gdf=gdf, colname=colname)
        if not isinstance(values, pd.Series):
            values = pd.Series(values, index=gdf.index)
        return pd.to_numeric(values, errors="coerce")

    def _fetch_aggregated_data(self, bbox=None, colnames=None, factor=1):
        cell_width = self.native_width * factor
        cell_height = self.native_height * factor

        aligned_bbox = self._align_bbox(bbox, cell_width, cell_height)

        cache_key = (
            factor,
            tuple(round(value, 4) for value in aligned_bbox),
            tuple(sorted(colnames or [])),
        )
        
        if cache_key in self._multiscale_cache:
            self._multiscale_cache.move_to_end(cache_key)
            return self._multiscale_cache[cache_key].copy()

        cropped_gdf, cropped_adata = self.crop(*aligned_bbox, crop_adata=True)
        if cropped_gdf.empty:
            empty = pd.DataFrame(columns=["x", "y", "width", "height", "obs_name"])
            self._cache_result(cache_key, empty)
            return empty.copy()

        coords = self.get_coordinates(cropped_gdf.copy())
        grid_x = np.floor(coords["x"].to_numpy() / cell_width).astype(np.int32)
        grid_y = np.floor(coords["y"].to_numpy() / cell_height).astype(np.int32)

        data = pd.DataFrame(
            {
                "grid_x": grid_x,
                "grid_y": grid_y,
                "x": (grid_x + 0.5) * cell_width,
                "y": (grid_y + 0.5) * cell_height,
            },
            index=coords.index,
        )
        data["width"] = cell_width
        data["height"] = cell_height
        data["obs_name"] = (
            "summary:"
            + data["grid_x"].astype(str)
            + ":"
            + data["grid_y"].astype(str)
            + f":{factor}"
        )

        grouped = data.groupby(["grid_x", "grid_y"], sort=False).agg(
            x=("x", "first"),
            y=("y", "first"),
            width=("width", "first"),
            height=("height", "first"),
            obs_name=("obs_name", "first"),
        )

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
            for colname in colnames:
                values = self._get_numeric_value(cropped_adata, cropped_gdf, colname)
                value_frame = pd.DataFrame(
                    {
                        "grid_x": grid_x,
                        "grid_y": grid_y,
                        colname: values.values,
                    },
                    index=cropped_gdf.index,
                )
                grouped[colname] = value_frame.groupby(["grid_x", "grid_y"], sort=False)[colname].sum(min_count=1)

        grouped.reset_index(drop=True, inplace=True)
        grouped.index = grouped["obs_name"]
        self._cache_result(cache_key, grouped)
        return grouped.copy()

    def fetch_data(self, bbox=None, colnames=None, return_coords=True, resolution_level=None) -> pd.DataFrame:
        factor = self.get_aggregation_factor(resolution_level)
        if factor <= 1:
            data = super().fetch_data(bbox=bbox, colnames=colnames, return_coords=return_coords)
            if return_coords:
                data["width"] = self.native_width
                data["height"] = self.native_height
            return data

        return self._fetch_aggregated_data(bbox=bbox, colnames=colnames, factor=factor)
    
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
