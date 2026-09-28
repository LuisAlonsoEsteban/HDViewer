
import geopandas as gpd
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse

class Pooler:

    join_result = None
    filtered_adata = None
    pooled_adata = None

    def __init__(self, 
                 bins: gpd.GeoDataFrame, 
                 adata: sc.AnnData, 
                 gdf: gpd.GeoDataFrame,
                 bins_index_col: str='index',
                 gdf_index_col: str='index',
                 keep_cols: list=["array_row", "array_col"]
                 ):
        
        self.bins = bins.copy()
        self.bins_index_col = bins_index_col
        assert self.bins_index_col in self.bins.columns, f"The placeholder column for indexes '{self.bins_index_col}' is not present in the bins GeoDataFrame"
        self.bins.rename(columns={self.bins_index_col: 'bins_index'}, inplace=True)

        self.adata = adata
        # assert self.bins[self.index_col].isin(self.adata.obs.index)

        self.gdf = gdf.copy()
        self.gdf_index_col= gdf_index_col
        assert self.gdf_index_col in self.gdf.columns, f"The placeholder column for indexes '{self.gdf_index_col}' is not present in the gdf"
        self.gdf.rename(columns={self.gdf_index_col: 'gdf_index'}, inplace=True)

        self.keep_cols = keep_cols
        assert all([col in self.gdf.columns for col in self.keep_cols]), "The columns in 'keep_cols' must be in the GeoDataFrame"


    def join(
            self, 
            predicate='within'
            ):
        
        if predicate not in self.bins.sindex.valid_query_predicates.intersection(self.gdf.sindex.valid_query_predicates):
            raise ValueError(f"Predicate '{predicate}' is not valid for the given geometries.")
        
        self.join_result = gpd.sjoin(self.bins, self.gdf, how='right', predicate=predicate)
        if (self.join_result["bins_index"].notna()).sum() == 0:
            raise ValueError("No geometries satisfy the predicate")

        self.join_result.set_index('bins_index', drop=False, inplace=True)
        self.join_result.index.name = None
        
        if len(self.keep_cols) > 0:
            for col in self.keep_cols:
                if col in self.bins.columns and col in self.gdf.columns:
                    self.join_result[col] = self.join_result[f'{col}_right']
                    self.join_result.drop(columns=[f'{col}_right', f'{col}_left'], inplace=True)

        return self.join_result
    
    def filter(self):
        filtered_obs_mask = self.adata.obs_names.isin(self.join_result['bins_index'])
        self.filtered_adata = self.adata[filtered_obs_mask, :]
        self.filtered_adata.obs.drop(columns=[col for col in self.keep_cols if col in self.filtered_adata.obs.columns], inplace=True)
        self.filtered_adata.obs = pd.merge(self.filtered_adata.obs, self.join_result[self.keep_cols + ['gdf_index']], left_index=True, right_index=True)

        return self.filtered_adata
    
    def pool(self):
        groupby_object = self.filtered_adata.obs.groupby(['gdf_index'] + self.keep_cols, observed=True)
        counts = self.filtered_adata.X if 'counts' not in self.filtered_adata.layers.keys() else self.filtered_adata.layers['counts']

        # Initialize a sparse matrix to store the summed gene counts for each spot
        summed_counts = sparse.lil_matrix((groupby_object.ngroups, counts.shape[1]))
        polygon_id = []
        metadata = {col: [] for col in self.keep_cols}
        
        if len(self.keep_cols) > 0:
            for row, ((geom, *cols), idx_) in enumerate(groupby_object.indices.items()):
                summed_counts[row] = counts[idx_].sum(0)
                polygon_id.append(geom)
                for i,col in enumerate(metadata.keys()):
                    metadata[col].append(cols[i])
        else:
            for row, (geom, idx_) in enumerate(groupby_object.indices.items()):
                summed_counts[row] = counts[idx_].sum(0)
                polygon_id.append(geom)
        
        summed_counts = summed_counts.tocsr()
        self.pooled_adata = sc.AnnData(X=summed_counts, obs=pd.DataFrame(metadata, index=polygon_id), var=self.filtered_adata.var)
        self.pooled_adata.var_names_make_unique()

        if 'array_row' in self.pooled_adata.obs.columns and 'array_col' in self.pooled_adata.obs.columns:
            self.pooled_adata.obsm = {'spatial': self.pooled_adata.obs[['array_row', 'array_col']].values}

        return self.pooled_adata
    
    def run(self, predicate='within'):
        self.join(predicate=predicate)
        self.filter()
        return self.pool()
    

class CellsPooler(Pooler):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def filter(self): 
        self.join_result['is_within_polygon'] = ~self.join_result['bins_index'].isna()
        isolated_bins = (~self.join_result['is_within_polygon']).sum()
        total_bins = self.bins.shape[0]
        print(f"Detected {isolated_bins} bins (out of {total_bins}, {(100 * isolated_bins / total_bins):.2f}%) that pointed to no cell")

        barcodes_in_overlaping_polygons = pd.unique(self.join_result[self.join_result.duplicated(subset=['bins_index'])]['bins_index'])
        overlapping_bins = barcodes_in_overlaping_polygons.shape[0]
        new_total_bins = total_bins - isolated_bins
        self.join_result['is_not_in_a_polygon_overlap'] = ~self.join_result['bins_index'].isin(barcodes_in_overlaping_polygons)
        print(f"Detected {overlapping_bins} bins (out of {new_total_bins}, {(100 * overlapping_bins / new_total_bins):.2f}%) that pointed to a cell jointly with another bin/s")
           
        self.barcodes_in_one_polygon = self.join_result[self.join_result['is_within_polygon'] & self.join_result['is_not_in_a_polygon_overlap']]

        filtered_obs_mask = self.adata.obs_names.isin(self.barcodes_in_one_polygon['bins_index'])
        self.filtered_adata = self.adata[filtered_obs_mask, :]
        self.filtered_adata.obs.drop(columns=[col for col in self.keep_cols if col in self.filtered_adata.obs.columns], inplace=True)
        self.filtered_adata.obs = pd.merge(self.filtered_adata.obs, self.join_result[self.keep_cols + ['gdf_index']], left_index=True, right_index=True)

        return self.filtered_adata

