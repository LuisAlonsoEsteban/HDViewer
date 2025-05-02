
import geopandas as gpd
import pandas as pd

import scanpy as sc
import scipy.sparse as sparse

class Comparison:

    sample1 = None
    sample2 = None

    join_result = None
    join_adata = None

    def __init__(self, sample1, sample2, bigger=1):
        self.sample1 = sample1
        self.sample2 = sample2
        self.bigger = bigger

        if not sample1.has_data() or not sample2.has_data():
            raise ValueError("Both samples must have both geometries and AnnData.")
        
    def join_and_pool(self, keep_cols=["array_row", "array_col"], predicate='within'):
        if predicate not in self.sample1.geometries.sindex.valid_query_predicates.intersection(self.sample2.geometries.sindex.valid_query_predicates):
            raise ValueError(f"Predicate '{predicate}' is not valid for the given geometries.")
        
        self.join_result = gpd.sjoin(self.sample1.geometries, self.sample2.geometries, how='left' if self.bigger == 1 else 'right', predicate=predicate)
        self.join_result.index = self.join_result['index']
        self.join_result[keep_cols] = self.join_result[[f'{col}_left' for col in keep_cols]]

        filtered_obs_mask = self.sample2.adata.obs_names.isin(self.join_result['index'])
        filtered_adata = self.sample2.adata[filtered_obs_mask, :]
        filtered_adata.obs.drop(keep_cols, axis=1, inplace=True)

        filtered_adata.obs = pd.merge(filtered_adata.obs, self.join_result[['id', 'geometry'] + keep_cols], left_index=True, right_index=True)

        groupby_object = filtered_adata.obs.groupby(['id'] + keep_cols, observed=True)
        counts = filtered_adata.X if 'counts' not in filtered_adata.layers.keys() else filtered_adata.layers['counts']

        # Initialize a sparse matrix to store the summed gene counts for each spot
        summed_counts = sparse.lil_matrix((groupby_object.ngroups, counts.shape[1]))
        polygon_id = []
        metadata = {col: [] for col in keep_cols}

        for row, ((geom, *cols), idx_) in enumerate(groupby_object.indices.items()):
            summed_counts[row] = counts[idx_].sum(0)
            polygon_id.append(geom)
            for i,col in enumerate(metadata.keys()):
                metadata[col].append(cols[i])

        summed_counts = summed_counts.tocsr()
        self.join_adata = sc.AnnData(X=summed_counts, obs=pd.DataFrame(metadata, index=polygon_id), var=self.filtered_adata.var)
        self.join_adata.obsm = {'spatial': self.join_adata.obs[['array_row', 'array_col']].values}

        return self.join_adata

    
