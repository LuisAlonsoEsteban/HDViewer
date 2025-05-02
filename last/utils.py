
from shapely.geometry import Polygon

def create_square(point, size=1):
    """Create a square polygon around the point with the given size.
    
    Args:
    point (shapely.geometry.Point): The central point of the square.
    size (float): The length of each side of the square.
    
    Returns:
    shapely.geometry.Polygon: The square polygon.
    """
    x, y = point.x, point.y
    half_size = size / 2
    square = Polygon([
        (x - half_size, y - half_size),
        (x - half_size, y + half_size),
        (x + half_size, y + half_size),
        (x + half_size, y - half_size)
    ])
    return square

import numpy as np

def compress_img(img):
    if img.dtype == 'uint8' and img.ndims == 3:
        return (img[:, :, 0].astype(np.uint32) << 16) | (img[:, :, 1].astype(np.uint32) << 8)  | (img[:, :, 2].astype(np.uint32))

import geopandas as gpd
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse

def join_and_pool(adata, geometries, new_geometries, keep_cols=["array_row", "array_col"], predicate='within'):
    if predicate not in geometries.sindex.valid_query_predicates.intersection(new_geometries.sindex.valid_query_predicates):
        raise ValueError(f"Predicate '{predicate}' is not valid for the given geometries.")
    
    join_result = gpd.sjoin(new_geometries, geometries, how='left', predicate=predicate)
    join_result.index = join_result['index']
    join_result[keep_cols] = join_result[[f'{col}_left' for col in keep_cols]]

    filtered_obs_mask = adata.obs_names.isin(join_result['index'])
    filtered_adata = adata[filtered_obs_mask, :]
    filtered_adata.obs.drop(keep_cols, axis=1, inplace=True)

    filtered_adata.obs = pd.merge(filtered_adata.obs, join_result[['id', 'geometry'] + keep_cols], left_index=True, right_index=True)

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
    join_adata = sc.AnnData(X=summed_counts, obs=pd.DataFrame(metadata, index=polygon_id), var=filtered_adata.var)
    join_adata.obsm = {'spatial': join_adata.obs[['array_row', 'array_col']].values}

    return join_adata
