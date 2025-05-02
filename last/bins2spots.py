
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point, Polygon

import scipy.sparse as sparse

import scanpy as sc

from tqdm import tqdm

from .converter import Converter


class Bins2Spots:

    def __init__(self, adata, tissue_positions, scale:str="px", pixels_per_bin:float=None):
        self.adata_hd = adata
        self.df_tissue_positions = tissue_positions

        self.scale = scale
        self.converter = Converter(pixels_per_bin=pixels_per_bin, df_tissue_positions=tissue_positions)

        self._bin_points = None
        self._spot_polygons = None
        self.adata_v2 = None


    @property
    def bin_points(self):
        if self._bin_points is None:
            self.create_bin_points()
        return self._bin_points
    
    @property
    def spot_polygons(self):
        if self._spot_polygons is None:
            self.create_spot_polygons()
        return self._spot_polygons
    
    def get_rotation(self, rad=False):
        return self.converter.visiumHD.get_rotation(self.df_tissue_positions, rad=rad)

    
    def create_bin_points(self, subset=False):
        geometry = [Point(y,x) for (x,y) in tqdm(zip(self.df_tissue_positions['pxl_row_in_fullres' if self.scale=="px" else 'array_row'], self.df_tissue_positions['pxl_col_in_fullres' if self.scale=="px" else 'array_row']))]
        self._bin_points = gpd.GeoDataFrame(self.df_tissue_positions, geometry=geometry)
        return self._bin_points[self._bin_points["in_tissue"]] if subset else self._bin_points

    def create_spot_polygons(self):
        spots = []
        metadata = {'id': [], 'array_row': [], 'array_col': [], 'pxl_row_in_fullres': [], 'pxl_col_in_fullres': []}
        
        nbins_rows, nbins_cols = pd.unique(self.bin_points["array_row"]).shape[0], pd.unique(self.bin_points["array_col"]).shape[0] 
        nspots_rows = int(self.converter.bins2spots(nbins_rows, axis=0))
        nspots_cols = int(self.converter.bins2spots(nbins_cols, axis=1))

        id = 0
        scaler_func = self.converter.spots2pixels if self.scale=="px" else self.converter.spots2bins
        scaled_radius = self.converter.get_pixels(self.converter.visiumV2.SPOT_DIAMETER / 2) if self.scale=="px" else self.converter.get_bins(self.converter.visiumV2.SPOT_DIAMETER / 2)
        max_scaled_x = scaler_func(nspots_rows)
        for i, x in enumerate(tqdm(range(nspots_rows))):
            for y in range(i % 2, nspots_cols, 2):
                metadata['id'].append(id)
                id += 1

                metadata['array_row'].append(x)
                metadata['array_col'].append(y)

                scaled_x = scaler_func(x, axis=0, shift=self.converter.visiumV2.SPOT_DIAMETER / 2)
                scaled_y = scaler_func(y, axis=1, shift=self.converter.visiumV2.SPOT_DIAMETER / 2)
                metadata['pxl_row_in_fullres'].append(scaled_x)
                # metadata['pxl_row_in_fullres'].append(max_scaled_x - scaled_x) # # Note that the direction of 'array_row' and 'pxl_row_in_fullres' are switched: as one increases, the other decreases and viceversa
                metadata['pxl_col_in_fullres'].append(scaled_y)

                spots.append(Point(scaled_y, scaled_x).buffer(scaled_radius))

        self._spot_polygons = gpd.GeoDataFrame(metadata, geometry=spots)

        return self._spot_polygons

    def join(self):
        points = self.bin_points
        polygons = self.spot_polygons

        result_spatial_join = gpd.sjoin(polygons, points, how='left', predicate='intersects')
        result_spatial_join.index = result_spatial_join['index']
        result_spatial_join[["array_row", "array_col"]] = result_spatial_join[["array_row_left", "array_col_left"]]
        result_spatial_join[["pxl_row_in_fullres", "pxl_col_in_fullres"]] = result_spatial_join[["pxl_row_in_fullres_left", "pxl_col_in_fullres_left"]]

        filtered_obs_mask = self.adata_hd.obs_names.isin(result_spatial_join['index'])
        filtered_adata = self.adata_hd[filtered_obs_mask, :]
        filtered_adata.obs.drop(['array_row', 'array_col', 'pxl_row_in_fullres','pxl_col_in_fullres'], axis=1, inplace=True)

        # Add the results of the point spatial join to the Anndata object
        filtered_adata.obs = pd.merge(filtered_adata.obs, result_spatial_join[['id', 'array_row', 'array_col', 'pxl_row_in_fullres','pxl_col_in_fullres', 'geometry']], left_index=True, right_index=True)

        return filtered_adata

    def pool_spots(self, filtered_adata):
        # Group the data by unique spot IDs
        groupby_object = filtered_adata.obs.groupby(['id', 'array_row', 'array_col', 'pxl_row_in_fullres','pxl_col_in_fullres'], observed=True)
        counts = filtered_adata.X if 'counts' not in filtered_adata.layers.keys() else filtered_adata.layers['counts']

        # Initialize a sparse matrix to store the summed gene counts for each spot
        summed_counts = sparse.lil_matrix((groupby_object.ngroups, counts.shape[1]))
        polygon_id = []
        coords = {'array_row': [], 'array_col': [], 'pxl_row_in_fullres': [], 'pxl_col_in_fullres': []}

        for row, ((polygons, arr_row, arr_col, pxl_row, pxl_col), idx_) in enumerate(groupby_object.indices.items()):
            summed_counts[row] = counts[idx_].sum(0)
            polygon_id.append(polygons)
            coords['array_row'].append(arr_row)
            coords['array_col'].append(arr_col)
            coords['pxl_row_in_fullres'].append(pxl_row)
            coords['pxl_col_in_fullres'].append(pxl_col)

        summed_counts = summed_counts.tocsr()
        self.adata_v2 = sc.AnnData(X=summed_counts, obs=pd.DataFrame(coords, index=polygon_id), var=filtered_adata.var)
        self.adata_v2.obsm = {'spatial': self.adata_v2.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values}

        return self.adata_v2
    
    def convert(self, normalize=True):
        self.pool_spots(self.join())
        if normalize:
            sc.pp.normalize_total(self.adata_hd, target_sum=1e4)
            sc.pp.normalize_total(self.adata_v2, target_sum=1e4)
        return self.adata_v2 
