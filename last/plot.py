
from .bins2spots import Bins2Spots
from .slide import Slide, Sample

import scanpy as sc
import matplotlib.colors as mcolors
from .utils import create_square
from shapely import affinity
import math

class PlotHD:

    def __init__(self, img, adata, df_tissue_positions=None, pixels_per_bin:float=None):
        self.b2s = Bins2Spots(adata, df_tissue_positions, pixels_per_bin=pixels_per_bin)
        self.slide = Slide(img=img)

    @staticmethod
    def get_cmap(color1, color2, name='custom_diverging'):
        return mcolors.LinearSegmentedColormap.from_list(
            name, 
            [color1, color2]
        )
    
    def clean(self):
        self.slide.clean_samples()
        self.slide.image.clean_cache()

    def get_whole_bbox(self):
        """
        Get the bounding box of the whole image.
        """

        bbox = [
            0,
            self.slide.image.img.shape[1],
            0, 
            self.slide.image.img.shape[0]
        ]

        return bbox
    
    def get_max_side(self, bbox):
        """
        Get the maximum side of the bounding box.
        """

        bbox_width = bbox[1] - bbox[0]
        bbox_height = bbox[3] - bbox[2]

        return max(bbox_width, bbox_height)

        
    def plot(
            self, 
            img=True, 
            cells=None, 
            bbox:list=None, 
            gene:str=None, 
            view_spots:bool=False, 
            view_not_sequenced:bool=False, 
            normalize:bool=True
        ):

        """
        Plot the Visium HD data.
        """

        self.clean()

        if bbox is None:
            bbox = self.get_whole_bbox()
            
        bins_subset = self.b2s.bin_points.copy()
        bins_subset.loc[:,'id'] = bins_subset.index

        angle = self.b2s.get_rotation(rad=False)
        transform = [
            lambda x: create_square(x, size=self.b2s.converter.visiumHD.pixels_per_bin),
            lambda x: affinity.rotate(x, angle=angle, origin="center")
        ]

        cmap = PlotHD.get_cmap([0,1,0,0], [0,1,0,1])

        max_bbox = self.get_max_side(bbox)

        plot_params = {
            'color_by': gene,
            'edgecolor': 'white' if img else 'grey',
            'linewidth': math.exp(-0.01 * max_bbox / self.b2s.converter.visiumHD.pixels_per_bin),
        }

        if gene is not None and normalize:
            sc.pp.normalize_total(self.b2s.adata_hd, target_sum=1e4)

        self.slide.add_sample(Sample(bins_subset, self.b2s.adata_hd, transform=transform, cmap=cmap, **plot_params))

        if view_spots:
            if self.b2s.adata_v2 is None:
                self.b2s.convert(normalize=normalize)

            plot_params = {
                'color': [0.2, 1, 0.2, 0.5],
                'edgecolor': "green",
                'linewidth_by': gene
            }

            self.slide.add_sample(Sample(self.b2s.spot_polygons, self.b2s.adata_v2, **plot_params))

        if cells is not None:
            plot_params = {
                'color': [0.2, 0.2, 0.2, 0.5],
                'edgecolor': "green",
                'linewidth': 0.7
            }
            
            self.slide.add_sample(Sample(cells, **plot_params))

        self.slide.plot(bbox=bbox, img=img, gene=gene, view_not_sequenced=view_not_sequenced)

        return self.slide