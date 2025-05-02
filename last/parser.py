
import geopandas as gpd
from bokeh.models import ColumnDataSource

class GeoDataFrameToColumnDataSource:
    def __init__(self, gdf):
        self.gdf = gdf

    def generate_cds(self, bbox=None, offset=0.0, width=None, height=None, initial_resolution=156543.03392804097):
        if bbox is not None:
            if self.gdf.geometry.geom_type[0] == 'Point':
                subset_gdf = self.gdf.cx[bbox[0]-offset:bbox[1]+offset, bbox[2]-offset:bbox[3]+offset].copy()
            else:
                subset_gdf = self.gdf.cx[bbox[0]:bbox[1], bbox[2]:bbox[3]].copy()
        else:
            subset_gdf = self.gdf

        data = {}

        if subset_gdf.geometry.geom_type.iloc[0] == 'Point':
            data['x'] = subset_gdf.geometry.x.values.tolist()
            data['mapped_x'] = [x / (width / 256) for x in data['x']]
            data['y'] = subset_gdf.geometry.y.values.tolist()
            data['mapped_y'] = [y / (height / 256) for y in data['y']]
        elif subset_gdf.geometry.geom_type.iloc[0] == 'Polygon':
            data['x'] = [geom.centroid.x for geom in subset_gdf.geometry]
            data['mapped_x'] = [x / (width / 256) for x in data['x']]
            data['y'] = [geom.centroid.y for geom in subset_gdf.geometry]
            data['mapped_y'] = [y / (height / 256) for y in data['y']]
            data['radius'] = [geom.buffer(0).exterior.distance(geom.centroid) for geom in subset_gdf.geometry]
            data['mapped_radius'] = [geom.buffer(0).exterior.distance(geom.centroid) / (width / 256) for geom in subset_gdf.geometry]

        # Create a ColumnDataSource with the x, y, and size coordinates
        return ColumnDataSource(data=data)