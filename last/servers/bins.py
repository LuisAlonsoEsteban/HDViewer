
from flask import Flask, request, jsonify
import geopandas as gpd

class BinsRenderer:

    def __init__(self, gdf):
        self.gdf = gdf

    def fetch_data(self, left, right, top, bottom):
        data = []
        cropped_gdf = self.gdf.cx[left:right, bottom:top].copy()
        data = {
            'x': [geom.x for geom in cropped_gdf.geometry if geom.geom_type == 'Point'], 
            'y': [geom.y for geom in cropped_gdf.geometry if geom.geom_type == 'Point']
        }
        # TODO: Add other data as needed
        return data


if __name__ == "__main__":
    gdf = gpd.read_file("/scratch/lalonsoeste/PhD/VisiumHD2Visium/dataframe.geojson")
    br = BinsRenderer(gdf)

    app = Flask("BinsRenderer")

    @app.route('/get_data')
    def get_data():
        left = request.args.get('left', type=float)
        right = request.args.get('right', type=float)
        top = request.args.get('top', type=float)
        bottom = request.args.get('bottom', type=float)

        data = br.fetch_data(left, right, top, bottom)
        return jsonify(data), {'Access-Control-Allow-Origin':'*'}

    app.run(debug=True, host='0.0.0.0', port=5001)
