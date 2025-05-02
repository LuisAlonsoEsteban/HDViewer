
from flask import Flask, send_file, abort
import os

app = Flask(__name__)

@app.route("/<int:zoom>/<int:x>/<int:y>")
def get_tile(zoom, x, y):
    # Assuming tiles are named based on a fixed format
    base_path = '/scratch/lalonsoeste/PhD/SpatialTranscriptomics/data/spatial/10x/VisiumHD/Mouse intestine/tiles'
    filename = f'z{zoom}_x{x}_y{y}.png'
    filepath = os.path.join(base_path, filename)

    if not os.path.exists(filepath):
        abort(404)

    return send_file(filepath, mimetype='image/png')

if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=5000)

def run(port=5000):
    app.run(debug=True, host='0.0.0.0', port=port)
