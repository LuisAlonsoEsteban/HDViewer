"""napari: run with the bench-napari env, on an X server (Xvfb + Mesa llvmpipe OpenGL).

Visium HD reaches napari through SpatialData (napari-spatialdata reads the same
Zarr store), so data come from the SpatialData cache built by
tool_spatialdata._zarr. Layers are added the way napari-spatialdata does:
the image as a multiscale RGB layer, bins as a Points layer (square symbols of
the bin size: napari's fastest option; a Shapes layer of millions of squares
does not load), and cells as a Shapes layer of polygons, coloured by a gene.

A view is rendered when `viewer.screenshot()` returns: napari draws
synchronously, including fetching the multiscale tiles the view needs.
"""
import os
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH.parent))
sys.path.insert(0, str(BENCH / "adapters"))

import warnings  # noqa: E402
warnings.filterwarnings("ignore")

from benchmarks.adapters import common  # noqa: E402
from benchmarks.harness.worker import main  # noqa: E402

CANVAS = (1100, 860)        # width, height: the map area of the browser runs


def _sd():
    import tool_spatialdata
    return tool_spatialdata


def _viewer(rec, phase="setup:start_viewer"):
    """A napari viewer; its start-up (Qt, OpenGL) is a setup phase except end to end."""
    with rec.phase(phase):
        import napari
        viewer = napari.Viewer(show=False)
        viewer.window._qt_window.resize(CANVAS[0] + 300, CANVAS[1] + 100)
        viewer.window._qt_viewer.canvas.size = (CANVAS[1], CANVAS[0])
    return viewer


def _shot(viewer):
    return viewer.screenshot(canvas_only=True, flash=False)


def _goto(viewer, box):
    x0, x1, y0, y1 = box
    height, width = viewer.window._qt_viewer.canvas.size
    viewer.camera.center = ((y0 + y1) / 2, (x0 + x1) / 2)
    viewer.camera.zoom = min(width / (x1 - x0), height / (y1 - y0))


def _add_image(viewer, sdata, spec):
    image = sdata.images[f"{spec['dataset_id']}_full_image"]
    levels = [image[k]["image"].data.transpose(1, 2, 0) for k in image.keys()]   # (y, x, c)
    return viewer.add_image(levels, multiscale=True, rgb=True, name="image")


def _gene(table, gene):
    import numpy as np
    column = table[:, gene].X
    return np.asarray(column.todense() if hasattr(column, "todense") else column).ravel()


def _session(spec, rec, viewer, sdata, layer=None, table=None, recolour=None):
    genes = spec["genes"]
    for name, box in common.views(spec):
        _goto(viewer, box)
        with rec.phase(f"view:{name}"):
            _shot(viewer)
    _goto(viewer, spec["viewports"]["zoom_250um"])
    _shot(viewer)
    with rec.phase("pan") as p:
        for box in spec["viewports"]["pan"]:
            _goto(viewer, box)
            _shot(viewer)
        p["steps"] = len(spec["viewports"]["pan"])
    if recolour is not None:
        with rec.phase("gene_switch"):
            recolour(_gene(table, genes[1]))
            _shot(viewer)
    if spec.get("snapshots"):
        from PIL import Image
        kind = "image" if layer is None else layer.name
        Image.fromarray(_shot(viewer)).save(Path(spec["workdir"]) / f"napari_{kind}.png")


def image_read(spec, rec):
    viewer = _viewer(rec)
    with rec.phase("read_image") as p:
        layers = viewer.open(spec["image_path"], plugin="napari")
        p["shape"] = list(layers[0].data.shape)


def image_render(spec, rec):
    sdata = _sd()._zarr(spec, rec, "bins", spec["bin_sizes"][-1])
    viewer = _viewer(rec)
    with rec.phase("add_layers"):
        _add_image(viewer, sdata, spec)
    _session(spec, rec, viewer, sdata)


def _bins_points(viewer, sdata, spec, bin_size):
    import shapely
    sd = _sd()
    shapes = sdata.shapes[sd._shapes_name(spec, bin_size)]
    table = sdata.tables[sd._table_name(bin_size)]
    instance = table.obs[table.uns["spatialdata_attrs"]["instance_key"]].to_numpy()
    centroids = shapely.centroid(shapes.geometry.loc[instance].values)
    yx = shapely.get_coordinates(centroids)[:, ::-1]
    side = common.bin_area_px(spec, bin_size)
    layer = viewer.add_points(yx, size=side, symbol="square", border_width=0,
                              features={"gene": _gene(table, spec["genes"][0])}, face_color="gene",
                              face_colormap="viridis", name="bins")
    return layer, table


def _recolour(layer):
    def apply(values):
        layer.features = {"gene": values}
        layer.face_color = "gene"
        layer.refresh()
    return apply


def bins_render(spec, rec):
    b = spec["bin_size"]
    sdata = _sd()._zarr(spec, rec, "bins", b)
    viewer = _viewer(rec)
    with rec.phase("add_layers") as p:
        _add_image(viewer, sdata, spec)
        layer, table = _bins_points(viewer, sdata, spec, b)
        p["n_obs"] = int(table.n_obs)
    _session(spec, rec, viewer, sdata, layer, table, _recolour(layer))


def _cells_shapes(viewer, sdata):
    import numpy as np
    import shapely
    table = sdata.tables["cells_table"]
    keep = np.asarray(table.X.sum(axis=1)).ravel() > 0
    table = table[keep]
    geoms = sdata.shapes["cells"].geometry.loc[table.obs["instance_id"].to_numpy()].values
    polygons = [shapely.get_coordinates(shapely.get_exterior_ring(g))[:-1, ::-1] for g in geoms]
    return polygons, table


def cells_render(spec, rec):
    import matplotlib
    import numpy as np
    sdata = _sd()._zarr(spec, rec, "cells")
    viewer = _viewer(rec)
    cmap = matplotlib.colormaps["viridis"]

    def colours(values):
        values = np.asarray(values, dtype=float)
        top = values.max() or 1.0
        return cmap(values / top)

    with rec.phase("add_layers") as p:
        _add_image(viewer, sdata, spec)
        polygons, table = _cells_shapes(viewer, sdata)
        layer = viewer.add_shapes(polygons, shape_type="polygon", edge_width=0,
                                  face_color=colours(_gene(table, spec["genes"][0])), name="cells")
        p["n_cells"] = len(polygons)

    def recolour(values):
        layer.face_color = colours(values)
        layer.refresh()

    _session(spec, rec, viewer, sdata, layer, table, recolour)


def e2e_first_view(spec, rec):
    """Read the Space Ranger outputs (as napari-spatialdata does), show image + bins."""
    import time
    b = spec["bin_size"]
    started = time.perf_counter()
    viewer = _viewer(rec, phase="start_viewer")
    sdata = _sd()._read(spec, [b], multiscale=True)
    _add_image(viewer, sdata, spec)
    _bins_points(viewer, sdata, spec, b)
    _goto(viewer, spec["viewports"]["overview"])
    _shot(viewer)
    rec.add("load_to_first_view", time.perf_counter() - started)


TASKS = {
    "image.read": image_read, "image.render": image_render,
    "bins.render": bins_render, "cells.render": cells_render,
    "e2e.first_view": e2e_first_view,
}

NOT_APPLICABLE = {
    "image.pyramid": "napari displays a precomputed pyramid (SpatialData/OME-Zarr); it does not build one",
    "image.query": "no data-query step: the image layer slices what the canvas shows",
    "bins.load": "reading is done by spatialdata-io (napari-spatialdata), see SpatialData",
    "bins.query": "all observations are held in the layer; no viewport query",
    "bins.aggregate": "no multi-resolution aggregation of bins",
    "bins.cluster": "no built-in analysis (users call scanpy, see SpatialData)",
    "bins.de": "no built-in analysis (users call scanpy, see SpatialData)",
    "cells.load": "reading is done by spatialdata/geopandas, see SpatialData",
    "cells.pool": "no bin-to-cell aggregation",
    "cells.query": "all observations are held in the layer; no viewport query",
    "cells.cluster": "no built-in analysis", "cells.de": "no built-in analysis",
}

if __name__ == "__main__":
    from benchmarks.harness.worker import NotApplicable

    def na(reason):
        def task(spec, rec):
            raise NotApplicable(reason)
        return task

    main({**{k: na(v) for k, v in NOT_APPLICABLE.items()}, **TASKS},
         versions=("napari", "vispy", "PySide6", "spatialdata"))
