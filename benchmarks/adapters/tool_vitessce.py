"""Vitessce: data preparation in Python (bench-sd env) and rendering in headless Chromium.

Vitessce reads SpatialData Zarr stores over HTTP. Following its SpatialData
examples, the Space Ranger outputs are converted once into a Vitessce-ready store:
  * bins rasterized to a label image linked to the table (spatialdata.rasterize_bins,
    as in Vitessce's Visium HD example), shown as a segmentation layer,
  * the expression table as CSC, so one gene is read without the whole matrix,
  * the full-resolution image as a multiscale (OME-Zarr) pyramid,
  * cells as a multiscale label image (Vitessce's SpatialData segmentation type).
The store and the built page (vitessce_app/dist) are served by a local static
server; views are set through the config's coordination values (spatialTargetX/Y,
spatialZoom, featureSelection), which Vitessce applies in place.
"""
import json
import math
import os
import sys
import threading
import time
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH.parent))
sys.path.insert(0, str(BENCH / "adapters"))

import warnings  # noqa: E402
warnings.filterwarnings("ignore")

from benchmarks.harness.browser import BrowserSession, free_port  # noqa: E402
from benchmarks.harness.worker import NotApplicable, main  # noqa: E402
from benchmarks.adapters import common  # noqa: E402

APP_DIST = BENCH / "adapters" / "vitessce_app" / "dist"
VIEW_PX = (1400, 900)        # Vitessce fills the window: the spatial view gets most of it


def _sd():
    import tool_spatialdata
    return tool_spatialdata


def _formats():
    """Zarr v2 SpatialData formats: Vitessce 3.5 reads .zattrs/.zarray metadata, not Zarr v3."""
    from spatialdata._io.format import (RasterFormatV01, ShapesFormatV01, SpatialDataContainerFormatV01,
                                        TablesFormatV01)
    return [SpatialDataContainerFormatV01(), RasterFormatV01(), ShapesFormatV01(), TablesFormatV01()]


# -- data preparation ----------------------------------------------------------------

def _labels_from_cells(spec, sdata_cells, image_shape):
    """Cell polygons -> uint32 label image (label = instance id + 1), full resolution."""
    import numpy as np
    import shapely
    from skimage.draw import polygon as fill
    labels = np.zeros(image_shape[:2], dtype=np.uint32)
    geoms = sdata_cells.geometry.values
    for i, g in enumerate(geoms):
        xy = shapely.get_coordinates(shapely.get_exterior_ring(g))
        rr, cc = fill(xy[:, 1], xy[:, 0], shape=labels.shape)
        labels[rr, cc] = sdata_cells.index[i] + 1
    return labels


def convert_bins(spec, rec, bin_size, out, prefix=""):
    """Space Ranger -> Vitessce-ready SpatialData store, as in Vitessce's Visium HD example:
    table as CSC, bins rasterized to a Labels element linked to the table, multiscale image."""
    from spatialdata import rasterize_bins, rasterize_bins_link_table_to_labels
    sdm = _sd()
    table_name, shapes_name = sdm._table_name(bin_size), sdm._shapes_name(spec, bin_size)
    labels_name = f"rasterized_{int(bin_size):03d}um"
    with rec.phase(prefix + "read_visium_hd"):
        sdata = sdm._read(spec, [bin_size], multiscale=True)
    with rec.phase(prefix + "table_to_csc"):
        sdata.tables[table_name].X = sdata.tables[table_name].X.tocsc()
    with rec.phase(prefix + "rasterize_bins_to_labels"):
        sdata[labels_name] = rasterize_bins(sdata, shapes_name, table_name, "array_col", "array_row",
                                            return_region_as_labels=True)
        rasterize_bins_link_table_to_labels(sdata, table_name=table_name, rasterized_labels_name=labels_name)
    keep = [sdm._image_name(spec), labels_name, table_name]
    with rec.phase(prefix + "write_zarr") as p:
        sdata.subset(keep).write(out, overwrite=True, sdata_formats=_formats())
        p["bytes_on_disk"] = sum(f.stat().st_size for f in Path(out).rglob("*") if f.is_file())


def convert_cells(spec, rec, out, prefix=""):
    """Cells -> label image + pooled table (pooling as SpatialData does), with the image."""
    import numpy as np
    import spatialdata as sd
    from spatialdata.models import Labels2DModel, TableModel
    from spatialdata.transformations import Identity
    sdm = _sd()
    cs = sdm._cs(spec)
    source = sdm._zarr(spec, rec, "cells")          # cached: image + cells + pooled table
    image = source.images[sdm._image_name(spec)]
    shape = image["scale0"]["image"].shape[1:]
    with rec.phase(prefix + "rasterize_labels") as p:
        labels = _labels_from_cells(spec, source.shapes["cells"], shape)
        p["shape"] = list(labels.shape)
    with rec.phase(prefix + "build_store"):
        labels_el = Labels2DModel.parse(labels, dims=("y", "x"), scale_factors=sdm.SCALE_FACTORS,
                                        chunks=(1024, 1024), transformations={cs: Identity()})
        table = source.tables["cells_table"].copy()
        table.obs["region"] = "cells_labels"
        table.obs["region"] = table.obs["region"].astype("category")
        table.obs["instance_id"] = table.obs["instance_id"].astype(np.int64) + 1
        table.uns.pop("spatialdata_attrs", None)
        table = TableModel.parse(table, region="cells_labels", region_key="region", instance_key="instance_id")
        table.X = table.X.tocsc()
        out_sdata = sd.SpatialData(images={sdm._image_name(spec): image}, labels={"cells_labels": labels_el},
                                   tables={"cells_table": table})
    with rec.phase(prefix + "write_zarr") as p:
        out_sdata.write(out, overwrite=True, sdata_formats=_formats())
        p["bytes_on_disk"] = sum(f.stat().st_size for f in Path(out).rglob("*") if f.is_file())


def _store(spec, rec, kind, bin_size=None):
    """Cached Vitessce store (built in 'setup' phases on first use, written atomically)."""
    import shutil
    name = f"{spec['dataset']}_{kind}_{bin_size or spec['cells_bin_size']}.zarr"
    path = Path(spec["cache_dir"]).resolve() / "vitessce" / name
    if not path.is_dir():
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + f".tmp{os.getpid()}")
        shutil.rmtree(tmp, ignore_errors=True)
        if kind == "bins":
            convert_bins(spec, rec, bin_size, tmp, prefix="setup:")
        else:
            convert_cells(spec, rec, tmp, prefix="setup:")
        tmp.rename(path)
    return path


# -- config -----------------------------------------------------------------------------

def _zoom_for(box):
    x0, x1, y0, y1 = box
    return math.log2(min(VIEW_PX[0] / (x1 - x0), VIEW_PX[1] / (y1 - y0)))


def make_config(spec, url, kind, bin_size=None, show_image=True):
    """spatialBeta view: RGB image layer + a segmentation layer (bins or cells labels) coloured by gene."""
    from vitessce import CoordinationLevel as CL
    from vitessce import SpatialDataWrapper, VitessceConfig, get_initial_coordination_scope_prefix
    sdm = _sd()
    cs = sdm._cs(spec)
    vc = VitessceConfig(schema_version="1.0.17", name="benchmark")
    obs_type = "bin" if kind == "bins" else "cell"
    common_args = dict(sdata_url=url, image_path=f"images/{sdm._image_name(spec)}", coordinate_system=cs)
    if kind == "bins":
        t = sdm._table_name(bin_size)
        wrapper = SpatialDataWrapper(**common_args, table_path=f"tables/{t}", obs_feature_matrix_path=f"tables/{t}/X",
                                     labels_path=f"labels/rasterized_{int(bin_size):03d}um",
                                     coordination_values={"obsType": obs_type})
    elif kind == "cells":
        wrapper = SpatialDataWrapper(**common_args, table_path="tables/cells_table",
                                     obs_feature_matrix_path="tables/cells_table/X",
                                     labels_path="labels/cells_labels",
                                     coordination_values={"obsType": obs_type})
    else:
        wrapper = SpatialDataWrapper(**common_args)
    dataset = vc.add_dataset(name="data").add_object(wrapper)
    spatial = vc.add_view("spatialBeta", dataset=dataset)
    controller = vc.add_view("layerControllerBeta", dataset=dataset)
    box = spec["viewports"]["overview"]
    vc.link_views_by_dict([spatial, controller], {
        "spatialTargetX": (box[0] + box[1]) / 2, "spatialTargetY": (box[2] + box[3]) / 2,
        "spatialZoom": _zoom_for(box)}, meta=False)
    if show_image:
        vc.link_views_by_dict([spatial, controller], {
            "imageLayer": CL([{"photometricInterpretation": "RGB"}])},
            scope_prefix=get_initial_coordination_scope_prefix("A", "image"))
    if kind in ("bins", "cells"):
        vc.link_views_by_dict([spatial, controller], {
            "segmentationLayer": CL([{"segmentationChannel": CL([{
                "spatialChannelOpacity": 1.0, "obsColorEncoding": "geneSelection",
                "featureValueColormapRange": [0, 0.5]}])}])},
            scope_prefix=get_initial_coordination_scope_prefix("A", "obsSegmentations"))
        vc.link_views([spatial, controller], ["obsType", "featureSelection"], [obs_type, [spec["genes"][0]]])
    vc.layout(spatial)
    config = vc.to_dict(base_url="")
    # vitessce-python 3.5 (the release compatible with spatialdata's Zarr v3) writes the
    # 1.0.17 schema; the page runs Vitessce JS 3.9, whose SpatialData labels support
    # (with coordinate transforms) is configured as 'obsSegmentations' + tablePath, as
    # vitessce-python 3.8 writes it.
    config["version"] = "1.0.18"
    if kind in ("bins", "cells"):
        options = config["datasets"][0]["files"][0]["options"]
        labels = options.pop("labels", None) or options.pop("obsSegmentations", None) or {}
        table = f"tables/{sdm._table_name(bin_size)}" if kind == "bins" else "tables/cells_table"
        options["obsSegmentations"] = dict(path=labels.get("path") or (
            f"labels/rasterized_{int(bin_size):03d}um" if kind == "bins" else "labels/cells_labels"),
            tablePath=table, coordinateSystem=cs)
        options.pop("obsSpots", None)
    return config


# -- serving and driving ----------------------------------------------------------------

class _Handler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()

    def log_message(self, *args):
        pass


class StaticServer:
    """Serves the built page at / and data folders at /data/<name>/ (one thread)."""

    def __init__(self, root):
        self.root, self.port = root, free_port()

    def __enter__(self):
        self.httpd = ThreadingHTTPServer(("127.0.0.1", self.port), partial(_Handler, directory=str(self.root)))
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.url = f"http://127.0.0.1:{self.port}"
        return self

    def __exit__(self, *exc):
        self.httpd.shutdown()


def _site(spec, store, config_dict):
    """A folder holding the page, the config and a link to the store."""
    site = Path(spec["workdir"]) / "site"
    site.mkdir(exist_ok=True)
    if not APP_DIST.exists():
        raise RuntimeError(f"Build the Vitessce page first: {APP_DIST} (see benchmarks/setup.sh)")
    for item in APP_DIST.iterdir():
        target = site / item.name
        if not target.exists():
            os.symlink(item, target)
    if store is not None:
        link = site / "data"
        if link.is_symlink() or link.exists():
            link.unlink()
        os.symlink(Path(store).resolve(), link)
    (site / "config.json").write_text(json.dumps(config_dict))
    return site


READY = """() => document.querySelectorAll('canvas').length > 0
  && !document.querySelector('[class*=loadingIndicatorContainer]')"""


def _wait_drawn(b):
    """Vitessce has mounted the deck.gl view and hidden its loading indicator
    (shown while it fetches and decodes data). Network quiet is checked after."""
    b.page.wait_for_timeout(50)
    b.page.wait_for_function(READY, timeout=3_600_000, polling=20)


def _act(b, script, arg):
    def run():
        b.page.evaluate(script, arg)
        _wait_drawn(b)
    return run


def _session(spec, rec, kind):
    bin_size = spec["bin_size"] if kind == "bins" else None
    if kind == "bins":
        store = _store(spec, rec, "bins", bin_size)
    elif kind == "cells":
        store = _store(spec, rec, "cells")
    else:                                   # the image of any bins store
        store = _store(spec, rec, "bins", spec["bin_sizes"][-1])
    with StaticServer(Path(spec["workdir"]) / "site") as server:
        config = make_config(spec, server.url + "/data", kind if kind != "image" else "image", bin_size)
        _site(spec, store, config)
        with BrowserSession(ignore=()) as b:
            def load():
                b.page.goto(f"{server.url}/index.html?config={server.url}/config.json")
                _wait_drawn(b)
            seconds, n, nbytes = b.timed(load, timeout=3600, quiet=1.0)
            rec.add("load_to_first_view", seconds, requests=n, bytes=nbytes)
            for name, box in common.views(spec):
                x, y, z = (box[0] + box[1]) / 2, (box[2] + box[3]) / 2, _zoom_for(box)
                seconds, n, nbytes = b.timed(_act(b, "([x, y, z]) => window.__setView(x, y, z)", [x, y, z]),
                                             quiet=1.0)
                rec.add(f"view:{name}", seconds, requests=n, bytes=nbytes)
            z = _zoom_for(spec["viewports"]["zoom_250um"])
            b.frames_start()
            steps, total_n, total_bytes = [], 0, 0
            for box in spec["viewports"]["pan"]:
                x, y = (box[0] + box[1]) / 2, (box[2] + box[3]) / 2
                seconds, n, nbytes = b.timed(_act(b, "([x, y, z]) => window.__setView(x, y, z)", [x, y, z]),
                                             quiet=1.0)
                steps.append(seconds)
                total_n += n
                total_bytes += nbytes
            frames = b.frames_stop()
            rec.add("pan", sum(steps), steps=len(steps), step_seconds=steps, requests=total_n, bytes=total_bytes,
                    **frames)
            if kind in ("bins", "cells"):
                seconds, n, nbytes = b.timed(_act(b, "(g) => window.__setGene(g)", spec["genes"][1]), quiet=1.0)
                rec.add("gene_switch", seconds, requests=n, bytes=nbytes)
            if spec.get("snapshots"):
                b.page.screenshot(path=str(Path(spec["workdir"]) / f"vitessce_{kind}.png"))
            rec.note(js_heap=b.js_heap(), page_errors=b.failed[:5])


def bins_load(spec, rec):
    convert_bins(spec, rec, spec["bin_size"], Path(spec["workdir"]) / "vitessce_bins.zarr")


def cells_load(spec, rec):
    convert_cells(spec, rec, Path(spec["workdir"]) / "vitessce_cells.zarr")


def e2e_first_view(spec, rec):
    """Convert the Space Ranger outputs for Vitessce, then open the first view."""
    out = Path(spec["workdir"]) / "e2e.zarr"
    convert_bins(spec, rec, spec["bin_size"], out, prefix="convert:")
    with StaticServer(Path(spec["workdir"]) / "site") as server:
        _site(spec, out, make_config(spec, server.url + "/data", "bins", spec["bin_size"]))
        with BrowserSession(ignore=()) as b:
            def load():
                b.page.goto(f"{server.url}/index.html?config={server.url}/config.json")
                _wait_drawn(b)
            seconds, n, nbytes = b.timed(load, timeout=3600, quiet=1.0)
            rec.add("open_first_view", seconds, requests=n, bytes=nbytes)
    convert = sum(p["seconds"] for p in rec.phases if p["name"].startswith("convert:"))
    rec.add("load_to_first_view", convert + seconds)


def _na(reason):
    def task(spec, rec):
        raise NotApplicable(reason)
    return task


TASKS = {
    "image.render": lambda spec, rec: _session(spec, rec, "image"),
    "bins.render": lambda spec, rec: _session(spec, rec, "bins"),
    "cells.render": lambda spec, rec: _session(spec, rec, "cells"),
    "bins.load": bins_load, "cells.load": cells_load, "e2e.first_view": e2e_first_view,
    "image.read": _na("the browser reads image tiles from the OME-Zarr pyramid; no whole-image read"),
    "image.pyramid": _na("uses the SpatialData/OME-Zarr pyramid (see SpatialData image.pyramid)"),
    "image.query": _na("tiles are fetched by the browser; measured in image.render"),
    "bins.query": _na("observations are fetched whole by the browser; measured in bins.render"),
    "bins.aggregate": _na("no multi-resolution aggregation of bins"),
    "bins.cluster": _na("visualisation only; no analysis"), "bins.de": _na("visualisation only; no analysis"),
    "cells.pool": _na("no bin-to-cell aggregation (uses SpatialData's)"),
    "cells.query": _na("measured in cells.render"),
    "cells.cluster": _na("visualisation only; no analysis"), "cells.de": _na("visualisation only; no analysis"),
}

if __name__ == "__main__":
    main(TASKS, versions=("vitessce", "spatialdata", "spatialdata-io", "playwright"))
