"""HDViewer (ours): headless tasks on the real backend (run with the project's .venv).

Every task calls the code paths the web application runs:
  hd_viewer.core.pipeline   (reading, alignment, tiles, bin geometry)
  hd_viewer.core.workspace  (map layers)
  hd_viewer.api.views/frames (viewport -> binary frame sent to the browser)
  hd_viewer.core.clustering / embedding, hd_viewer.analysis.de
Rendering tasks are in tool_hdviewer_browser.py.
"""
import sys
from pathlib import Path
from types import SimpleNamespace

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH.parent))           # repo root: hd_viewer + benchmarks

from benchmarks.harness.worker import NotApplicable, main  # noqa: E402
from benchmarks.adapters import common  # noqa: E402

VIEW_PX = (1100, 860)       # map area of the 1400x900 window used by the browser runs


def _dataset(spec, bin_size=None):
    from hd_viewer.startup.data.VisiumHD import VisiumHD
    bin_size = bin_size or spec["bin_size"]
    return VisiumHD(img_path=spec["image_path"], adata_path=common.matrix_path(spec, bin_size),
                    positions_path=common.positions_path(spec, bin_size), bins_size=bin_size)


def _load_bins(spec, rec, bin_size=None, prefix=""):
    """The app's bins load (without tiles): returns (dataset, transform, b2s, sample)."""
    from hd_viewer.core import pipeline
    from hd_viewer.core.workspace import build_sample
    dataset = _dataset(spec, bin_size)
    with rec.phase(prefix + "read_counts_positions"):
        pipeline.read_inputs(dataset)
    with rec.phase(prefix + "fit_alignment"):
        transform = pipeline.fit_transform(dataset)
        processor, _ = pipeline.make_image_processor(dataset, transform, spec["workdir"])
    with rec.phase(prefix + "build_bin_geometry"):
        b2s = pipeline.build_b2s(dataset, transform, processor)
        squares = b2s.bin_squares
    with rec.phase(prefix + "build_layer") as p:
        sample = build_sample("bins", squares, b2s.adata_hd, b2s)
        p["n_obs"] = int(sample.n_geometries)
    return dataset, transform, b2s, sample


def _map_box(transform, box):
    """Full-resolution pixel box [x0, x1, y0, y1] -> HDViewer map coordinates."""
    import numpy as np
    x0, x1, y0, y1 = box
    corners = transform.points(np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], dtype=float))
    return [corners[:, 0].min(), corners[:, 0].max(), corners[:, 1].min(), corners[:, 1].max()]


# -- image --------------------------------------------------------------------------

def image_read(spec, rec):
    dataset = _dataset(spec)
    with rec.phase("read_image") as p:
        dataset.read_img()
        p["shape"] = list(dataset.img.shape)


def image_pyramid(spec, rec):
    """Align and tile the image as on first load (the cache is in a fresh work dir)."""
    from hd_viewer.core import pipeline
    dataset = _dataset(spec)
    with rec.phase("fit_alignment"):
        transform = pipeline.fit_transform(dataset)
    with rec.phase("read_image"):
        dataset.read_img()
    processor, single_image = pipeline.make_image_processor(dataset, transform, spec["workdir"])
    with rec.phase("warp_and_tile") as p:
        processor.create_tiles(overwrite=True, single_image=single_image)
        files = list(Path(processor.tiles_path).glob("z*"))
        p.update(tiles=len(files), max_zoom=processor.max_zoom, single_image=single_image,
                 bytes_on_disk=sum(f.stat().st_size for f in files))


def image_query(spec, rec):
    """Tiles a view needs: resolving and reading them from the cache (what the server sends)."""
    import math
    from hd_viewer.core import pipeline
    from hd_viewer.core.tiles import resolve_tile
    from hd_viewer.startup.data.image.image import register_tiles
    dataset = _dataset(spec)
    transform = pipeline.fit_transform(dataset)
    processor, single_image = pipeline.make_image_processor(dataset, transform, spec["workdir"])
    with rec.phase("setup_tiles"):
        processor.create_tiles(overwrite=False, single_image=single_image)
    key = register_tiles(processor.tiles_path)
    size = transform.output_size
    for name, box in common.views(spec):
        x0, x1, y0, y1 = _map_box(transform, box)
        with rec.phase(f"view:{name}") as p:
            if single_image:
                found = resolve_tile(key, 0, 0, 0)
                total = Path(found[0]).stat().st_size if found else 0
                n = 1
            else:
                # Zoom whose tiles show about one image pixel per screen pixel
                zoom = min(processor.max_zoom, max(0, round(math.log2(VIEW_PX[0] * size / max(x1 - x0, 1)
                                                                      / processor.tile_size))))
                tile_units = size / 2 ** zoom
                xs = range(max(0, int(x0 // tile_units)), min(2 ** zoom, int(x1 // tile_units) + 1))
                ys = range(max(0, int(y0 // tile_units)), min(2 ** zoom, int(y1 // tile_units) + 1))
                total = n = 0
                for tx in xs:
                    for ty in ys:
                        found = resolve_tile(key, zoom, tx, ty)
                        if found:
                            total += len(Path(found[0]).read_bytes())
                            n += 1
                p["zoom"] = zoom
            p.update(tiles=n, bytes=total)


# -- bins ---------------------------------------------------------------------------

def bins_load(spec, rec):
    _load_bins(spec, rec)


def _views_for(sample, kind="bins"):
    from hd_viewer.api.views import build_views
    loaded = SimpleNamespace(samples={kind: (kind, sample, kind.title())})
    return build_views(loaded)[kind]


def _query(spec, rec, sample, transform, kind):
    """Viewport -> data -> binary frame, coloured by the benchmark gene.

    Each named view uses a fresh layer view (a cold fetch of that region); the
    pan runs on one view, so it reuses what is already loaded, as the app does.
    """
    from hd_viewer.api import frames
    gene = spec["genes"][0]

    def show(view, box):
        x0, x1, y0, y1 = _map_box(transform, box)
        ops = view.update_viewport(x0, x1, y0, y1, *VIEW_PX)
        payload = frames.encode(*view.frame_parts(ops))
        return len(payload), len(view.rows.sink.data.get("x", ())) if hasattr(view.rows.sink, "data") else None

    for name, box in common.views(spec):
        view = _views_for(sample, kind)
        view.colour_by_column(gene)
        clear = getattr(sample, "_multiscale_cache", None)
        if clear is not None:
            clear.clear()
        with rec.phase(f"view:{name}") as p:
            p["bytes"], p["rows"] = show(view, box)
            p["level"] = str(view.rows.resolution_level) if hasattr(view.rows, "resolution_level") else None
    view = _views_for(sample, kind)
    view.colour_by_column(gene)
    show(view, spec["viewports"]["zoom_250um"])
    with rec.phase("pan") as p:
        sizes = [show(view, box)[0] for box in spec["viewports"]["pan"]]
        p.update(steps=len(sizes), bytes=sum(sizes))
    with rec.phase("gene_switch") as p:
        view.colour_by_column(spec["genes"][1])
        p["bytes"] = show(view, spec["viewports"]["pan"][-1])[0]


def bins_query(spec, rec):
    _, transform, _, sample = _load_bins(spec, rec, prefix="setup:")
    _query(spec, rec, sample, transform, "bins")


def bins_aggregate(spec, rec):
    """Multi-resolution: 2 µm bins summed into 8 and 16 µm blocks (factors 4 and 8)."""
    if spec["bin_size"] != 2:
        raise NotApplicable("aggregation is measured from 2 µm bins")
    _, transform, _, sample = _load_bins(spec, rec, prefix="setup:")
    gene = spec["genes"][0]
    bbox = _map_box(transform, spec["viewports"]["overview"])
    for factor in (4, 8):
        sample._multiscale_cache.clear()
        with rec.phase(f"grid:x{factor}"):
            sample.aggregation_grid(factor)
        with rec.phase(f"aggregate_counts:x{factor}") as p:
            agg = sample._get_aggregated_sample(factor)
            p["n_obs"] = int(agg.adata.n_obs)
        with rec.phase(f"fetch_overview:x{factor}") as p:
            data = sample.fetch_data(bbox=bbox, colnames=[gene], resolution_level=f"agg:{factor}")
            p["rows"] = len(data)


def _cluster(spec, rec, sample, prefix=""):
    from hd_viewer.core import embedding
    from hd_viewer.core.clustering import ClusteringState
    p = spec["analysis"]
    state = ClusteringState(sample)
    with rec.phase(prefix + "normalize_log1p_hvg"):
        state.prepare_adata_processed()
    with rec.phase(prefix + "pca"):
        state._ensure_pca(p["n_pcs"])
    with rec.phase(prefix + "neighbors"):
        state._ensure_neighbors(p["n_pcs"], p["n_neighbors"])
    with rec.phase(prefix + "leiden") as ph:
        key = state.run_clustering(resolution=p["resolution"], n_pcs=p["n_pcs"], n_neighbors=p["n_neighbors"])
        ph["clusters"] = int(sample.adata.obs[key].nunique())
    with rec.phase(prefix + "umap"):
        embedding.compute(sample, state, "umap", p["n_pcs"], p["n_neighbors"])
    return state


def _de(spec, rec, sample, positions):
    from hd_viewer.analysis.de import DEAnalysis
    groups = common.quadrant_groups(positions)
    masks = {g: (groups == g) for g in sorted(set(groups))}
    with rec.phase("wilcoxon_vs_rest") as p:
        result = DEAnalysis(sample).run(masks, mode=0)
        p["groups"] = len(masks)
    return result


def _bins_positions(sample, dataset):
    """Full-res (x, y) of each observation of the layer, in its adata order."""
    pos = dataset.df_tissue_positions.loc[sample.adata.obs_names]
    return pos["pxl_col_in_fullres"].to_numpy(), pos["pxl_row_in_fullres"].to_numpy()


def bins_cluster(spec, rec):
    _, _, _, sample = _load_bins(spec, rec, prefix="setup:")
    _cluster(spec, rec, sample)


def bins_de(spec, rec):
    dataset, _, _, sample = _load_bins(spec, rec, prefix="setup:")
    _de(spec, rec, sample, _bins_positions(sample, dataset))


# -- cells --------------------------------------------------------------------------

def _read_cells(spec, rec, transform, prefix=""):
    from hd_viewer.startup.data.image.segmentation.service import read_polygons
    with rec.phase(prefix + "read_polygons") as p:
        cells = read_polygons(spec["cells_path"])
        p["n_cells"] = len(cells)
    with rec.phase(prefix + "to_map_coordinates"):
        cells.geometry = transform.image_geometry(cells.geometry.values)
        cells["index"] = [f"cell_{i + 1}" for i in range(len(cells))]
        cells.index = cells["index"].values
    return cells


def cells_load(spec, rec):
    from hd_viewer.core import pipeline
    dataset = _dataset(spec, spec["cells_bin_size"])
    transform = pipeline.fit_transform(dataset)
    _read_cells(spec, rec, transform)


def _cells_layer(spec, rec, prefix=""):
    from hd_viewer.core.workspace import build_sample
    from hd_viewer.startup.operations.pooler import CellsPooler
    dataset, transform, b2s, _ = _load_bins(spec, rec, spec["cells_bin_size"], prefix="setup:bins:")
    cells = _read_cells(spec, rec, transform, prefix="setup:")
    with rec.phase(prefix + "pool_bins_into_cells") as p:
        adata = CellsPooler(b2s.bin_squares, b2s.adata_hd, cells, keep_cols=[]).run(predicate="intersects")
        p["n_cells_with_counts"] = int(adata.n_obs)
    with rec.phase(prefix + "build_layer"):
        sample = build_sample("cells", cells, adata, b2s)
    return dataset, transform, sample


def cells_pool(spec, rec):
    _cells_layer(spec, rec)


def cells_query(spec, rec):
    _, transform, sample = _cells_layer(spec, rec, prefix="setup:")
    _query(spec, rec, sample, transform, "cells")


def _cell_positions(sample, transform):
    import numpy as np
    centroids = sample.gdf.loc[sample.adata.obs_names].geometry.centroid
    xy = transform.points(np.column_stack([centroids.x, centroids.y]), inverse=True)
    return xy[:, 0], xy[:, 1]


def cells_cluster(spec, rec):
    _, _, sample = _cells_layer(spec, rec, prefix="setup:")
    _cluster(spec, rec, sample)


def cells_de(spec, rec):
    _, transform, sample = _cells_layer(spec, rec, prefix="setup:")
    _de(spec, rec, sample, _cell_positions(sample, transform))


TASKS = {
    "image.read": image_read, "image.pyramid": image_pyramid, "image.query": image_query,
    "bins.load": bins_load, "bins.query": bins_query, "bins.aggregate": bins_aggregate,
    "bins.cluster": bins_cluster, "bins.de": bins_de,
    "cells.load": cells_load, "cells.pool": cells_pool, "cells.query": cells_query,
    "cells.cluster": cells_cluster, "cells.de": cells_de,
}

if __name__ == "__main__":
    main(TASKS, versions=("scanpy", "anndata", "geopandas", "shapely", "numpy"))
