"""SpatialData (spatialdata, spatialdata-io, spatialdata-plot): run with the bench-sd env.

The documented SpatialData workflow is followed:
  read Space Ranger outputs with spatialdata_io.visium_hd -> write/read Zarr ->
  query with spatialdata.bounding_box_query -> plot with spatialdata-plot
  (matplotlib, rendered off-screen at the browser's map size) -> analyse the
  table with scanpy -> aggregate bins with spatialdata.aggregate / rasterize_bins.
Analysis parameters match HDViewer's (min cells, target sum, HVG, PCs, k, resolution).

Rendering, query and analysis tasks start from a Zarr store that is built once
per dataset/bin size in the shared cache (a 'setup' phase, not part of the
task's figures), because that is how SpatialData users work after the first read.
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH.parent))

import warnings  # noqa: E402
warnings.filterwarnings("ignore")

from benchmarks.adapters import common  # noqa: E402
from benchmarks.harness.worker import NotApplicable, main  # noqa: E402

FIG_INCHES, DPI = (11.0, 8.6), 100       # 1100 x 860 px, the map area of the browser runs
SCALE_FACTORS = [2, 2, 2, 2, 2, 2]


def _cs(spec):
    """Coordinate system of the full-resolution pixels. spatialdata-io names it after the
    dataset; it is renamed 'global' (spatialdata's default, which spatialdata-plot's
    datashader path assumes)."""
    return "global"


def _shapes_name(spec, bin_size):
    return f"{spec['dataset_id']}_square_{int(bin_size):03d}um"


def _table_name(bin_size):
    return f"square_{int(bin_size):03d}um"


def _read(spec, bin_sizes, multiscale=False, **kwargs):
    from spatialdata_io import visium_hd
    if multiscale:
        kwargs["image_models_kwargs"] = dict(scale_factors=SCALE_FACTORS, chunks=(3, 1024, 1024))
    sdata = visium_hd(spec["sd_root"], dataset_id=spec["dataset_id"], bin_size=list(bin_sizes), **kwargs)
    sdata.rename_coordinate_systems({spec["dataset_id"]: _cs(spec)})
    return sdata


def _read_cells(spec):
    import geopandas as gpd
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import Identity
    gdf = gpd.read_file(spec["cells_path"]).set_crs(None, allow_override=True)
    gdf = gdf[["geometry"]].reset_index(drop=True)
    return ShapesModel.parse(gdf, transformations={_cs(spec): Identity()})


def _pool(spec, sdata, cells_name="cells"):
    import spatialdata as sd
    b = spec["cells_bin_size"]
    table = sdata.tables[_table_name(b)]
    return sd.aggregate(values_sdata=sdata, values=_shapes_name(spec, b), by_sdata=sdata, by=cells_name,
                        value_key=list(table.var_names), agg_func="sum", target_coordinate_system=_cs(spec),
                        table_name=_table_name(b))


def _zarr(spec, rec, kind, bin_size=None):
    """Cached Zarr store: kind 'bins' (bins + multiscale image) or 'cells' (cells + pooled table)."""
    import spatialdata as sd
    name = f"{spec['dataset']}_{kind}_{bin_size or spec['cells_bin_size']}.zarr"
    path = Path(spec["cache_dir"]) / "spatialdata" / name
    if not _complete(path):
        with rec.phase("setup:build_zarr_cache"):
            if kind == "bins":
                sdata = _read(spec, [bin_size], multiscale=True)
            else:
                sdata = _read(spec, [spec["cells_bin_size"]], multiscale=True)
                sdata["cells"] = _read_cells(spec)
                pooled = _pool(spec, sdata)
                table = pooled.tables[_table_name(spec["cells_bin_size"])]
                del sdata.tables[_table_name(spec["cells_bin_size"])]
                del sdata.shapes[_shapes_name(spec, spec["cells_bin_size"])]
                sdata.tables["cells_table"] = table
            write_atomic(sdata, path)
    with rec.phase("setup:read_zarr"):
        return sd.read_zarr(path)


def _complete(path):
    """A finished store (written atomically, so a present folder is complete)."""
    return Path(path).is_dir()


def write_atomic(sdata, path, **kwargs):
    """Write next to `path`, then rename: an interrupted run never leaves a partial cache."""
    import shutil
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    shutil.rmtree(tmp, ignore_errors=True)
    sdata.write(tmp, overwrite=True, **kwargs)
    shutil.rmtree(path, ignore_errors=True)
    tmp.rename(path)


def _image_name(spec):
    return f"{spec['dataset_id']}_full_image"


def _bbox(sdata, box, cs):
    import spatialdata as sd
    x0, x1, y0, y1 = box
    kwargs = dict(filter_table=True) if isinstance(sdata, sd.SpatialData) else {}
    return sd.bounding_box_query(sdata, axes=("x", "y"), min_coordinate=[x0, y0], max_coordinate=[x1, y1],
                                 target_coordinate_system=cs, **kwargs)


def _gene_values(table, gene):
    import numpy as np
    column = table[:, gene].X
    return np.asarray(column.todense() if hasattr(column, "todense") else column).ravel()


# -- image --------------------------------------------------------------------------

def _imread_model(spec, **kwargs):
    from dask_image.imread import imread
    from spatialdata.models import Image2DModel
    from spatialdata.transformations import Identity
    data = imread(spec["image_path"])
    if data.ndim == 4:
        data = data[0]
    return Image2DModel.parse(data, dims=("y", "x", "c"), transformations={_cs(spec): Identity()}, **kwargs)


def image_read(spec, rec):
    """As spatialdata-io loads the full-resolution image (dask), then into memory."""
    with rec.phase("read_lazy"):
        image = _imread_model(spec)
    with rec.phase("read_image") as p:
        array = image.data.compute()
        p["shape"] = list(array.shape)


def image_pyramid(spec, rec):
    """Multiscale image written to Zarr (what napari/Vitessce/spatialdata-plot then read)."""
    import spatialdata as sd
    with rec.phase("build_multiscale_and_write_zarr") as p:
        image = _imread_model(spec, scale_factors=SCALE_FACTORS, chunks=(3, 1024, 1024))
        path = Path(spec["workdir"]) / "image.zarr"
        sd.SpatialData(images={"image": image}).write(path, overwrite=True)
        p["bytes_on_disk"] = sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def _image_level(sdata_image, box, width_px=FIG_INCHES[0] * DPI):
    """The pyramid level a viewer would use for this view (about 1 image px per screen px)."""
    import math
    scales = list(sdata_image.keys())
    want = max(0, int(math.floor(math.log2(max((box[1] - box[0]) / width_px, 1)))))
    return scales[min(want, len(scales) - 1)]


def image_query(spec, rec):
    sdata = _zarr(spec, rec, "bins", spec["bin_sizes"][-1])
    image = sdata.images[_image_name(spec)]
    for name, box in common.views(spec):
        with rec.phase(f"view:{name}") as p:
            cropped = _bbox(image, box, _cs(spec))
            level = _image_level(cropped, box)
            pixels = cropped[level]["image"].data.compute()
            p.update(level=level, shape=list(pixels.shape))


def _render(sdata, spec, box, image=True, shapes=None, color=None, table=None, save=None):
    """One view: query the box, then spatialdata-plot draws image and shapes (Agg canvas)."""
    import matplotlib.pyplot as plt
    import spatialdata_plot  # noqa: F401  (registers the .pl accessor)
    view = _bbox(sdata, box, _cs(spec))
    fig, ax = plt.subplots(figsize=FIG_INCHES, dpi=DPI)
    plot = view
    if image:
        plot = plot.pl.render_images(_image_name(spec))
    if shapes and shapes in view.shapes and len(view.shapes[shapes]):
        plot = plot.pl.render_shapes(shapes, color=color, table_name=table)
    plot.pl.show(coordinate_systems=_cs(spec), ax=ax)
    fig.canvas.draw()
    if save:
        fig.savefig(save)
    plt.close(fig)


def _render_session(spec, rec, sdata, shapes=None, table=None):
    genes = spec["genes"]
    for name, box in common.views(spec):
        with rec.phase(f"view:{name}"):
            _render(sdata, spec, box, shapes=shapes, color=genes[0] if shapes else None, table=table)
    with rec.phase("pan") as p:
        for box in spec["viewports"]["pan"]:
            _render(sdata, spec, box, shapes=shapes, color=genes[0] if shapes else None, table=table)
        p["steps"] = len(spec["viewports"]["pan"])
    if shapes:
        with rec.phase("gene_switch"):
            _render(sdata, spec, spec["viewports"]["pan"][-1], shapes=shapes, color=genes[1], table=table)
    if spec.get("snapshots"):
        _render(sdata, spec, spec["viewports"]["pan"][-1], shapes=shapes, color=genes[1] if shapes else None,
                table=table, save=str(Path(spec["workdir"]) / f"spatialdata_{'image' if not shapes else 'cells' if shapes == 'cells' else 'bins'}.png"))


def image_render(spec, rec):
    sdata = _zarr(spec, rec, "bins", spec["bin_sizes"][-1])
    _render_session(spec, rec, sdata)


# -- bins ---------------------------------------------------------------------------

def bins_load(spec, rec):
    import spatialdata as sd
    with rec.phase("read_visium_hd") as p:
        sdata = _read(spec, [spec["bin_size"]])
        p["n_obs"] = int(sdata.tables[_table_name(spec["bin_size"])].n_obs)
    path = Path(spec["workdir"]) / "bins.zarr"
    with rec.phase("write_zarr"):
        sdata.write(path, overwrite=True)
    with rec.phase("read_zarr"):
        sd.read_zarr(path)


def _query(spec, rec, sdata, shapes, table_name):
    sub = sdata.subset([shapes, table_name], filter_tables=True)
    genes = spec["genes"]

    def show(box, gene):
        view = _bbox(sub, box, _cs(spec))
        if view is None or table_name not in view.tables:
            return 0
        values = _gene_values(view.tables[table_name], gene)
        geoms = view.shapes[shapes].geometry
        return len(values) + 0 * len(geoms)

    for name, box in common.views(spec):
        with rec.phase(f"view:{name}") as p:
            p["rows"] = show(box, genes[0])
    with rec.phase("pan") as p:
        p["rows"] = sum(show(box, genes[0]) for box in spec["viewports"]["pan"])
    with rec.phase("gene_switch") as p:
        p["rows"] = show(spec["viewports"]["pan"][-1], genes[1])


def bins_query(spec, rec):
    b = spec["bin_size"]
    sdata = _zarr(spec, rec, "bins", b)
    _query(spec, rec, sdata, _shapes_name(spec, b), _table_name(b))


def bins_render(spec, rec):
    b = spec["bin_size"]
    sdata = _zarr(spec, rec, "bins", b)
    _render_session(spec, rec, sdata, shapes=_shapes_name(spec, b), table=_table_name(b))


def bins_aggregate(spec, rec):
    """2 µm bins -> 8 and 16 µm (spatialdata.aggregate), and the 2 µm raster (rasterize_bins)."""
    import spatialdata as sd
    if spec["bin_size"] != 2:
        raise NotApplicable("aggregation is measured from 2 µm bins")
    with rec.phase("setup:read_visium_hd"):
        sdata = _read(spec, [2, 8, 16])
    table = sdata.tables[_table_name(2)]
    for size in (8, 16):
        with rec.phase(f"aggregate_counts:to_{size}um") as p:
            out = sd.aggregate(values_sdata=sdata, values=_shapes_name(spec, 2), by_sdata=sdata,
                               by=_shapes_name(spec, size), value_key=list(table.var_names), agg_func="sum",
                               target_coordinate_system=_cs(spec), table_name=_table_name(2))
            p["n_obs"] = int(out.tables[_table_name(2)].n_obs)
    with rec.phase("rasterize_bins:2um") as p:
        table.X = table.X.tocsc()
        raster = sd.rasterize_bins(sdata, _shapes_name(spec, 2), _table_name(2), "array_col", "array_row",
                                   value_key=[spec["genes"][0]])
        p["shape"] = list(raster.shape)


def _cluster(spec, rec, adata):
    import scanpy as sc
    a = spec["analysis"]
    with rec.phase("normalize_log1p_hvg"):
        adata = adata.copy()
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable].copy()
    with rec.phase("pca"):
        sc.pp.pca(adata, n_comps=a["n_pcs"])
    with rec.phase("neighbors"):
        sc.pp.neighbors(adata, n_neighbors=a["n_neighbors"], n_pcs=a["n_pcs"])
    with rec.phase("leiden") as p:
        sc.tl.leiden(adata, resolution=a["resolution"], flavor="igraph", n_iterations=2, directed=False,
                     random_state=0)
        p["clusters"] = int(adata.obs["leiden"].nunique())
    with rec.phase("umap"):
        sc.tl.umap(adata)


def _de(spec, rec, adata, xy):
    import scanpy as sc
    adata = adata.copy()
    adata.obs["quadrant"] = common.quadrant_groups(xy)
    with rec.phase("wilcoxon_vs_rest") as p:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, "quadrant", method="wilcoxon", reference="rest", pts=True)
        p["groups"] = int(adata.obs["quadrant"].nunique())


def _bins_table(spec, rec):
    with rec.phase("setup:read_visium_hd"):
        sdata = _read(spec, [spec["bin_size"]])
    return sdata.tables[_table_name(spec["bin_size"])]


def bins_cluster(spec, rec):
    _cluster(spec, rec, _bins_table(spec, rec))


def bins_de(spec, rec):
    table = _bins_table(spec, rec)
    _de(spec, rec, table, (table.obs["pxl_col_in_fullres"].to_numpy() if "pxl_col_in_fullres" in table.obs else
                           table.obsm["spatial"][:, 0],
                           table.obs["pxl_row_in_fullres"].to_numpy() if "pxl_row_in_fullres" in table.obs else
                           table.obsm["spatial"][:, 1]))


# -- cells --------------------------------------------------------------------------

def cells_load(spec, rec):
    with rec.phase("read_polygons") as p:
        cells = _read_cells(spec)
        p["n_cells"] = len(cells)


def cells_pool(spec, rec):
    with rec.phase("setup:read_visium_hd"):
        sdata = _read(spec, [spec["cells_bin_size"]])
    with rec.phase("setup:read_polygons"):
        sdata["cells"] = _read_cells(spec)
    with rec.phase("pool_bins_into_cells") as p:
        pooled = _pool(spec, sdata)
        p["n_cells_with_counts"] = int((pooled.tables[_table_name(spec["cells_bin_size"])].X.sum(axis=1) > 0).sum())


def cells_query(spec, rec):
    sdata = _zarr(spec, rec, "cells")
    _query(spec, rec, sdata, "cells", "cells_table")


def cells_render(spec, rec):
    sdata = _zarr(spec, rec, "cells")
    _render_session(spec, rec, sdata, shapes="cells", table="cells_table")


def _cells_table(spec, rec):
    sdata = _zarr(spec, rec, "cells")
    table = sdata.tables["cells_table"]
    centroids = sdata.shapes["cells"].geometry.centroid
    ids = table.obs["instance_id"].to_numpy()
    xy = (centroids.loc[ids].x.to_numpy(), centroids.loc[ids].y.to_numpy())
    return table, xy


def cells_cluster(spec, rec):
    table, _ = _cells_table(spec, rec)
    _cluster(spec, rec, table[table.X.sum(axis=1).A1 > 0] if hasattr(table.X, "toarray") else table)


def cells_de(spec, rec):
    import numpy as np
    table, (x, y) = _cells_table(spec, rec)
    keep = np.asarray(table.X.sum(axis=1)).ravel() > 0
    _de(spec, rec, table[keep], (x[keep], y[keep]))


def e2e_first_view(spec, rec):
    """Read the Space Ranger outputs and draw the first (overview) plot."""
    b = spec["bin_size"]
    with rec.phase("load_to_first_view"):
        sdata = _read(spec, [b], multiscale=True)
        _render(sdata, spec, spec["viewports"]["overview"], shapes=_shapes_name(spec, b), color=spec["genes"][0],
                table=_table_name(b))


TASKS = {
    "image.read": image_read, "image.pyramid": image_pyramid, "image.query": image_query,
    "image.render": image_render,
    "bins.load": bins_load, "bins.query": bins_query, "bins.render": bins_render, "bins.aggregate": bins_aggregate,
    "bins.cluster": bins_cluster, "bins.de": bins_de,
    "cells.load": cells_load, "cells.pool": cells_pool, "cells.query": cells_query, "cells.render": cells_render,
    "cells.cluster": cells_cluster, "cells.de": cells_de,
    "e2e.first_view": e2e_first_view,
}

if __name__ == "__main__":
    main(TASKS, versions=("spatialdata", "spatialdata-io", "spatialdata-plot", "scanpy", "anndata", "dask"))
