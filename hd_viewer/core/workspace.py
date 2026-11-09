"""Samples (map layers' data) built from the outputs of the load pipeline."""
import pandas as pd

from ..settings.config import sample_kwargs
from ..startup.sample import Bins, Cells, Spots


def sequenced_geometries(gdf, adata, on="index"):
    """The rows of `gdf` with expression counts (their `on` key is in `adata`).

    Geometries without counts (bins outside the tissue, empty pseudospots, cells
    with no bins) are never shown or analysed, so layers drop them when they are
    built instead of indexing, fetching and sending them on every view.
    """
    if adata is None:
        return gdf
    keys = gdf.index if on == "index" and on not in gdf.columns else gdf[on]
    keep = pd.Index(adata.obs_names.astype(str)).get_indexer(pd.Index(keys).astype(str)) >= 0
    return gdf if keep.all() else gdf.iloc[keep]


def build_sample(sample_type, polygons, adata, b2s):
    """The Spots, Bins or Cells sample of a loaded dataset (sequenced geometries only)."""
    polygons = sequenced_geometries(polygons, adata, "index")
    spatial_transform = getattr(b2s, "spatial_transform", None)
    if sample_type == "spots":
        radius = (spatial_transform.pitch * b2s.converter.VisiumV2.SPOT_DIAMETER / (2 * b2s.bins_size)
                  if spatial_transform is not None else
                  b2s.converter.get_pixels(b2s.converter.VisiumV2.SPOT_DIAMETER) / 2 * b2s.tile_scale_xfactor)
        sample = Spots(polygons, adata, radius=radius, **sample_kwargs())
        sample.spatial_transform = spatial_transform
    elif sample_type == "bins":
        sample = Bins(
            polygons, adata, on='index', **sample_kwargs("bins"),
            pixels_per_bin=b2s.converter.visiumHD.pixels_per_bin,
            spatial_transform=spatial_transform,
            bin_size=b2s.bins_size,
            tile_scale_xfactor=b2s.tile_scale_xfactor,
            tile_scale_yfactor=b2s.tile_scale_yfactor,
        )
    elif sample_type == "cells":
        sample = Cells(polygons, adata, **sample_kwargs("cells"))
        sample.spatial_transform = spatial_transform
    else:
        raise ValueError(f"Unknown sample type: {sample_type}")
    return sample


def micrometres_per_unit(bins):
    """µm per full-resolution image pixel, from a Bins sample (None if unknown)."""
    if bins.pixels_per_bin and bins.bin_size > 0:
        return float(bins.bin_size) / bins.pixels_per_bin
    return None
