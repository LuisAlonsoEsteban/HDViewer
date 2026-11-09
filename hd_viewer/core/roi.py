"""Regions of interest (ROIs), without any UI.

`ROIGeometry` builds ROI shapes (box, lasso, radius around a point, a barcode or a
group's bins), computes their native-observation membership, resolves barcodes
and prepares outlines for display. The HTTP API holds one per layer view.
"""
import threading
import uuid
from dataclasses import dataclass, field

import numpy as np
import shapely
from shapely.affinity import scale
from shapely.geometry import Point, Polygon, box

from .. import settings
from ..startup.utils.progress import report


@dataclass(eq=False)
class ROI:
    name: str
    geometry: object
    kind: str
    color: str
    id: str = field(default_factory=lambda: uuid.uuid4().hex)
    mask: np.ndarray | None = None
    visible: bool = True
    error: str | None = None
    seed_group: str | None = None
    seed_count: int | None = None
    excluded_seed_mask: np.ndarray | None = None
    seed_barcode: str | None = None

    @property
    def count(self):
        return None if self.mask is None else int(self.mask.sum())


@dataclass
class _ROIOverlay:
    """Bounded display data; never used for membership or saved geometry."""
    geometry: object
    xs: object = None
    ys: object = None
    raster: object = None
    bounds: object = None
    images: dict = field(default_factory=dict)

    def image(self, color):
        if color not in self.images:
            from PIL.ImageColor import getrgb
            rgba = np.zeros((*self.raster.shape, 4), dtype=np.uint8)
            rgba[self.raster > 0, :3] = getrgb(color)[:3]
            rgba[self.raster == 1, 3] = 15  # Same light fill as vector ROIs.
            rgba[self.raster == 2, 3] = 255
            self.images[color] = rgba.view(np.uint32).reshape(self.raster.shape)
        return self.images[color]


class ROIGeometry:
    def __init__(self, sample=None):
        self._sample = sample
        self.rois = []
        # Object whose `.data` holds px_x, px_y, um_per_px (calibration of radii)
        self.coordinate_scale = None
        self._counter = 1
        self._barcode_cache = None
        self._overlay_cache = {}
        self._group_union_cache = None
        self._group_union_lock = threading.Lock()

    @property
    def sample(self):
        return self._sample

    @staticmethod
    def validate_geometry(geometry):
        if geometry.is_empty or geometry.geom_type not in ('Polygon', 'MultiPolygon') or not geometry.is_valid or geometry.area <= 0:
            raise ValueError("Draw a nonempty, non-self-intersecting region with positive area.")
        if not np.isfinite(shapely.get_coordinates(geometry)).all():
            raise ValueError("ROI coordinates must be finite.")
        return geometry

    @classmethod
    def box_geometry(cls, x0, y0, x1, y1):
        if not np.isfinite([x0, y0, x1, y1]).all() or x0 == x1 or y0 == y1:
            raise ValueError("A bounding box needs finite coordinates and positive width and height.")
        return cls.validate_geometry(box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))

    @classmethod
    def lasso_geometry(cls, xs, ys):
        if len(xs) != len(ys) or len(xs) < 3 or not np.isfinite([xs, ys]).all():
            raise ValueError("A lasso needs at least three finite vertices.")
        return cls.validate_geometry(Polygon(zip(xs, ys)))

    @classmethod
    def radius_geometry(cls, x, y, radius, units_x=1., units_y=1.):
        if not np.isfinite([x, y, radius, units_x, units_y]).all() or min(radius, units_x, units_y) <= 0:
            raise ValueError("Seed coordinates must be finite and the radius and calibration must be positive.")
        circle = Point(x, y).buffer(radius, quad_segs=settings.ROI_BUFFER_QUAD_SEGS)
        return cls.validate_geometry(scale(circle, xfact=1 / units_x, yfact=1 / units_y, origin=(x, y)))

    def seed_radius_geometry(self, x, y, radius, units_x=1., units_y=1.):
        transform = getattr(self.sample, 'spatial_transform', None)
        if transform is None:
            return self.radius_geometry(x, y, radius, units_x, units_y)
        seed = transform.points([x, y], inverse=True)
        circle = self.radius_geometry(*seed, radius, units_x, units_y)
        return self.validate_geometry(transform.geometry(circle))

    def radius_scale(self):
        data = None if self.coordinate_scale is None else self.coordinate_scale.data
        fx = (data['px_x'][0] if data else None) or 1 / getattr(self.sample, 'tile_scale_xfactor', 1.)
        fy = (data['px_y'][0] if data else None) or 1 / getattr(self.sample, 'tile_scale_yfactor', 1.)
        um = data['um_per_px'][0] if data else None
        if getattr(self.sample, 'spatial_transform', None) is not None:
            fx = fy = 1.
        return ((fx * um, fy * um), 'µm') if um is not None and um > 0 else ((fx, fy), 'full-res px')

    def compute_membership(self, geometry):
        """Intersect full native geometries, then align to sequenced observations."""
        with self.sample.fetch_lock:
            gdf = self.sample.gdf
            candidates = gdf.sindex.query(geometry, predicate='intersects')
            names = gdf.iloc[candidates][self.sample.on].astype(str)
            return np.asarray(self.sample.adata.obs_names.isin(names), dtype=bool)

    def compute_group_radius(self, group, radius, factors, mask=None, cancelled=None, exclude_group=False):
        """Union native-bin radius neighbourhoods, optionally subtracting the group."""
        if not hasattr(self.sample, 'AGGREGATION_FACTORS'):
            raise ValueError("Group-bin seeds are available on bins layers.")
        self.radius_geometry(0., 0., radius, *factors)  # Validate radius and calibration.
        with self.sample.fetch_lock:
            key = f'{settings.GROUP_MASK_PREFIX}{group}'
            if key not in self.sample.adata.obs:
                raise ValueError("The seed group no longer exists.")
            if mask is None:
                mask = self.sample.adata.obs[key].fillna(False).to_numpy(dtype=bool)
            names = self.sample.adata.obs_names[np.asarray(mask, dtype=bool)]
            if not len(names):
                raise ValueError("The seed group contains no sequenced bins.")
            gdf = self.sample.gdf
            rows = gdf.loc[gdf[self.sample.on].astype(str).isin(names)].drop_duplicates(self.sample.on)
            group_geometries = rows.geometry.to_numpy().copy()
            centers = shapely.centroid(group_geometries)
            coordinates = np.column_stack((shapely.get_x(centers), shapely.get_y(centers)))
            if len(rows) != len(names) or not np.isfinite(coordinates).all():
                raise ValueError("Every grouped bin must have valid geometry for its seed.")
        transform = getattr(self.sample, "spatial_transform", None)
        if transform is not None:
            coordinates = transform.points(coordinates, inverse=True)
        coordinates *= np.asarray(factors)
        levels = []
        # Bound temporary geometry memory; nearby circles merge within each chunk.
        for start in range(0, len(coordinates), settings.ROI_SEED_BATCH_SIZE):
            if cancelled is not None and cancelled():
                raise InterruptedError("Group preview cancelled.")
            points = shapely.points(coordinates[start:start + settings.ROI_SEED_BATCH_SIZE])
            region = shapely.union_all(shapely.buffer(points, radius, quad_segs=settings.ROI_BUFFER_QUAD_SEGS))
            level = 0
            # Merge similarly sized batches instead of repeatedly unioning an
            # ever-growing region, especially costly for disjoint group bins.
            while level < len(levels) and levels[level] is not None:
                if cancelled is not None and cancelled():
                    raise InterruptedError("Group preview cancelled.")
                region = shapely.union_all([levels[level], region])
                levels[level] = None
                level += 1
            if level == len(levels):
                levels.append(region)
            else:
                levels[level] = region
            report(f"Buffered {min(start + settings.ROI_SEED_BATCH_SIZE, len(coordinates)):,} / {len(coordinates):,} group-bin seeds")

        region = shapely.union_all([region for region in levels if region is not None])
        region = scale(region, xfact=1 / factors[0], yfact=1 / factors[1], origin=(0, 0))
        if transform is not None:
            region = transform.geometry(region)

        if exclude_group:
            if cancelled is not None and cancelled():
                raise InterruptedError("Group preview cancelled.")

            with self._group_union_lock:
                cached = self._group_union_cache
            if (cached is not None and cached[0] is gdf and cached[1] is gdf.geometry.values
                    and cached[2].equals(names)):
                report("Reusing source-group footprint")
                footprint = cached[3]
            else:
                report("Combining source-group bin geometries")
                parts = []
                for start in range(0, len(group_geometries), settings.ROI_SEED_BATCH_SIZE):
                    self._check_cancelled(cancelled)
                    parts.append(shapely.union_all(group_geometries[start:start + settings.ROI_SEED_BATCH_SIZE]))
                
                while len(parts) > 1:
                    merged = []
                    for start in range(0, len(parts), 2):
                        self._check_cancelled(cancelled)
                        merged.append(shapely.union_all(parts[start:start + 2]))
                    parts = merged

                footprint = parts[0]
                
                self._check_cancelled(cancelled)
                with self._group_union_lock:
                    self._group_union_cache = (gdf, gdf.geometry.values, names.copy(), footprint)

            self._check_cancelled(cancelled)
            report("Subtracting source-group footprint")
            region = region.difference(footprint)
            if region.is_empty or region.area <= 0:
                raise ValueError("Subtracting the group leaves no ROI. Increase the radius or choose Full ROI.")

        if cancelled is not None and cancelled():
            raise InterruptedError("Group preview cancelled.")

        return self.validate_geometry(region), len(names)

    def _barcode_index(self):
        """Sorted identifiers and native positions, built only when first used."""
        gdf = self.sample.gdf
        identifiers = gdf[self.sample.on]
        cached = self._barcode_cache
        if cached is None or cached[0] is not gdf or cached[1] is not identifiers._values:
            valid = identifiers.notna().to_numpy()
            names = identifiers.astype(str).to_numpy()[valid]
            order = np.argsort(names, kind='stable')
            cached = self._barcode_cache = (gdf, identifiers._values, names[order], np.flatnonzero(valid)[order])
        return cached[2], cached[3]

    def _barcode_suggestions(self, prefix):
        prefix = prefix.strip()
        if not prefix:
            return []
        with self.sample.fetch_lock:
            names, _ = self._barcode_index()
            start = names.searchsorted(prefix, side='left')
            candidates = names[start:start + settings.ROI_BARCODE_MAX_COMPLETIONS]
            return list(dict.fromkeys(name for name in candidates if name.startswith(prefix)))

    def barcode_seed_center(self, barcode):
        """Resolve an exact native identifier independently of displayed rows."""
        barcode = barcode.strip()
        if not barcode:
            raise ValueError('Enter a native barcode or bin index.')
        with self.sample.fetch_lock:
            names, rows = self._barcode_index()
            start, end = names.searchsorted(barcode, side='left'), names.searchsorted(barcode, side='right')
            if start == end:
                raise ValueError(f"Barcode '{barcode}' was not found in this layer. Choose an exact barcode from the suggestions.")
            if end - start != 1:
                raise ValueError(f"Barcode '{barcode}' identifies more than one native bin.")
            centre = shapely.centroid(self.sample.gdf.geometry.iloc[rows[start]])
            x, y = float(shapely.get_x(centre)), float(shapely.get_y(centre))
            if not np.isfinite([x, y]).all():
                raise ValueError(f"Barcode '{barcode}' has no valid seed centre.")
            return x, y

    @staticmethod
    def _check_cancelled(cancelled):
        if cancelled is not None and cancelled():
            raise InterruptedError("Group preview cancelled.")

    @classmethod
    def _prepare_overlay(cls, geometry, cancelled=None):
        """Prepare large shapes off the document thread without exporting their vertices."""
        cls._check_cancelled(cancelled)
        if shapely.get_num_coordinates(geometry) <= settings.ROI_OVERLAY_MAX_VERTICES:
            parts = [geometry] if geometry.geom_type == 'Polygon' else list(geometry.geoms)
            xs, ys = [], []
            for poly in parts:
                cls._check_cancelled(cancelled)
                rings = [poly.exterior, *poly.interiors]
                coords = [np.asarray(ring.coords) for ring in rings]
                xs.append([xy[:, 0].tolist() for xy in coords])
                ys.append([xy[:, 1].tolist() for xy in coords])
            return _ROIOverlay(geometry, xs=xs, ys=ys)

        from scipy.ndimage import binary_erosion
        report("Rasterizing complex ROI preview; membership keeps the exact geometry")
        x0, y0, x1, y1 = geometry.bounds
        dx, dy = x1 - x0, y1 - y0
        size = max(4, int(settings.ROI_OVERLAY_MAX_SIZE))
        width = min(size, max(4, int(np.ceil(size * dx / max(dx, dy)))))
        height = min(size, max(4, int(np.ceil(size * dy / max(dx, dy)))))
        # A transparent margin keeps the outermost outline visible.
        x0, x1 = x0 - dx / width, x1 + dx / width
        y0, y1 = y0 - dy / height, y1 + dy / height
        xs = x0 + (np.arange(width) + .5) * (x1 - x0) / width
        ys = y0 + (np.arange(height) + .5) * (y1 - y0) / height
        mask = np.empty((height, width), dtype=bool)
        shapely.prepare(geometry)
        # Row batches keep stale rasterization cancellable and memory bounded.
        for start in range(0, height, 32):
            cls._check_cancelled(cancelled)
            mask[start:start + 32] = shapely.intersects_xy(geometry, xs[None, :], ys[start:start + 32, None])
        raster = mask.astype(np.uint8)
        raster[mask & ~binary_erosion(mask)] = 2
        return _ROIOverlay(geometry, raster=raster, bounds=(x0, y0, x1, y1))

    def membership(self, roi):
        """Native observations inside `roi` (minus its excluded seed bins)."""
        mask = self.compute_membership(roi.geometry)
        # Intersects includes boundary contact. Remove the captured source bins
        # explicitly even when they touch the edge of a subtracted polygon.
        return mask if roi.excluded_seed_mask is None else mask & ~roi.excluded_seed_mask

    def next_name(self):
        while any(roi.name == f"ROI {self._counter}" for roi in self.rois):
            self._counter += 1
        return f"ROI {self._counter}"

    def next_colour(self):
        return settings.GROUP_PALETTE[(self._counter - 1) % len(settings.GROUP_PALETTE)]
