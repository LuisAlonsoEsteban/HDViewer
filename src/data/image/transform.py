"""One invertible coordinate system shared by images and spatial layers."""


from dataclasses import dataclass
import hashlib
import json

import numpy as np
import shapely


def fit_grid_basis(positions, *, max_points=4096, max_condition_number=1e10):
    sampled = positions.iloc[np.linspace(0, len(positions) - 1, min(max_points, len(positions)), dtype=int)]
    design = np.column_stack((sampled['array_col'], sampled['array_row'], np.ones(len(sampled))))
    pixels = sampled[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy(dtype=float)

    if not len(sampled) or not np.isfinite(design).all() or not np.isfinite(pixels).all():
        raise ValueError('Cannot align image: tissue positions must be finite.')

    basis, _, rank, _ = np.linalg.lstsq(design, pixels, rcond=None)

    if rank < 3 or np.linalg.cond(basis[:2]) > max_condition_number:
        raise ValueError('Cannot align image: the tissue grid is singular or underdetermined.')
    
    return basis


@dataclass(frozen=True)
class SpatialTransform:
    """Full-res pixels -> plot: uniform image scale, y reversal and padding.

    Images keep their original orientation and pixels. Only bin glyphs rotate.
    Matrices are stored as tuples so calibration remains immutable.
    """
    coefficients: tuple
    inverse_coefficients: tuple
    grid_coefficients: tuple
    output_size: int
    image_scale: float
    pitch: float
    fingerprint: str

    @property
    def forward(self):
        return np.array(self.coefficients).reshape(3, 3)

    @property
    def inverse(self):
        return np.array(self.inverse_coefficients).reshape(3, 3)

    @property
    def grid_forward(self):
        return np.array(self.grid_coefficients).reshape(3, 3)

    @property
    def bin_width(self):
        return float(np.linalg.norm(self.grid_forward[:2, 0]))

    @property
    def bin_height(self):
        return float(np.linalg.norm(self.grid_forward[:2, 1]))

    @property
    def angle(self):
        edge = self.grid_forward[:2, 0]
        return float(np.arctan2(edge[1], edge[0]))

    def rectangles(self, centres, factor=1):
        """The same rotated rectangle footprints used by Bokeh and native ROIs."""
        angle = self.angle
        rotation = np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])
        offsets = np.array([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]])
        offsets = (offsets * [self.bin_width * factor, self.bin_height * factor]) @ rotation.T
        return shapely.polygons(np.asarray(centres)[:, None, :] + offsets)

    @property
    def bounds(self):
        return (0., 0., float(self.output_size), float(self.output_size))

    @staticmethod
    def _points(points, matrix):
        points = np.asarray(points, dtype=float)
        return points @ matrix[:2, :2].T + matrix[:2, 2]

    def points(self, points, inverse=False):
        return self._points(points, self.inverse if inverse else self.forward)

    def image_points(self, points):
        return self.points(np.asarray(points) / self.image_scale)

    def array_points(self, points):
        return self._points(points, self.grid_forward)

    def geometry(self, geometry, inverse=False):
        return shapely.transform(geometry, lambda xy: self.points(xy, inverse=inverse))

    def image_geometry(self, geometry):
        return shapely.transform(geometry, self.image_points)

    @classmethod
    def fit(cls, positions, image_shape, image_scale=1., *, max_points=4096, max_condition_number=1e10):
        if not np.isfinite(image_scale) or image_scale <= 0:
            raise ValueError('Cannot align image: image scale must be positive.')

        for name in ('array_col', 'array_row', 'pxl_col_in_fullres', 'pxl_row_in_fullres'):
            if not np.isfinite(positions[name].to_numpy(dtype=float)).all():
                raise ValueError('Cannot display image: tissue positions must be finite.')

        basis = fit_grid_basis(positions, max_points=max_points, max_condition_number=max_condition_number)

        array_to_pixel = np.eye(3)
        array_to_pixel[:2, :] = basis.T

        forward = np.diag([image_scale, -image_scale, 1.])

        height, width = image_shape[:2]
        points = np.array([[0,0],[width,0],[width,-height],[0,-height]], dtype=float)
        centres = cls._points(positions[['pxl_col_in_fullres', 'pxl_row_in_fullres']], forward)

        edges = forward[:2, :2] @ basis[:2].T
        bin_width, bin_height = np.linalg.norm(edges, axis=0)
        angle = np.arctan2(edges[1,0], edges[0,0])

        half = np.array([abs(np.cos(angle))*bin_width + abs(np.sin(angle))*bin_height,
                         abs(np.sin(angle))*bin_width + abs(np.cos(angle))*bin_height]) / 2
        low = np.minimum(points.min(axis=0), centres.min(axis=0) - half)
        high = np.maximum(points.max(axis=0), centres.max(axis=0) + half)
        translation = -np.floor(low)
        forward[:2, 2] = translation
        array_to_plot = forward @ array_to_pixel

        side = max(1, int(np.ceil((high + translation).max())))
        inverse = np.linalg.inv(forward)
        
        pitch = float((bin_width + bin_height) / 2)
        
        spec = dict(version=2, forward=forward.tolist(), size=side, image_scale=image_scale,
                    image_shape=list(image_shape), resampling='none', exterior='transparent')
        fingerprint = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:20]

        return cls(
                tuple(forward.ravel()), 
                tuple(inverse.ravel()), 
                tuple(array_to_plot.ravel()),
                side, 
                float(image_scale), 
                pitch, 
                fingerprint
            )

    def warp_image(self, image):
        """Pad the original raster without rotation, stretching or resampling."""
        from PIL import Image

        original = Image.fromarray(image).convert('RGBA')
        canvas = Image.new('RGBA', (self.output_size, self.output_size), (0, 0, 0, 0))

        left = int(round(self.forward[0, 2]))
        top = int(round(self.output_size - self.forward[1, 2]))
        
        canvas.paste(original, (left, top))
        return np.asarray(canvas)
