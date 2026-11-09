"""Steps that turn a Visium HD dataset into map layers, without any UI.

Each step reports progress through a callback and returns plain objects; the
Bokeh loading screen and the HTTP API run the same steps (each with its own
threading and progress display). `LoadPipeline` runs them all for one load.
"""
import hashlib
import threading
from dataclasses import dataclass, field
from pathlib import Path

from .. import settings
from ..settings.config import bins_to_spots_kwargs, grid_fit_kwargs, image_kwargs, segmentation_kwargs
from ..startup.data.VisiumHD import VisiumHD


def _noop(*args, **kwargs):
    pass


# -- steps --------------------------------------------------------------------

def image_key(img_path):
    """Identifies the image file, so tiles and segmentations are cached per image."""
    path = Path(img_path).resolve()
    stat = path.stat()
    return hashlib.sha1(f"{path}|{stat.st_size}|{stat.st_mtime_ns}".encode()).hexdigest()[:12]


def fit_transform(dataset):
    from ..startup.data.image.transform import SpatialTransform
    if dataset.df_tissue_positions is None:
        dataset.read_positions()
    return SpatialTransform.fit(dataset.df_tissue_positions, dataset.image_shape, dataset.image_scale,
                                **grid_fit_kwargs())


def make_image_processor(dataset, transform, base_dir):
    """(ImageProcessor, single_image): small images are served as one picture."""
    from ..startup.data.image.image import ImageProcessor

    def load_image():
        if dataset.img is None:
            dataset.read_img()
        return dataset.img

    cache_key = f"padded-v2-{image_key(dataset.img_path)}-{transform.fingerprint}"
    processor = ImageProcessor(dataset.img, Path(base_dir) / 'tiles' / cache_key,
                               spatial_transform=transform, image_loader=load_image, **image_kwargs())
    return processor, processor.mem <= settings.TILE_MEM_THRESHOLD


def create_tiles(processor, single_image, progress=_noop):
    """True if tiles were generated, False if a cached set was reused."""
    def on_progress(done, total):
        progress(int(100 * done / total) if total else 0, "Creating padded image tiles")
    return processor.create_tiles(overwrite=False, progress_callback=on_progress, single_image=single_image)


def read_inputs(dataset, annotations_path=None, progress=_noop):
    """Expression matrix, bin positions and optional spatial annotations."""
    if dataset.adata is None or dataset.df_tissue_positions is None:
        progress(None, "Reading expression data")
        if dataset.adata is None:
            dataset.read_adata()
        if dataset.df_tissue_positions is None:
            dataset.read_positions()
    if annotations_path:
        from ..startup.data.singlecell import add_annotations
        progress(None, "Reading spatial annotations")
        add_annotations(dataset.adata, annotations_path)


def build_b2s(dataset, transform, processor, progress=_noop):
    """Bins2Spots: native bins, their geometry and the pseudospot converter.
    `progress(done, total)` follows the pseudospot geometry stage."""
    from ..startup.operations.bins2spots import Bins2Spots
    return Bins2Spots(
        dataset.adata, dataset.df_tissue_positions,
        bins_size=dataset.bins_size,
        pixels_per_bin=dataset.pixels_per_bin,
        spatial_transform=transform,
        tile_scale_xfactor=processor.scale_xfactor * dataset.image_scale,
        tile_scale_yfactor=processor.scale_yfactor * dataset.image_scale,
        progress_callback=lambda step, done, total: progress(done, total),
        **bins_to_spots_kwargs(),
    )


def pseudospots(b2s, progress=_noop):
    """(polygons, adata) of Visium v2-like 55 µm pseudospots."""
    progress(33, "Building spot geometry")
    polygons = b2s.spot_polygons
    progress(66, "Aggregating expression")
    return polygons, b2s.convert(normalize=False)


@dataclass
class SegmentationInput:
    image_path: Path
    image_key: str
    work_dir: Path
    um_per_px: float | None


def prepare_segmentation_image(service, dataset, base_dir):
    """Write the image once per dataset for every segmentation plugin run."""
    key = image_key(dataset.img_path)
    work_dir = Path(base_dir) / "segmentation"
    if dataset.img is None:
        dataset.read_img()
    # OpenCV (png/tif) reads BGR; plugins always receive RGB
    image_path = service.prepare_image(dataset.img, key, work_dir,
                                       bgr=dataset.img_format in ("png", "tif", "tiff"))
    # Plugins see the given image: µm per pixel of that image (not of the full-res one)
    um_per_px = (float(dataset.bins_size) / (dataset.pixels_per_bin * dataset.image_scale)
                 if dataset.pixels_per_bin else None)
    return SegmentationInput(image_path, key, work_dir, um_per_px)


def segment(service, plugin, segmentation_input, params, transform, b2s, owner=None, progress=_noop):
    """(cells, adata): cell outlines in map coordinates and their pooled expression."""
    from ..startup.operations.pooler import CellsPooler

    def on_progress(fraction, message):
        progress(5 + int(80 * fraction), message or "Segmenting")

    cells, out_dir = service.run(plugin.name, segmentation_input.image_path, segmentation_input.image_key,
                                 segmentation_input.work_dir, params=params,
                                 um_per_px=segmentation_input.um_per_px, on_progress=on_progress, owner=owner)
    print(f"[segmentation] {plugin.title}: {len(cells)} cells ({out_dir})")
    if cells.empty:
        raise RuntimeError(f"{plugin.title} found no cells")
    # Segment the original image; transform only the resulting outlines.
    cells.geometry = transform.image_geometry(cells.geometry.values)
    cells["index"] = [f"{plugin.name}_{i + 1}" for i in range(len(cells))]
    cells.index = cells["index"].values

    progress(85, f"Aggregating expression of {len(cells)} cells")
    adata = CellsPooler(b2s.bin_squares, b2s.adata_hd, cells, keep_cols=[]).run(predicate="intersects")
    return cells, adata


def read_deconvolution_reference(path, annotations_path=None, name=None):
    from ..startup.data.singlecell import read_reference
    try:
        return read_reference(path, annotations_path)
    except Exception as exc:
        raise ValueError(f"Cannot load single-cell reference {name or Path(path).name}: {exc}") from exc


_segmentation_service = None
_segmentation_lock = threading.Lock()


def segmentation_service():
    """One queue (max_parallel) and runner per process, shared by all sessions."""
    global _segmentation_service
    with _segmentation_lock:
        if _segmentation_service is None:
            from ..startup.data.image.segmentation import SegmentationService
            _segmentation_service = SegmentationService(**segmentation_kwargs())
        return _segmentation_service


# -- a whole load ---------------------------------------------------------------

@dataclass
class LoadConfig:
    """What to load. Paths are server paths (an example dataset or session uploads)."""
    img_path: str
    adata_path: str
    positions_path: str
    bins_size: int
    base_dir: Path
    pseudospots: bool = True
    segmentation: list = field(default_factory=list)       # [(plugin name, params)]
    spatial_annotations_path: str | None = None
    reference_path: str | None = None
    reference_annotations_path: str | None = None
    reference_label_column: str = ""
    owner: str | None = None

    @classmethod
    def from_example(cls, config, **options):
        return cls(img_path=config["image_path"], adata_path=config["expression_path"],
                   positions_path=config["positions_path"], bins_size=config["bins_size"],
                   base_dir=Path(config["image_path"]).resolve().parent,
                   spatial_annotations_path=config.get("spatial_annotations_path"),
                   reference_path=config.get("reference_path"),
                   reference_annotations_path=config.get("reference_annotations_path"),
                   reference_label_column=config.get("reference_label_column", ""), **options)

    @property
    def tasks(self):
        return ["tiles", "bins", *(["spots"] if self.pseudospots else []),
                *(f"cells:{name}" for name, _ in self.segmentation)]


@dataclass
class LoadedDataset:
    dataset: VisiumHD
    transform: object
    image_processor: object
    b2s: object
    samples: dict = field(default_factory=dict)        # layer id -> (sample type, Sample, title)
    reference: object = None
    reference_name: str | None = None
    reference_label_column: str = ""
    failures: dict = field(default_factory=dict)       # task -> message (non-fatal)
    base_dir: Path | None = None                       # dataset folder: deconvolution results cache here


class LoadPipeline:
    """Run a whole load in the calling thread (call it from a worker thread).

    `progress(task, pct, status)` is called from this and helper threads.
    Tiles, bins and pseudospots are fatal on failure; a failed segmentation
    method is reported and skipped.
    """

    def __init__(self, config, progress=_noop, segmentation=None):
        self.config = config
        self.progress = progress
        self._segmentation = segmentation

    def _task(self, task):
        return lambda pct=None, status="Processing": self.progress(task, pct, status)

    def run(self):
        from .workspace import build_sample

        config = self.config
        dataset = VisiumHD(img_path=config.img_path, adata_path=config.adata_path,
                           positions_path=config.positions_path, bins_size=config.bins_size)

        reference, reference_name = None, None
        if config.reference_path:
            reference_name = Path(config.reference_path).name
            reference = read_deconvolution_reference(config.reference_path, config.reference_annotations_path)

        tiles = self._task("tiles")
        tiles(None, "Checking image metadata and cached tiles")
        transform = fit_transform(dataset)
        processor, single_image = make_image_processor(dataset, transform, config.base_dir)
        generated = create_tiles(processor, single_image, tiles)
        tiles(100, "Complete" if generated else "Using cached tiles")

        bins = self._task("bins")
        read_inputs(dataset, config.spatial_annotations_path, bins)
        bins(25, "Preparing bins")
        spots = self._task("spots")
        b2s = build_b2s(dataset, transform, processor,
                        lambda done, total: config.pseudospots and spots(33 + int(33 * done / total) if total else 33,
                                                                         "Building spot geometry"))
        bins(50, "Building bin geometry")
        loaded = LoadedDataset(dataset, transform, processor, b2s, reference=reference,
                               reference_name=reference_name, reference_label_column=config.reference_label_column,
                               base_dir=Path(config.base_dir))
        loaded.samples["bins"] = ("bins", build_sample("bins", b2s.bin_squares, b2s.adata_hd, b2s), "Bins")
        bins(100, "Complete")

        if config.pseudospots:
            polygons, adata = pseudospots(b2s, spots)
            loaded.samples["spots"] = ("spots", build_sample("spots", polygons, adata, b2s), "Spots")
            spots(100, "Complete")

        if config.segmentation:
            service = self._segmentation or segmentation_service()
            segmentation_input = prepare_segmentation_image(service, dataset, config.base_dir)
            for name, params in config.segmentation:
                task = self._task(f"cells:{name}")
                plugin = service.plugins[name]
                try:
                    cells, adata = segment(service, plugin, segmentation_input, params, transform, b2s,
                                           owner=config.owner, progress=task)
                except Exception as exc:
                    loaded.failures[f"cells:{name}"] = str(exc)
                    task(None, "Failed")
                    continue
                loaded.samples[f"cells_{name}"] = ("cells", build_sample("cells", cells, adata, b2s),
                                                   f"Cells · {plugin.title}")
                task(100, "Complete")
        return loaded
