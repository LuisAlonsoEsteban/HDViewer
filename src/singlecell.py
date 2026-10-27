"""Read uploaded or local single-cell references and barcode-aligned annotations."""

from io import BytesIO
from pathlib import Path
from tempfile import NamedTemporaryFile

import anndata as ad
import h5py
import pandas as pd


def _source(contents):
    return BytesIO(contents) if isinstance(contents, bytes) else Path(contents)


def read_reference(contents, annotations=None):
    """Read uploaded bytes or a local path without copying local count matrices."""
    with h5py.File(_source(contents), "r") as file:
        is_anndata = "obs" in file and "var" in file
    if is_anndata:
        reference = ad.read_h5ad(_source(contents))
    else:
        import scanpy as sc

        # Scanpy's public 10x reader takes a filename, including for legacy H5.
        if isinstance(contents, bytes):
            with NamedTemporaryFile(suffix=".h5") as file:
                file.write(contents)
                file.flush()
                reference = sc.read_10x_h5(file.name)
        else:
            reference = sc.read_10x_h5(contents)

    if reference.n_obs == 0 or reference.n_vars == 0:
        raise ValueError("The reference must contain cells and genes.")

    reference.var_names_make_unique()

    if annotations is not None:
        add_annotations(reference, annotations)

    return reference


def add_annotations(adata, annotations):
    """Align uploaded or local CSV/Parquet annotations to observations by barcode."""
    if isinstance(annotations, bytes):
        is_parquet = annotations.startswith(b"PAR1")
    else:
        with open(annotations, "rb") as file:
            is_parquet = file.read(4) == b"PAR1"

    if is_parquet:
        metadata = pd.read_parquet(_source(annotations))
        if "barcode" in metadata.columns:
            metadata = metadata.set_index("barcode")
        elif isinstance(metadata.index, pd.RangeIndex) and len(metadata.columns):
            metadata = metadata.set_index(metadata.columns[0])
    else:
        metadata = pd.read_csv(_source(annotations), index_col=0, dtype=str, keep_default_na=False)

    metadata.index = metadata.index.astype(str)
    if metadata.empty or not metadata.index.is_unique or not adata.obs_names.is_unique:
        raise ValueError("Annotation file must contain unique cell barcodes and annotation columns.")

    missing = adata.obs_names.difference(metadata.index)
    if len(missing):
        raise ValueError(f"Annotation file is missing {len(missing):,} reference cell barcodes.")
        
    metadata = metadata.reindex(adata.obs_names)
    for name in metadata.columns:
        adata.obs[name] = metadata[name]
