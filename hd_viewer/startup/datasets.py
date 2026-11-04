"""Example datasets offered on the loading screen.

Paths are relative to ``settings.DATA_DIR`` (``<repo>/data`` or ``$VHD_DATA_DIR``).
Datasets whose files are missing on this machine are not offered.
"""
from pathlib import Path

from ..settings import DATA_DIR

_DATASETS = {
    "Mouse brain (2um)": {
        "expression_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/binned_outputs/square_002um/filtered_feature_bc_matrix.h5"),
        "positions_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/binned_outputs/square_002um/spatial/tissue_positions.parquet"),
        "image_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_tissue_image.tif"),
        "bins_size": 2,
    },
    "Mouse brain (8um)": {
        "expression_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"),
        "positions_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/binned_outputs/square_008um/spatial/tissue_positions.parquet"),
        "image_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_tissue_image.tif"),
        "bins_size": 8,
    },
    "Mouse brain (16um)": {
        "expression_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/binned_outputs/square_016um/filtered_feature_bc_matrix.h5"),
        "positions_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/binned_outputs/square_016um/spatial/tissue_positions.parquet"),
        "image_path": str(DATA_DIR / "Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_tissue_image.tif"),
        "bins_size": 16,
    },
    "Human CRC": {
        "expression_path": str(DATA_DIR / "CRC_P1/VisiumHD/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"),
        "positions_path": str(DATA_DIR / "CRC_P1/VisiumHD/binned_outputs/square_008um/spatial/tissue_positions.parquet"),
        "image_path": str(DATA_DIR / "CRC_P1/VisiumHD/Visium_HD_Human_Colon_Cancer_P1_tissue_image.btf"),
        "bins_size": 8,
        "reference_path": str(DATA_DIR / "CRC_P1/Chromium/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5"),
        "reference_annotations_path": str(DATA_DIR / "CRC_P1/Chromium/SingleCell_MetaData.csv"),
        "reference_label_column": "Level1",
    },
}


def _is_available(config):
    return all(Path(value).exists() for key, value in config.items() if key.endswith("_path"))


AVAILABLE_DATASETS = {name: config for name, config in _DATASETS.items() if _is_available(config)}
