from fastapi import APIRouter

from ...startup.datasets import AVAILABLE_DATASETS

router = APIRouter()


@router.get("/datasets")
def datasets():
    """Example datasets present on this server."""
    return [dict(id=name, title=name, bins_size=config["bins_size"],
                 has_reference=bool(config.get("reference_path")))
            for name, config in AVAILABLE_DATASETS.items()]
