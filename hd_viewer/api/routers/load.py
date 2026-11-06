"""Load a dataset (example or uploads) with the shared LoadPipeline.

The pipeline runs in a worker thread; progress, failures and readiness are
published as session events (see routers.events).
"""
import threading
import traceback

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field

from ...core.pipeline import LoadConfig, LoadPipeline, segmentation_service
from ...plugins.jobs import JOBS
from ...startup.datasets import AVAILABLE_DATASETS
from . import current_session
from .uploads import FIELDS, _field

router = APIRouter()


class SegmentationChoice(BaseModel):
    plugin: str
    params: dict = Field(default_factory=dict)


class LoadRequest(BaseModel):
    source: str                         # "example" | "upload"
    dataset_id: str | None = None
    bins_size: int | None = None        # uploads: bin size of the files
    pseudospots: bool = True
    segmentation: list[SegmentationChoice] = Field(default_factory=list)


def _config(session, body):
    options = dict(pseudospots=body.pseudospots, owner=session.sid)
    if body.segmentation:
        try:
            service = segmentation_service()
        except Exception as exc:
            raise HTTPException(422, f"Segmentation is unavailable: {exc}")
        segmentation = []
        for choice in body.segmentation:
            plugin = service.plugins.get(choice.plugin)
            if plugin is None:
                raise HTTPException(422, f"Unknown segmentation method: {choice.plugin}")
            try:
                segmentation.append((plugin.name, plugin.resolve_params(choice.params)))
            except ValueError as exc:
                raise HTTPException(422, f"{plugin.title}: {exc}")
        options["segmentation"] = segmentation

    if body.source == "example":
        if body.dataset_id not in AVAILABLE_DATASETS:
            raise HTTPException(422, "Please choose an example dataset.")
        return LoadConfig.from_example(AVAILABLE_DATASETS[body.dataset_id], **options)

    if body.source != "upload":
        raise HTTPException(422, "source must be 'example' or 'upload'.")
    fields = {key: _field(session, key) for key in FIELDS}
    missing = [key for key, (_, required) in FIELDS.items() if required and not fields[key].complete]
    if missing:
        raise HTTPException(422, f"Complete the required uploads first: {', '.join(missing)}.")
    reference, annotations = fields["reference"], fields["reference_annotations"]
    if not reference.complete and (annotations.complete or reference.status == reference.UPLOADING):
        raise HTTPException(422, "The single-cell reference has not finished uploading.")
    if reference.complete and annotations.status == annotations.UPLOADING:
        raise HTTPException(422, "The reference annotations have not finished uploading.")
    if body.bins_size not in (2, 8, 16):
        raise HTTPException(422, "Choose the bin size of the uploaded files (2, 8 or 16 µm).")
    return LoadConfig(img_path=str(fields["image"].path), adata_path=str(fields["expression"].path),
                      positions_path=str(fields["positions"].path), bins_size=body.bins_size,
                      base_dir=session.directory,
                      reference_path=str(reference.path) if reference.complete else None,
                      reference_annotations_path=str(annotations.path) if annotations.complete else None,
                      **options)


def cancel_load(session):
    if session.load.get("state") == "running":
        session.load["cancelled"] = True
        JOBS.cancel_owner(session.sid)


@router.get("/load")
def load_status(session=Depends(current_session)):
    return session.load


@router.post("/load", status_code=202)
def start_load(body: LoadRequest, session=Depends(current_session)):
    with session.lock:
        if session.load.get("state") == "running":
            raise HTTPException(409, "A dataset is already loading.")
        config = _config(session, body)
        session.load = dict(state="running", tasks={task: dict(pct=0, status="Waiting") for task in config.tasks},
                            error=None, cancelled=False)
        session.loaded, session.views, session.comparison = None, {}, {}

    def progress(task, pct, status):
        entry = session.load["tasks"].setdefault(task, dict(pct=0, status="Waiting"))
        if pct is not None:
            entry["pct"] = pct
        entry["status"] = status
        session.publish(dict(type="progress", task=task, pct=entry["pct"], status=status))

    def run():
        from .layers import open_workspace
        try:
            loaded = LoadPipeline(config, progress=progress).run()
            if session.load.get("cancelled"):
                raise RuntimeError("Loading was cancelled.")
            for task, message in loaded.failures.items():
                session.publish(dict(type="failed", task=task, message=message, fatal=False))
            with session.lock:
                session.loaded = loaded
                try:
                    open_workspace(session)
                except Exception:
                    session.loaded, session.views = None, {}
                    raise
            session.load["state"] = "ready"
            session.publish(dict(type="ready"))
        except Exception as exc:
            traceback.print_exc()
            session.load.update(state="failed", error=str(exc))
            session.publish(dict(type="failed", task=None, message=str(exc), fatal=True))

    threading.Thread(target=run, daemon=True, name=f"load-{session.sid[:6]}").start()
    return session.load
