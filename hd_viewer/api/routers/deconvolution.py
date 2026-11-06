"""Deconvolution of a layer against the dataset's single-cell reference.

Methods are plugins run by the shared service (core.deconvolution); a run is a
background job (events `job` with name "deconvolution", separate from
clustering/DE since it may take hours) whose result becomes `deconv_<run>:<type>`
columns, shown with the layer's usual colour-by (PATCH style colour_by).
"""
import threading
import traceback

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field

from ...core import deconvolution as core
from . import current_session
from .layers import _loaded, _view
from .plugins import describe

router = APIRouter()


def _service():
    try:
        return core.default_service(), None
    except Exception as exc:
        return None, str(exc)


def _reference(session):
    loaded = _loaded(session)
    reference = getattr(loaded, "reference", None)
    columns = [] if reference is None else [str(column) for column in reference.obs.columns]
    preferred = getattr(loaded, "reference_label_column", "")
    label = preferred if preferred in columns else (columns[0] if columns else "")
    return reference, columns, label


def _unavailable(service, error, reference, columns):
    """Why no method can run on this dataset, or ''."""
    if service is None:
        return f"Deconvolution plugins unavailable: {error}"
    if reference is None:
        return ("No single-cell reference was provided for this dataset. Deconvolution is unavailable; "
                "the other layer tools remain available.")
    if not columns:
        return ("The reference contains counts but no .obs annotation columns. Upload a barcode-indexed "
                "annotation CSV or Parquet file alongside it in the data loading step.")
    return ""


def _runs(view):
    return [dict(id=run, method=record.get("method", ""), cell_types=list(record.get("cell_types", [])),
                 fraction_columns=list(record.get("fraction_columns", [])),
                 abundance_columns=list(record.get("abundance_columns", [])),
                 abundance_units=record.get("abundance_units", ""))
            for run, record in core.runs(view.sample).items()]


@router.get("/layers/{layer_id}/deconvolution")
def state(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    service, error = _service()
    reference, columns, label = _reference(session)
    with view.lock:
        runs = _runs(view)
    return dict(
        unavailable=_unavailable(service, error, reference, columns),
        reference=None if reference is None else dict(name=getattr(session.loaded, "reference_name", None)
                                                      or "loaded reference",
                                                      cells=reference.n_obs, genes=reference.n_vars),
        columns=columns, label_column=label, runner=service.runner.name if service else None,
        methods=[dict(describe(plugin), citation=plugin.citation, unavailable=service.unavailable_reason(name))
                 for name, plugin in service.plugins.items()] if service else [],
        runs=runs, job=view.deconvolution_job)


class RunRequest(BaseModel):
    method: str
    label_column: str
    params: dict = Field(default_factory=dict)


def _validate(session, body):
    service, error = _service()
    reference, columns, _ = _reference(session)
    reason = _unavailable(service, error, reference, columns)
    if reason:
        raise HTTPException(422, reason)
    if body.method not in service.plugins:
        raise HTTPException(422, f"Unknown deconvolution method: {body.method}")
    reason = service.unavailable_reason(body.method)
    if reason:
        raise HTTPException(422, reason)
    try:
        params, _ = service.validate_params(body.method, reference, body.label_column, body.params)
    except (ValueError, TypeError) as exc:
        raise HTTPException(422, str(exc))
    return service, reference, params


@router.post("/layers/{layer_id}/deconvolution/check")
def check(layer_id: str, body: RunRequest, session=Depends(current_session)):
    """Whether the method can run with these settings (the form's live hint)."""
    _view(session, layer_id)
    _validate(session, body)
    return dict(ok=True)


@router.post("/layers/{layer_id}/deconvolution", status_code=202)
def run(layer_id: str, body: RunRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    service, reference, params = _validate(session, body)
    loaded = session.loaded
    plugin = service.plugins[body.method]
    with view.lock:
        if view.deconvolution_job["state"] == "running":
            raise HTTPException(409, "A deconvolution is already running on this layer.")
        view.deconvolution_job = dict(name="deconvolution", state="running", message=f"Deconvolution ({plugin.title})")
    # As the Bokeh UI: cells deconvolve at the dataset's bin size, other layers at their own
    sample = view.sample
    dataset_size = getattr(getattr(loaded, "dataset", None), "bins_size", None)
    size = dataset_size if view.kind == "cells" else (getattr(sample, "bin_size", None) or dataset_size)
    namespace = core.cache_namespace(sample, view.kind, view.id, size)
    base_dir = getattr(loaded, "base_dir", None)
    work_dir = None if base_dir is None else base_dir / "deconvolution"
    target = core.snapshot_target(sample, body.method)

    def publish(**changes):
        view.deconvolution_job.update(changes)
        session.publish(dict(type="job", layer=layer_id, **view.deconvolution_job))

    def work():
        try:
            result = service.run(body.method, target, reference, body.label_column, params, work_dir=work_dir,
                                 cache_namespace=namespace, owner=session.sid,
                                 on_progress=lambda fraction, message: publish(message=str(message)))
            with view.lock:
                run_id = core.publish_result(sample, result, body.label_column,
                                             reference_name=getattr(loaded, "reference_name", "") or "")
            publish(state="done", run=run_id,
                    message=f"Deconvolved {len(result.fractions):,} observations into "
                            f"{len(result.fractions.columns)} cell types")
        except Exception as exc:
            traceback.print_exc()
            publish(state="failed", message=str(exc))

    threading.Thread(target=work, daemon=True, name=f"deconvolution-{layer_id}").start()
    return view.deconvolution_job
