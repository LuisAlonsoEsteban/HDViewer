"""Per-layer analyses: clustering (a background job) and custom groups.

Results recolour the layer; the browser is told with a `layer_stale` event (or
the response) and re-reads the layer's data (GET state + viewport).
"""
import threading
import traceback

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import Response
from pydantic import BaseModel, Field

from ... import settings
from ...core.clustering import SUPPORTED_ALGORITHMS
from ...core.groups import GROUPS_COL, OPERATORS
from ...settings.config import clustering_kwargs
from ...startup.utils.progress import reporting
from .. import frames
from . import current_session
from .layers import _view

router = APIRouter()
FRAME = "application/octet-stream"


def _frame(view, ops):
    return Response(frames.encode(*view.frame_parts(ops)), media_type=FRAME)


# -- clustering ---------------------------------------------------------------

class ClusteringRequest(BaseModel):
    algorithm: str = settings.CLUSTER_ALGORITHM
    resolution: float = settings.CLUSTER_SLIDER_DEFAULT
    n_clusters: int | None = None
    eps: float | None = None
    min_samples: int | None = None
    min_cluster_size: int | None = None
    covariance_type: str | None = None
    linkage: str | None = None


def _clusters(view):
    clustering = view.clustering
    return dict(column=clustering.cluster_col, clusters=clustering.cluster_summary() if clustering.is_clustered else [],
                summary=clustering._result_summary() if clustering.is_clustered else "",
                shown=view.colour_by.get("kind") == "clusters", job=view.job)


@router.get("/clustering/options")
def clustering_options():
    defaults = clustering_kwargs()
    return dict(algorithms=list(SUPPORTED_ALGORITHMS), defaults=defaults,
                resolution=dict(min=settings.CLUSTER_SLIDER_MIN, max=settings.CLUSTER_SLIDER_MAX,
                                step=settings.CLUSTER_SLIDER_STEP, value=settings.CLUSTER_SLIDER_DEFAULT))


@router.get("/layers/{layer_id}/clusters")
def clusters(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return _clusters(view)


@router.post("/layers/{layer_id}/clustering", status_code=202)
def run_clustering(layer_id: str, body: ClusteringRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    if body.algorithm not in SUPPORTED_ALGORITHMS:
        raise HTTPException(422, f"Unknown clustering algorithm: {body.algorithm}")
    options = {key: value for key, value in body.model_dump().items()
               if value is not None and key not in ("resolution",)}
    with view.lock:
        if view.job["state"] == "running":
            raise HTTPException(409, "An analysis is already running on this layer.")
        view.job = dict(name="clustering", state="running", message="Starting…")
        level = view.rows.effective_resolution_level

    def publish(**changes):
        view.job.update(changes)
        session.publish(dict(type="job", layer=layer_id, **view.job))

    def work():
        clustering = view.clustering
        try:
            with reporting(lambda message: publish(message=str(message))):
                column = clustering.compute(body.resolution, resolution_level=level, **options)
            with view.lock:
                view.drop_comparisons(clusters=True)
                clustering._build_cluster_state()
                view.colour_by_categories(column, clustering._color_map, "clusters")
            session.publish(dict(type="layer_stale", layer=layer_id))
            publish(message=f"{clustering._result_summary()} · marker genes")
            with reporting(lambda message: publish(message=str(message))):
                clustering.compute_markers()
            publish(state="done", message=f"{clustering._result_summary()}, marker genes ready")
        except Exception as exc:
            traceback.print_exc()
            publish(state="failed", message=str(exc))

    threading.Thread(target=work, daemon=True, name=f"clustering-{layer_id}").start()
    return view.job


class ClusterDisplay(BaseModel):
    visible: list[str] | None = None        # labels to show
    colours: dict[str, str] | None = None   # label -> colour
    show: bool | None = None                # colour the layer by the clusters again


@router.patch("/layers/{layer_id}/clusters")
def display_clusters(layer_id: str, body: ClusterDisplay, session=Depends(current_session)):
    view = _view(session, layer_id)
    clustering = view.clustering
    with view.lock:
        if not clustering.is_clustered:
            raise HTTPException(409, "Run clustering first.")
        ops = []
        if body.show:
            ops += view.colour_by_categories(clustering.cluster_col, clustering._color_map, "clusters")
        if body.colours:
            clustering._color_map.update({label: colour for label, colour in body.colours.items()
                                          if label in clustering._cluster_labels})
            if view.colour_by.get("kind") == "clusters":
                ops += view.recolour_categories(clustering._color_map)
        if body.visible is not None:
            clustering._active_clusters = set(body.visible) & set(clustering._cluster_labels)
        return _frame(view, ops)


# -- groups -------------------------------------------------------------------

class Condition(BaseModel):
    column: str
    op: str
    value: str
    logic: str = "AND"


class GroupRequest(BaseModel):
    name: str = ""
    conditions: list[Condition] = Field(default_factory=list)


class GroupDisplay(BaseModel):
    visible: bool | None = None
    colour: str | None = None


def _groups(view):
    return dict(groups=view.groups.group_summary(), operators=OPERATORS,
                shown=view.colour_by.get("kind") == "groups")


def _show_groups(view):
    groups = view.groups
    if groups._groups:
        groups.refresh_active_labels()
        return view.colour_by_categories(GROUPS_COL, groups._color_map, "groups")
    obs = view.sample.adata.obs
    if GROUPS_COL in obs:
        obs.drop(columns=[GROUPS_COL], inplace=True)
    return view.colour_by_column(view.sample.adata.var_names[0])


@router.get("/layers/{layer_id}/groups")
def groups(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return _groups(view)


@router.post("/layers/{layer_id}/groups/check")
def check_condition(layer_id: str, body: Condition, session=Depends(current_session)):
    view = _view(session, layer_id)
    try:
        return dict(hint=view.groups.check_condition(body.column, body.op, body.value))
    except ValueError as exc:
        raise HTTPException(422, str(exc))


@router.post("/layers/{layer_id}/groups")
def create_group(layer_id: str, body: GroupRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    conditions = [(c.column, c.op, c.value, c.logic if c.logic in ("AND", "OR") else "AND") for c in body.conditions]
    try:
        prepared = view.groups.prepare_condition_group(body.name, conditions)   # heavy: outside the lock
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    with view.lock:
        if prepared[0] in view.groups._groups:
            raise HTTPException(422, f"Group '{prepared[0]}' already exists.")
        view.groups.commit_condition_group(*prepared)
        view.drop_comparisons(groups_changed=True)
        _show_groups(view)
        return _groups(view)


@router.patch("/layers/{layer_id}/groups/{name}")
def display_group(layer_id: str, name: str, body: GroupDisplay, session=Depends(current_session)):
    view = _view(session, layer_id)
    groups = view.groups
    with view.lock:
        if name not in groups._groups:
            raise HTTPException(404, f"Unknown group: {name}")
        ops = []
        if body.colour:
            groups._color_map[name] = body.colour
            if view.colour_by.get("kind") == "groups":
                ops += view.recolour_categories(groups._color_map)
        if body.visible is not None:
            (groups._active_groups.add if body.visible else groups._active_groups.discard)(name)
            groups.refresh_active_labels()
            if view.colour_by.get("kind") == "groups":
                view.rows.update_colname(GROUPS_COL, GROUPS_COL, attrs=("color",))   # new labels
                ops += view.rows.sink.take()
        return _frame(view, ops)


@router.delete("/layers/{layer_id}/groups/{name}")
def delete_group(layer_id: str, name: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    groups = view.groups
    with view.lock:
        if name not in groups._groups:
            raise HTTPException(404, f"Unknown group: {name}")
        view.selected -= groups.members(name)
        if groups.forget_group(name):
            groups.store_groups()
        view.drop_comparisons(groups=[name])
        _show_groups(view)
        return _groups(view)


@router.post("/layers/{layer_id}/groups/show")
def show_groups(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        if not view.groups._groups:
            raise HTTPException(409, "Create a group first.")
        return _frame(view, _show_groups(view))
