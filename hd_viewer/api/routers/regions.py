"""Regions of interest of a layer: drawing, membership, groups and selection.

Shapes are in map coordinates (those of the map figure). Membership is exact
(native geometries); outlines sent for display are simplified when large.
"""
import numpy as np
import shapely
from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel

from ... import settings
from ...core.roi import ROI
from . import current_session
from .analysis import _groups, _show_groups
from .layers import _view

router = APIRouter()


class Shape(BaseModel):
    kind: str                                  # box | lasso | point | barcode | group
    name: str = ""
    box: list[float] | None = None             # x0, y0, x1, y1
    xs: list[float] | None = None              # lasso vertices
    ys: list[float] | None = None
    centre: list[float] | None = None          # point seed (map coordinates)
    barcode: str | None = None
    group: str | None = None
    radius: float | None = None                # in the calibrated unit (µm when known)
    exclude_group: bool = False


def _outline(rois, geometry):
    """Rings for a Bokeh multi_polygons row, simplified to a bounded size."""
    limit = settings.ROI_OVERLAY_MAX_VERTICES
    tolerance = 0.0
    x0, y0, x1, y1 = geometry.bounds
    while shapely.get_num_coordinates(geometry) > limit:
        tolerance = tolerance * 2 if tolerance else max(x1 - x0, y1 - y0) / 2000
        geometry = geometry.simplify(tolerance, preserve_topology=True)
    overlay = rois._prepare_overlay(geometry)
    return overlay.xs, overlay.ys


def _describe(view, roi):
    xs, ys = _outline(view.rois, roi.geometry)
    return dict(id=roi.id, name=roi.name, kind=roi.kind, colour=roi.color, count=roi.count,
                visible=roi.visible, error=roi.error, group=view.groups.group_for_roi(roi.id),
                seed_group=roi.seed_group, seed_count=roi.seed_count, seed_barcode=roi.seed_barcode,
                group_excluded=roi.excluded_seed_mask is not None, xs=xs, ys=ys)


def _list(view):
    factors, unit = view.rois.radius_scale()
    return dict(rois=[_describe(view, roi) for roi in view.rois.rois], unit=unit,
                seed_groups=list(view.groups._groups) if hasattr(view.sample, "AGGREGATION_FACTORS") else [])


def _group_geometry(view, shape, factors, preview):
    """Union of the group bins' radius neighbourhoods (slow on large groups).

    A new preview cancels the one still computing; saving reuses the last
    preview when the group's members have not changed since."""
    key = f"{settings.GROUP_MASK_PREFIX}{shape.group}"
    state = view.roi_preview
    with view.lock:
        if key not in view.sample.adata.obs:
            raise ValueError("The seed group no longer exists.")
        mask = view.sample.adata.obs[key].fillna(False).to_numpy(dtype=bool, copy=True)
        request = (shape.group, shape.radius, shape.exclude_group, tuple(factors))
        cached = state["result"]
        if state["key"] == request and np.array_equal(cached[2], mask):
            return cached
        if preview:
            state["generation"] += 1
        generation = state["generation"]
    cancelled = (lambda: state["generation"] != generation) if preview else None
    try:
        geometry, count = view.rois.compute_group_radius(shape.group, shape.radius, factors, mask=mask,
                                                         cancelled=cancelled, exclude_group=shape.exclude_group)
    except InterruptedError:
        raise HTTPException(409, "Replaced by a newer preview.")
    result = (geometry, (shape.group, count, mask if shape.exclude_group else None), mask)
    with view.lock:
        state.update(key=request, result=result)
    return result


def _geometry(view, shape, preview=False):
    """(geometry, seed info or None) of a shape; seed info is (group, count, excluded mask)."""
    rois = view.rois
    try:
        if shape.kind == "box" and shape.box and len(shape.box) == 4:
            return rois.box_geometry(*shape.box), None
        if shape.kind == "lasso" and shape.xs is not None:
            return rois.lasso_geometry(shape.xs, shape.ys or []), None
        if shape.radius is None:
            raise ValueError("Choose a radius.")
        factors, _ = rois.radius_scale()
        if shape.kind == "point" and shape.centre and len(shape.centre) == 2:
            return rois.seed_radius_geometry(*shape.centre, shape.radius, *factors), None
        if shape.kind == "barcode":
            centre = rois.barcode_seed_center(shape.barcode or "")
            return rois.seed_radius_geometry(*centre, shape.radius, *factors), None
        if shape.kind == "group" and shape.group:
            geometry, seed, _ = _group_geometry(view, shape, factors, preview)
            return geometry, seed
    except (ValueError, KeyError, TypeError) as exc:
        raise HTTPException(422, str(exc))
    raise HTTPException(422, "Unknown or incomplete ROI shape.")


@router.get("/layers/{layer_id}/rois")
def list_rois(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return _list(view)


@router.post("/layers/{layer_id}/rois/preview")
def preview(layer_id: str, shape: Shape, session=Depends(current_session)):
    view = _view(session, layer_id)
    geometry, seed = _geometry(view, shape, preview=True)
    xs, ys = _outline(view.rois, geometry)
    return dict(xs=xs, ys=ys, bounds=list(geometry.bounds), seeds=None if seed is None else seed[1])


@router.post("/layers/{layer_id}/rois")
def create_roi(layer_id: str, shape: Shape, session=Depends(current_session)):
    view = _view(session, layer_id)
    geometry, seed = _geometry(view, shape)
    rois = view.rois
    with view.lock:
        name = shape.name.strip() or rois.next_name()
        if any(roi.name == name for roi in rois.rois):
            raise HTTPException(422, "An ROI with that name already exists in this layer.")
        roi = ROI(name, geometry, shape.kind, rois.next_colour())
        roi.seed_barcode = shape.barcode if shape.kind == "barcode" else None
        if seed is not None:
            roi.seed_group, roi.seed_count, roi.excluded_seed_mask = seed
        rois._counter += 1
    roi.mask = rois.membership(roi)                     # exact, possibly slow: outside the lock
    with view.lock:
        rois.rois.append(roi)
        return _list(view)


class ROIDisplay(BaseModel):
    visible: bool | None = None


def _roi(view, roi_id):
    roi = next((roi for roi in view.rois.rois if roi.id == roi_id), None)
    if roi is None:
        raise HTTPException(404, "Unknown ROI.")
    return roi


@router.patch("/layers/{layer_id}/rois/{roi_id}")
def display_roi(layer_id: str, roi_id: str, body: ROIDisplay, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        roi = _roi(view, roi_id)
        if body.visible is not None:
            roi.visible = body.visible
        return _list(view)


@router.delete("/layers/{layer_id}/rois/{roi_id}")
def delete_roi(layer_id: str, roi_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        roi = _roi(view, roi_id)
        linked = view.groups.group_for_roi(roi.id)
        if linked:
            view.selected -= view.groups.members(linked)
            if view.groups.forget_group(linked):
                view.groups.store_groups()
            view.drop_comparisons(groups=[linked])
            if view.colour_by.get("kind") == "groups" or not view.groups._groups:
                _show_groups(view)
        if roi.mask is not None:
            view.selected -= set(view.sample.adata.obs_names[roi.mask])
        view.rois.rois = [item for item in view.rois.rois if item is not roi]
        return dict(**_list(view), groups=_groups(view)["groups"], recoloured=bool(linked))


@router.post("/layers/{layer_id}/rois/{roi_id}/group")
def roi_group(layer_id: str, roi_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        roi = _roi(view, roi_id)
        if roi.mask is None or not roi.mask.any():
            raise HTTPException(422, "This ROI contains no sequenced observations.")
        name = view.groups.add_fixed_group(f"ROI: {roi.name}", roi.mask, origin_roi_id=roi.id)
        view.drop_comparisons(groups_changed=True)
        _show_groups(view)
        return dict(**_list(view), group=name, groups=_groups(view)["groups"])


@router.post("/layers/{layer_id}/rois/{roi_id}/select")
def select_roi(layer_id: str, roi_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        roi = _roi(view, roi_id)
        if roi.mask is not None:
            view.selected = set(view.sample.adata.obs_names[roi.mask])
        return dict(row_version=view.rows.row_version, rows=view.selected_positions(), selected=len(view.selected))


@router.post("/layers/{layer_id}/rois/{roi_id}/unselect")
def unselect_roi(layer_id: str, roi_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        roi = _roi(view, roi_id)
        if roi.mask is not None:
            view.selected -= set(view.sample.adata.obs_names[roi.mask])
        return dict(row_version=view.rows.row_version, rows=view.selected_positions(), selected=len(view.selected))


@router.post("/layers/{layer_id}/rois/{roi_id}/focus")
def focus_roi(layer_id: str, roi_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        x0, y0, x1, y1 = _roi(view, roi_id).geometry.bounds
    px, py = settings.ROI_FOCUS_PADDING * (x1 - x0), settings.ROI_FOCUS_PADDING * (y1 - y0)
    return dict(x0=x0 - px, x1=x1 + px, y0=y0 - py, y1=y1 + py)


@router.get("/layers/{layer_id}/barcodes")
def barcodes(layer_id: str, q: str = Query("", max_length=200), session=Depends(current_session)):
    view = _view(session, layer_id)
    if not hasattr(view.sample, "AGGREGATION_FACTORS"):
        return []
    return view.rois._barcode_suggestions(q)

