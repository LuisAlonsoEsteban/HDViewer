"""Custom groups of observations, without any UI.

A group is defined by gene/metadata conditions (combined with AND/OR) or by a
fixed membership (e.g. a region of interest). Memberships are stored in
``adata.obs['grp_<name>']`` and the shown groups' labels in
``adata.obs[GROUPS_COL]``. The Bokeh UI subclasses `GroupsState`
(board.layout.layer.groups) and the HTTP API uses it directly.
"""
import numpy as np
import pandas as pd

from .. import settings
from ..startup.utils.progress import report

OPERATORS = [">", "<", ">=", "<=", "=="]
GROUPS_COL = settings.GROUPS_COL
GROUP_PALETTE = settings.GROUP_PALETTE


class GroupsState:
    def __init__(self, sample=None):
        self._sample = sample
        self._groups: dict = {}   # group_name → [conditions]
        self._active_groups: set = set()
        self._fixed_masks = {}
        self._group_origins = {}
        self._change_listeners = []
        self._color_map: dict = {}
        # Bumped whenever group masks are (re)written: the active-groups label
        # column only needs rebuilding when masks or the active set change
        self._mask_version = 0
        self._groups_col_signature = None

    @property
    def _base_sample(self):
        return self._sample

    @property
    def sample(self):
        return self._base_sample

    @property
    def has_groups(self) -> bool:
        return GROUPS_COL in self.sample.adata.obs.columns and len(self._groups) > 0

    # -- listeners -------------------------------------------------------------

    def add_change_listener(self, callback):
        if callback not in self._change_listeners:
            self._change_listeners.append(callback)

    def remove_change_listener(self, callback):
        if callback in self._change_listeners:
            self._change_listeners.remove(callback)

    def _notify_group_change(self):
        for callback in list(self._change_listeners):
            callback()

    # -- conditions ------------------------------------------------------------

    def check_condition(self, col, op, val):
        """Validate one condition; returns a hint (gene range) or raises ValueError."""
        adata = self.sample.adata
        if not col or not str(val).strip():
            raise ValueError("Gene and value are required.")
        if op not in OPERATORS:
            raise ValueError(f"Unknown operator: {op}")
        if col in adata.var_names:
            try:
                float(val)
            except ValueError:
                raise ValueError("Gene expression conditions require a numeric value.") from None
            expr = adata.X[:, adata.var_names.get_loc(col)]
            expr = expr.toarray().ravel() if hasattr(expr, "toarray") else np.asarray(expr).ravel()
            return f"{col}: max={float(expr.max()):.2f}, mean={float(expr.mean()):.3f}"
        if col in adata.obs.columns:
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                try:
                    float(val)
                except ValueError:
                    raise ValueError("Numeric metadata requires a numeric value.") from None
            elif op != "==":
                raise ValueError("Only == is supported for non-numeric metadata values.")
            return ""
        raise ValueError("Unknown gene or metadata column.")

    def _build_mask(self, conditions, adata) -> pd.Series:
        obs = adata.obs
        var_idx = {g: i for i, g in enumerate(adata.var_names)}
        mask = None

        for col, op, val, logic in conditions:
            report(f"Evaluating condition: {col} {op} {val}")
            if col in var_idx:
                expr = adata.X[:, var_idx[col]]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                series = pd.Series(expr, index=obs.index)
                try:
                    val_f = float(val)
                except ValueError as exc:
                    raise ValueError(f"Gene condition '{col}' requires a numeric value.") from exc
                ops = {'>': series > val_f, '<': series < val_f,
                       '>=': series >= val_f, '<=': series <= val_f}
                cond_mask = ops.get(op, series == val_f)
            elif col in obs.columns:
                series = obs[col]
                if pd.api.types.is_numeric_dtype(series):
                    try:
                        val_f = float(val)
                    except ValueError as exc:
                        raise ValueError(f"Column '{col}' requires a numeric value for operator '{op}'.") from exc
                    ops = {'>': series > val_f, '<': series < val_f,
                           '>=': series >= val_f, '<=': series <= val_f}
                    cond_mask = ops.get(op, series == val_f)
                else:
                    if op not in ("==",):
                        raise ValueError(f"Operator '{op}' is only supported for numeric columns; use == for '{col}'.")
                    cond_mask = series.astype(str) == str(val)
            else:
                continue

            if mask is None:
                mask = cond_mask
            elif logic == "AND":
                mask = mask & cond_mask
            else:
                mask = mask | cond_mask

        if mask is None:
            return pd.Series([False] * len(obs), index=obs.index)
        return mask

    # -- membership --------------------------------------------------------------

    def group_for_roi(self, roi_id):
        return next((name for name, origin in self._group_origins.items() if origin == roi_id), None)

    def _count_spots_bool(self, cname: str) -> int:
        col_key = f"{settings.GROUP_MASK_PREFIX}{cname}"
        obs = self.sample.adata.obs
        if col_key not in obs.columns:
            return 0
        return int(obs[col_key].astype(bool).sum())

    def _prepare_groups(self, groups, fixed_masks=None):
        """Prepare membership and labels without changing observations (thread-safe)."""
        adata = self.sample.adata
        fixed_masks = self._fixed_masks if fixed_masks is None else fixed_masks
        masks = {}
        labels = pd.Series(None, index=adata.obs_names, dtype=object)
        for gname, conditions in groups.items():
            report(f"Evaluating group '{gname}' ({adata.n_obs:,} observations)...")
            key = f"{settings.GROUP_MASK_PREFIX}{gname}"
            if key in adata.obs:
                mask = adata.obs[key].astype(bool)
            elif gname in fixed_masks:
                mask = fixed_masks[gname].reindex(adata.obs_names, fill_value=False)
            else:
                mask = self._build_mask(conditions, adata)
            masks[gname] = mask
            labels.loc[mask] = gname
            report(f"Group '{gname}': {int(mask.sum()):,} members")
        report("Preparing group labels...")
        return masks, labels.astype('category')

    def store_groups(self, prepared=None):
        """Write memberships, colours and labels into the sample's observations."""
        adata = self.sample.adata
        masks, labels = prepared if prepared is not None else self._prepare_groups(self._groups)
        self._mask_version += 1
        for gname, mask in masks.items():
            adata.obs[f"{settings.GROUP_MASK_PREFIX}{gname}"] = mask.values
            if gname not in self._color_map:
                self._color_map[gname] = GROUP_PALETTE[len(self._color_map) % len(GROUP_PALETTE)]
        adata.obs[GROUPS_COL] = labels
        if hasattr(self.sample, 'clear_multiscale_cache'):
            self.sample.clear_multiscale_cache()

    def _apply_groups_to_adata(self, prepared=None):
        """Store the groups (the Bokeh UI also recolours its map here)."""
        self.store_groups(prepared)

    def _rollback(self, name):
        self._groups.pop(name, None)
        self._fixed_masks.pop(name, None)
        self._group_origins.pop(name, None)
        self._active_groups.discard(name)
        self.sample.adata.obs.drop(columns=[f"{settings.GROUP_MASK_PREFIX}{name}"], errors="ignore", inplace=True)

    def prepare_condition_group(self, name, conditions):
        """Validate and compute a new condition group (thread-safe); pass the
        result to `commit_condition_group`."""
        name = (name or "").strip() or "Group A"
        if name in self._groups:
            raise ValueError(f"Group '{name}' already exists.")
        if not conditions:
            raise ValueError("Add at least one condition first.")
        conditions = [tuple(condition) for condition in conditions]
        for col, op, val, _ in conditions:
            self.check_condition(col, op, val)
        return name, conditions, self._prepare_groups({**self._groups, name: conditions})

    def commit_condition_group(self, name, conditions, prepared):
        self._groups[name] = conditions
        self._active_groups.add(name)
        try:
            self._apply_groups_to_adata(prepared)
        except Exception:
            self._rollback(name)
            raise
        self._notify_group_change()
        return name

    def add_fixed_group(self, name, mask, origin_roi_id=None, prepared=None):
        """A group with a fixed membership (e.g. a ROI); returns its unique name."""
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (self.sample.adata.n_obs,) or not mask.any():
            raise ValueError("A group needs a nonempty mask aligned with the sample observations.")
        if origin_roi_id is not None and self.group_for_roi(origin_roi_id):
            return self.group_for_roi(origin_roi_id)
        base, suffix = name, 2
        while name in self._groups:
            name = f"{base} ({suffix})"
            suffix += 1
        self._fixed_masks[name] = pd.Series(mask.copy(), index=self.sample.adata.obs_names.copy())
        self._groups[name] = []
        self._group_origins[name] = origin_roi_id
        self._active_groups.add(name)
        try:
            self._apply_groups_to_adata(prepared)
        except Exception:
            self._rollback(name)
            raise
        return name

    def forget_group(self, name):
        """Remove a group's state and membership column; True if groups remain."""
        report(f"Removing group '{name}'...")
        del self._groups[name]
        self._mask_version += 1
        self._fixed_masks.pop(name, None)
        self._group_origins.pop(name, None)
        self._active_groups.discard(name)
        self.sample.adata.obs.drop(columns=[f"{settings.GROUP_MASK_PREFIX}{name}"], errors="ignore", inplace=True)
        if hasattr(self.sample, 'clear_multiscale_cache'):
            self.sample.clear_multiscale_cache()
        return bool(self._groups)

    def members(self, name):
        key = f"{settings.GROUP_MASK_PREFIX}{name}"
        obs = self.sample.adata.obs
        if key not in obs:
            return set()
        return set(self.sample.adata.obs_names[obs[key].astype(bool)])

    def _rebuild_custom_groups_col(self):
        adata = self.sample.adata
        self._mask_version += 1
        col = pd.Series(None, index=adata.obs.index, dtype=object)
        for gname in self._groups:
            col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
            if col_key in adata.obs.columns:
                col[adata.obs[col_key].astype(bool)] = gname
        adata.obs[GROUPS_COL] = col.astype('category')

    def refresh_active_labels(self):
        """``adata.obs[GROUPS_COL]`` holds the shown groups only (rebuilt only when
        masks or the shown set changed: it spans every native observation)."""
        adata = self.sample.adata
        signature = (id(adata), self._mask_version, tuple(self._groups), tuple(sorted(self._active_groups)))
        if signature == self._groups_col_signature:
            return
        col = pd.Series(None, index=adata.obs.index, dtype=object)
        for gname in self._groups:
            if gname in self._active_groups:
                col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
                if col_key in adata.obs.columns:
                    col[adata.obs[col_key].astype(bool).to_numpy()] = gname
        adata.obs[GROUPS_COL] = col
        self._groups_col_signature = signature

    def visible_mask(self, obs_names):
        """Displayed rows (by observation name) belonging to a shown group."""
        adata = self.sample.adata
        names = pd.Series(obs_names, dtype=object)
        visible = pd.Series(False, index=names.index)
        for gname in self._active_groups:
            col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
            if col_key in adata.obs.columns:
                visible |= names.map(adata.obs[col_key].astype(bool)).fillna(False).astype(bool)
        return visible.to_numpy(dtype=bool)

    def group_summary(self):
        return [dict(name=name, colour=self._color_map.get(name, "#cccccc"), size=self._count_spots_bool(name),
                     visible=name in self._active_groups, origin=self._group_origins.get(name),
                     conditions=[list(condition) for condition in conditions])
                for name, conditions in self._groups.items()]
