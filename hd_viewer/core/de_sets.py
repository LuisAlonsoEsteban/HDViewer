"""The groups and clusters a differential-expression analysis can compare.

Names shown to users: group names, and "Cluster <label>" for clusters (with a
" (clustering)" suffix when a group already uses that name). Masks align with
the sample's native observations.
"""
from .. import settings

CLUSTER_PREFIX = settings.DE_CLUSTER_PREFIX
SOURCES = ("groups", "clusters")


def _clustered(clustering):
    return clustering if clustering is not None and getattr(clustering, "is_clustered", False) else None


def de_options(source, groups=None, clustering=None):
    """Names of the sets available for DE, in display order."""
    if source == "groups":
        return list(getattr(groups, "_groups", {})) if groups is not None else []
    clustering = _clustered(clustering)
    if clustering is None:
        return []
    names = []
    group_names = getattr(groups, "_groups", {})
    for label in getattr(clustering, "_cluster_labels", []):
        name = f"{CLUSTER_PREFIX}{label}"
        names.append(f"{name} (clustering)" if name in group_names else name)
    return names


def cluster_label(name):
    return name[len(CLUSTER_PREFIX):].removesuffix(" (clustering)")


def de_masks(sample, source, names, groups=None, clustering=None):
    """Boolean membership per set, aligned with the sample's observations."""
    if sample is None:
        return {}
    obs = sample.adata.obs
    if source == "groups":
        if groups is None:
            return {}
        return {name: obs[f"{settings.GROUP_MASK_PREFIX}{name}"].astype(bool).to_numpy()
                for name in names if f"{settings.GROUP_MASK_PREFIX}{name}" in obs.columns}
    labels = obs[_clustered(clustering).cluster_col].astype(str)
    return {name: (labels == cluster_label(name)).to_numpy() for name in names}


def de_colours(source, groups=None, clustering=None):
    if source == "groups":
        return dict(getattr(groups, "_color_map", {})) if groups is not None else {}
    colours = _clustered(clustering)._color_map
    return {name: colours.get(cluster_label(name)) for name in de_options("clusters", groups, clustering)}


def de_inputs(sample, chosen, groups=None, clustering=None):
    """(names, masks, colours, source) for chosen [(source, name)] pairs."""
    names, masks, colours = [], {}, {}
    for source in SOURCES:
        picked = [name for kind, name in chosen if kind == source]
        if picked:
            names += picked
            masks.update(de_masks(sample, source, picked, groups, clustering))
            colours.update(de_colours(source, groups, clustering))
    sources = {source for source, _ in chosen}
    return names, masks, colours, (sources.pop() if len(sources) == 1 else "mixed")
