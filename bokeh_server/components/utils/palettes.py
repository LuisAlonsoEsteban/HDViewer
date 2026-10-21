"""Colours of clusters and custom groups.

The two palettes never share a colour, so a group and a cluster can be told apart
wherever both appear (e.g. a DE comparing a group with a cluster). Greys are left
out: they read as the grey of not-sequenced observations.
"""

from bokeh.palettes import Category20, Category20b, Category20c, Dark2, Muted, Set1

_GREYS = {"#7f7f7f", "#c7c7c7", "#636363", "#969696", "#bdbdbd", "#d9d9d9",
          "#999999", "#666666"}


def _without_greys(colors):
    return [c for c in colors if c.lower() not in _GREYS]


# 54 colours: Leiden often finds more than 20 clusters
CLUSTER_PALETTE = _without_greys([*Category20[20], *Category20b[20], *Category20c[20]])

# 24 bold colours, none of them in CLUSTER_PALETTE
GROUP_PALETTE = _without_greys([*Set1[9], *Dark2[8], *Muted[9]])

assert len(set(CLUSTER_PALETTE)) == len(CLUSTER_PALETTE)
assert len(set(GROUP_PALETTE)) == len(GROUP_PALETTE)
assert not {c.lower() for c in CLUSTER_PALETTE} & {c.lower() for c in GROUP_PALETTE}


def palette_color(palette, index):
    """Colour number `index`; colours repeat only past the end of the palette."""
    return palette[index % len(palette)]
