"""Colours of clusters and custom groups.

The two palettes never share a colour, so a group and a cluster can be told apart
wherever both appear (e.g. a DE comparing a group with a cluster). Greys are left
out: they read as the grey of not-sequenced observations.
"""

from ...settings import CLUSTER_PALETTE, GROUP_PALETTE

assert len(set(CLUSTER_PALETTE)) == len(CLUSTER_PALETTE)
assert len(set(GROUP_PALETTE)) == len(GROUP_PALETTE)
assert not {c.lower() for c in CLUSTER_PALETTE} & {c.lower() for c in GROUP_PALETTE}

def palette_color(palette, index):
    """Colour number `index`; colours repeat only past the end of the palette."""
    return palette[index % len(palette)]
