"""Compact viewport lookup over distinct aggregation-block centres."""
import numpy as np
import shapely


class BlockViewportIndex:
    def __init__(self, grid):
        # One positional row map per level; the spatial tree contains only one
        # point per block, rather than repeating a point for every native bin.
        self.grid = grid
        cols = grid.grid_x.to_numpy()
        rows = grid.grid_y.to_numpy()
        self.order = np.lexsort((rows, cols))
        ordered_cols, ordered_rows = cols[self.order], rows[self.order]
        boundaries = np.flatnonzero((ordered_cols[1:] != ordered_cols[:-1])
                                   | (ordered_rows[1:] != ordered_rows[:-1])) + 1
        self.offsets = np.r_[0, boundaries, len(grid)] if len(grid) else np.array([0])
        first = self.order[self.offsets[:-1]]
        centres = shapely.points(grid.x.to_numpy()[first], grid.y.to_numpy()[first])
        self.tree = shapely.STRtree(centres)

    def query(self, left, right, bottom, top):
        blocks = self.tree.query(shapely.box(left, bottom, right, top))
        if not len(blocks):
            return np.empty(0, dtype=np.intp)
        if len(blocks) == len(self.offsets) - 1:
            return np.arange(len(self.order))
        starts = self.offsets[blocks]
        lengths = self.offsets[blocks + 1] - starts
        # Expand selected blocks to their complete native membership without a
        # sample-sized boolean mask. Restore native order for categorical ties,
        # groupby ordering, observation alignment, and selection.
        ends = np.cumsum(lengths)
        positions = np.repeat(starts - (ends - lengths), lengths) + np.arange(ends[-1])
        return np.sort(self.order[positions])
