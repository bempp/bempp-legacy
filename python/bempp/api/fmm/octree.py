"""Python interface for the Octree structure."""


class Octree(object):
    """Octree data structure from a given grid."""

    def __init__(self, grid, levels):
        """Initialize from a given grid."""
        import bempp.core.fmm.octree

        self._grid = grid
        self._impl = bempp.core.fmm.octree.Octree(grid._impl, levels)


    @property
    def grid(self):
        """Return the grid associated with the octree."""
        return self._grid

    @property
    def bounding_box(self):
        """Return grid bounding box."""
        return self._impl.bounding_box

    @property
    def levels(self):
        """Return the number of levels."""
        return self._impl.levels

    def parent(self, n):
        """Return Morton index of parent node."""
        return self._impl.parent(n)

    def children(self, n):
        """Return iterator over the children of node n."""
        return range(*self._impl.children(n))

    def nodes_per_side(self, level):
        """Return number of nodes per side."""
        return self._impl.nodes_per_side(level)

    def nodes_per_level(self, level):
        """Return number of nodes per level."""
        return self._impl.nodes_per_level(level)

    def cube_width(self, level):
        """ Get the width of a cube on a given level."""
        return self._impl.cube_width(level)

    def extended_cube_width(self, level):
        """ Get the extended width of a cube on a given level."""
        return self._impl.extended_cube_width(level)


