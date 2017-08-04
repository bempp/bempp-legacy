"""Define the interaction list for a given octree node."""


class InteractionList(object):

    def __init__(self, octree, node_index, level):
        """Define an interaction list for a node index at a given level."""
        from bempp.core.fmm import interaction_list

        self._impl = interaction_list.InteractionList(octree._impl, node_index, level)

    def __iter__(self):
        """Return self as iterator"""
        return self

    def next(self):
        """Next element"""
        return self._impl.next()

    def __next__(self):
        """Next element"""
        return self.next()

    



