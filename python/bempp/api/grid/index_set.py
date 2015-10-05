"""Implement the interface to index sets."""

class IndexSet(object):
    """Query the index set of a grid view."""

    def __init__(self, impl):

        self._impl = impl

    def entity_index(self, entity):
        """Return the index of a given entity."""
        return self._impl.entity_index(entity._impl)

    def sub_entity_index(self, element, i, codim):
        """Return the subentity index of an element.

        This method returns the index of a given subentity
        of an element.

        Parameters
        ----------
        element : bempp.api.grid.entity
            Element for which to compute a subentity index.
        i : int
            Number of the subentity.
        codim : int
            Codimension of the subentity.

        Returns
        -------
        id : int
            Index o the subentity.

        Examples
        --------
        The following code returns the index of the first vertex
        of a given element.

        >>> id_set.sub_entity_id(element, 0, 2)

        """
        if element.codimension != 0:
            return ValueError("`Element` must be an entity of codimension0.")
        return self._impl.sub_entity_index(element._impl, i, codim)
