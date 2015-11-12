"""Implement interfaces to id sets."""


class IdSet(object):
    """Query the id set of a grid."""

    def __init__(self, impl):
        self._impl = impl

    def entity_id(self, entity):
        """Return the id of an entity."""
        return self._impl.entity_id(entity._impl)

    def sub_entity_id(self, element, i, codim):
        """Return the subentity id of an element.

        This method returns the id of a given subentity
        of an element.

        Parameters
        ----------
        element : bempp.api.grid.entity
            Element for which to compute a subentity id.
        i : int
            Number of the subentity.
        codim : int
            Codimension of the subentity.

        Returns
        -------
        id : int
            Id o the subentity.

        Examples
        --------
        The following code returns the id of the first vertex.

        >>> id_set.sub_entity_id(element, 0, 2)

        """
        if element.codimension != 0:
            return ValueError("`Element` must be an entity of codimension0.")
        return self._impl.sub_entity_id(element._impl, i, codim)

