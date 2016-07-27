"""Implement the interface to index sets."""


class IndexSet(object):
    """Query the index set of a grid view."""

    def __init__(self, impl, view):

        import numpy as np

        self._impl = impl
        self._element_indices = {}
        self._vertex_indices = {}
        self._edge_indices = {}
        self._sub_entity_indices = {}

        for vertex in view.entity_iterator(2):
            self._vertex_indices[vertex] = self._impl.entity_index(vertex._impl)

        for edge in view.entity_iterator(1):
            self._edge_indices[edge] = self._impl.entity_index(edge._impl)

        # Get all the indices and cache them
        for element in view.entity_iterator(0):

            self._element_indices[element] = self._impl.entity_index(element._impl)

            edge_indices = np.zeros(3, dtype='int')
            vertex_indices = np.zeros(3, dtype='int')

            for i in range(3):
                edge_indices[i] = self._impl.sub_entity_index(element._impl, i, 1)
                vertex_indices[i] = self._impl.sub_entity_index(element._impl, i, 2)

            self._sub_entity_indices[element] = {2:vertex_indices, 1:edge_indices}


    def entity_index(self, entity):
        """Return the index of a given entity."""
        
        if entity.codimension == 0:
            return self._element_indices[entity]
        if entity.codimension == 1:
            return self._edge_indices[entity]
        if entity.codimension == 2:
            return self._vertex_indices[entity]
        raise ValueError("Unknonw codimension.")

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
            return ValueError("`Element` must be an entity of codimension 0.")
        return self._sub_entity_indices[element][codim][i]
