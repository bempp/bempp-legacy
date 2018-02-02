"""Implementation of a view onto a grid."""


class GridView(object):
    """Provides interfaces to query the entities of a grid."""

    def __init__(self, impl):
        """Constructor. Not to be called directly."""
        self._impl = impl
        self._index_set = None
        self._reverse_element_map = None
        self._reverse_vertex_map = None
        self._reverse_edge_map = None
        self._elements = None
        self._vertices = None
        self._edges = None
        self._connectivity = None
        self._minimum_diameter = None
        self._maximum_diameter = None
        self._reference_edges = None
        self._marked_list = None

    def _create_connectivity_matrices(self):
        """
        Create connectivity matrices between entities.

        This function creates sparse matrices that map vertices to
        elements and edges to elements.

        """
        import numpy as np
        from scipy.sparse import csc_matrix

        entity_count = {0: self.entity_count(0),
                        1: self.entity_count(1),
                        2: self.entity_count(2)}

        edge_indices = []
        edge_element_indices = []
        vertex_indices = []
        vertex_element_indices = []

        ind_set = self.index_set()

        for element in self.entity_iterator(0):

            index = ind_set.entity_index(element)

            for i in range(3):
                # Loop over vertices
                vertex_index = ind_set.sub_entity_index(element, i, 2)
                vertex_indices.append(vertex_index)
                vertex_element_indices.append(index)

            for j in range(3):
                # Loop over edges
                edge_index = ind_set.sub_entity_index(element, j, 1)
                edge_indices.append(edge_index)
                edge_element_indices.append(index)

        # Now create the matrices

        self._connectivity = {}

        self._connectivity['vertices'] = csc_matrix(
            (np.ones(len(vertex_indices), dtype='int64'),
             (vertex_indices, vertex_element_indices)),
            shape=(entity_count[2], entity_count[0]),
            dtype=np.int64)

        self._connectivity['edges'] = csc_matrix(
            (np.ones(len(edge_indices), dtype='int64'),
             (edge_indices, edge_element_indices)),
            shape=(entity_count[1], entity_count[0]),
            dtype=np.int64)

    def entity_count(self, codimension):
        """Return the number of entities of a given codimension."""
        return self._impl.entity_count(codimension)

    def index_set(self):
        """Return an IndexSet object for the GridView."""
        if self._index_set is None:
            from bempp.api.grid.index_set import IndexSet
            self._index_set = IndexSet(self._impl.index_set())
        return self._index_set

    def entity_iterator(self, codimension):
        """Return an entity iterator for a given codimension."""
        from bempp.api.grid.entity_iterator import EntityIterator
        if codimension == 0:
            if self._elements is None:
                self._elements = list(
                    EntityIterator(codimension, self._impl.entity_iterator(
                        codimension)))
            return iter(self._elements)
        if codimension == 1:
            if self._edges is None:
                self._edges = list(
                    EntityIterator(codimension, self._impl.entity_iterator(
                        codimension)))
            return iter(self._edges)
        if codimension == 2:
            if self._vertices is None:
                self._vertices = list(
                    EntityIterator(codimension, self._impl.entity_iterator(
                        codimension)))
            return iter(self._vertices)
        raise ValueError("Unknown codimension.")

    def element_from_index(self, index):
        """Map a given index to the associated element."""
        if self._reverse_element_map is None:
            self._reverse_element_map = [None] * self.entity_count(0)
            ind_set = self.index_set()
            for element in self.entity_iterator(0):
                elem_index = ind_set.entity_index(element)
                self._reverse_element_map[elem_index] = element

        return self._reverse_element_map[index]

    def vertex_from_index(self, index):
        """Map a given index to the associated vertex."""
        if self._reverse_vertex_map is None:
            self._reverse_vertex_map = [None] * self.entity_count(2)
            ind_set = self.index_set()
            for vertex in self.entity_iterator(2):
                vertex_index = ind_set.entity_index(vertex)
                self._reverse_vertex_map[vertex_index] = vertex

        return self._reverse_vertex_map[index]
    
    def edge_from_index(self, index):
        """Map a given index to the associated edge."""
        if self._reverse_edge_map is None:
            self._reverse_edge_map = [None] * self.entity_count(1)
            ind_set = self.index_set()
            for edge in self.entity_iterator(1):
                edge_index = ind_set.entity_index(edge)
                self._reverse_edge_map[edge_index] = edge

        return self._reverse_edge_map[index]

    
    @property
    def reference_edges(self):
        """
        Returns for each triangle the local index of the associated reference edge.

        The reference edge determines along which edge a triangle is next refined
        when it is marked for refinement.

        If the grid has not been created from the refinement of another grid
        the reference edge is just the longest edge in a triangle. Otherwise,
        it is derived from the newest vertex bisection rule.
        """

        if self._reference_edges is None:
            # The grid has never been refined.
            import numpy as np

            index_set = self.index_set()
            
            self._reference_edges = np.zeros(
                self.entity_count(0), dtype='uint32')
            for element in self.entity_iterator(0):
                index = index_set.entity_index(element)
                volumes = np.zeros(3, dtype='float64')
                for i, edge in enumerate(element.sub_entity_iterator(1)):
                    volumes[i] = edge.geometry.volume
                self._reference_edges[index] = np.argmax(volumes)
        return self._reference_edges
                                
    @property
    def vertex_to_element_matrix(self):
        """
        Return the vertex to element matrix.

        The vertex to element matrix is a sparse matrix A, where A[i, j]
        is 1 if vertex i is associated with element j, otherwise A[i, j] = 0.

        """
        if self._connectivity is None:
            self._create_connectivity_matrices()
        return self._connectivity['vertices']

    @property
    def edge_to_element_matrix(self):
        """
        Return the edge to element matrix.

        The dge to element matrix is a sparse matrix A, where A[i, j]
        is 1 if edge i is associated with element j, otherwise A[i, j] = 0.

        """
        if self._connectivity is None:
            self._create_connectivity_matrices()
        return self._connectivity['edges']

    @property
    def dim(self):
        """Return the dimension of the grid."""
        return self._impl.dim

    @property
    def dim_world(self):
        """Return the dimension of the space containing the grid."""
        return self._impl.dim_world

    @property
    def vertices(self):
        """Return a (3 x n_vertices) array with the vertices of the grid."""
        return self._impl.vertices

    @property
    def elements(self):
        """Return a (3 x n_elements) array with the elements of the grid."""
        return self._impl.elements

    @property
    def domain_indices(self):
        """Return a list of domain indices."""
        return self._impl.domain_indices

    @property
    def minimum_element_diameter(self):
        """Return minimum element diameter."""
        if self._minimum_diameter is None:
            self._minimum_diameter = self._impl.minimum_element_diameter
        return self._minimum_diameter

    @property
    def maximum_element_diameter(self):
        """Return maximum element diameter."""
        if self._maximum_diameter is None:
            self._maximum_diameter = self._impl.maximum_element_diameter
        return self._maximum_diameter

    def get_mark(self, element):
        """Return the mark status of a given element."""
        if self._marked_list is None:
            return False
        index = self.index_set().entity_index(element)
        return self._marked_list[index]
        
    def mark(self, element):
        """Mark an element for refinement."""
        if self._marked_list is None:
            import numpy as np
            self._marked_list = set()
        index = self.index_set().entity_index(element)
        self._marked_list.add(index)

    def clear_marks(self):
        """Clear all marks."""
        self._marked_list.clear()

    def refine(self):
        """Refine grid."""
        import numpy as np
        import bempp.api

        bempp.api.log("Grid refinement is experimental and may contain bugs.")
        
        # Create the initial refinement list
        number_of_elements = self.entity_count(0)
        number_of_vertices = self.entity_count(2)
        index_set = self.index_set()

        refinement_list = -1 * np.ones((number_of_elements, 3),
                                       dtype='int')
                                           
        marked = self._marked_list.copy()

        # Now start the refinement loop
        finished = False

        
        edge_connectivity = self.edge_to_element_matrix.tolil().rows
        vertex_number = number_of_vertices
        new_vertices = []
       
        while len(marked) != 0:
            index = marked.pop()
            element = self.element_from_index(index)
            # Check if reference edge is not yet refined:
            if refinement_list[index, self.reference_edges[index]] == -1:
                # Add reference edge to refinement list
                refinement_list[index, self.reference_edges[index]] = vertex_number
                edge_index = index_set.sub_entity_index(
                    element, self.reference_edges[index], 1)
                # Get the coordinate of the mid-point of the reference edge
                edge = self.edge_from_index(edge_index)
                new_vertices.append(
                    .5 * np.sum(edge.geometry.corners, axis=1))
                # Do nothing more if there is no neighboring element at reference edge
                if len(edge_connectivity[edge_index]) == 1:
                    vertex_number += 1
                    continue
                # Get the neighbor
                neighbors = edge_connectivity[edge_index]
                other = neighbors[0] if neighbors[0] != index else neighbors[1]
                # If reference edge not yet marked, add to refinement list
                if refinement_list[other, self.reference_edges[other]] == -1:
                    marked.add(other)

                other_element = self.element_from_index(other)
                for other_edge_index, other_edge in \
                    enumerate(other_element.sub_entity_iterator(1)):
                    if index_set.entity_index(other_edge) == edge_index:
                        refinement_list[other, other_edge_index] = vertex_number
                vertex_number += 1
        # Now create the new vertices and elements.

        vertices = np.zeros((3, number_of_vertices + len(new_vertices)),
                            dtype='float64')
        for i, vertex in enumerate(self.vertices.T):
            vertices[:, i] = vertex

        for i, vertex in enumerate(new_vertices):
            vertices[:, i + number_of_vertices] = vertex

        elements = []
        domain_indices = []
        new_reference_edges = []


        # Mapping of edges to the corresponding node pairs
        edge_to_local_vertices = [[0, 1], [2, 0], [1, 2]]
        
        for i, element_vertices in enumerate(self.elements.T):
            element = self.element_from_index(i)
            ref_index = self.reference_edges[i]
            if refinement_list[i, ref_index] != -1:
                # Element was refined
                if refinement_list[i, (ref_index + 1) % 3] != -1:
                    # Reference edge plus first neighbor to refine
                    v1 = edge_to_local_vertices[ref_index]
                    elements.append(
                        [index_set.sub_entity_index(element, v1[0], 2),
                         refinement_list[i, ref_index],
                         refinement_list[i, (ref_index + 1) % 3]])
                    new_reference_edges.append(0)
                    domain_indices.append(element.domain)
                    elements.append(
                        [refinement_list[i, (ref_index + 1) % 3],
                         refinement_list[i, ref_index],
                         index_set.sub_entity_index(element, 2 - ref_index, 2)])
                    new_reference_edges.append(2)
                    domain_indices.append(element.domain)
                    if refinement_list[i, (ref_index + 2) % 3] == -1:
                        elements.append(
                            [index_set.sub_entity_index(element, 2 - ref_index, 2),
                             refinement_list[i, ref_index],
                             index_set.sub_entity_index(element, v1[1], 2)])
                        new_reference_edges.append(1)
                        domain_indices.append(element.domain)
                        continue
                if refinement_list[i, (ref_index + 2) % 3] != -1:
                    # Reference edge plus second neighbor to refine
                    v1 = edge_to_local_vertices[ref_index]
                    elements.append(
                        [index_set.sub_entity_index(element, v1[1], 2),
                         refinement_list[i, (ref_index + 2) % 3],
                         refinement_list[i, ref_index]])
                    new_reference_edges.append(1)
                    domain_indices.append(element.domain)
                    elements.append(
                        [refinement_list[i, (ref_index + 2) % 3],
                         index_set.sub_entity_index(element, 2 - ref_index, 2),
                         refinement_list[i, ref_index]])
                    new_reference_edges.append(2)
                    domain_indices.append(element.domain)
                    if refinement_list[i, (ref_index + 1) % 3] == -1:
                        elements.append(
                            [index_set.sub_entity_index(element, 2 - ref_index, 2),
                             index_set.sub_entity_index(element, v1[0], 2),
                             refinement_list[i, ref_index]])
                        new_reference_edges.append(0)
                        domain_indices.append(element.domain)
                        continue
                if refinement_list[i, ref_index] != -1:
                    # Only reference edge to refine
                    v1 = edge_to_local_vertices[ref_index]
                    elements.append(
                        [refinement_list[i, ref_index],
                         index_set.sub_entity_index(element, 2 - ref_index, 2),
                         index_set.sub_entity_index(element, v1[0], 2)])
                    new_reference_edges.append(2)
                    domain_indices.append(element.domain)
                    elements.append(
                        [refinement_list[i, ref_index],
                         index_set.sub_entity_index(element, v1[1], 2),
                         index_set.sub_entity_index(element, 2 - ref_index, 2)])
                    new_reference_edges.append(2)
                    domain_indices.append(element.domain)
                    continue
            else:
                # Element is not refined
                elements.append(element_vertices)
                new_reference_edges.append(self.reference_edges[i])
                domain_indices.append(element.domain)
        elements = np.array(elements, dtype='uint32').T
            
        grid = bempp.api.grid_from_element_data(vertices, elements, domain_indices)

        # Now fix the reference edges
        grid.leaf_view.reference_edges[:] = new_reference_edges

        return grid

                    
                    
                    
                    
                    
            
            
            
        
        
                                
                        
                        
                    
                    
                
                
                
            
            
