import bempp as bpp

class ParallelInterface:
    """ An interface used to hold a grid and matching spaces.

        Attributes

        grid : grid
            The grid to use in parallel calculations

        Spaces : dict
            Dictionary of space-name, space pairs. The spaces must be defined
            over the above grid.

    """

    def __init__(self, grid, spaces):
        # Make sure that grid matches the one that is created non locally
        # Dune may not always create grids with nodes in the same order as the
        # original grid if duplicated. So duplicate the original grid too.
        self.__grid = bpp.grid_from_element_data(grid.leaf_view.vertices,
                                                 grid.leaf_view.elements,
                                                 grid.leaf_view.domain_indices)
        self.__spaces = {}
        for spacename in spaces.keys():
            origspace = spaces[spacename]
            self.spaces[spacename] = bpp.function_space(self.grid,
                                                        origspace.kind,
                                                        origspace.order,
                                                        origspace.domains,
                                                        origspace.closed)

    @property
    def spaces(self):
        """ Dictionary of spaces used for parallel calculations"""
        return self.__spaces

    @property
    def grid(self):
        """ Grid used for parallel calculations. """
        return self.__grid


    def __setstate__(self, state):
        gridstate = state['grid']
        vertices = gridstate['vertices'],
        elements = gridstate['elements'],
        domain_indices = gridstate['domain_indices']
        vertices = vertices[0]
        elements = elements[0]
        self.__grid = bpp.grid_from_element_data(vertices,elements,domain_indices)
        self.__spaces = {}
        for spacename in state['spaces'].keys():
            spacestate = state['spaces'][spacename]
            self.__spaces[spacename] =  bpp.function_space(self.grid,
                                                           spacestate['kind'],
                                                           spacestate['order'],
                                                           spacestate['domains'],
                                                           spacestate['closed'])

    def __getstate__(self):
        state = dict()
        state['grid'] = self.grid.__getstate__()
        state['spaces'] = dict()
        for spacename in self.spaces.keys():
            state['spaces'][spacename] = self.spaces[spacename].__getstate__()
        return state
