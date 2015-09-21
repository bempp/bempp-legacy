import bempp as bpp

class ParallelInterface:
    
    def __init__(self, grid, spaces):
        # Make sure that grid matches the one that is 
        # Created non locally
        self.grid = bpp.grid_from_element_data(grid.leaf_view.vertices, 
                                               grid.leaf_view.elements,
                                               grid.leaf_view.domain_indices)
        self.spaces = {}
        for spacename in spaces.keys():
            origspace = spaces[spacename]
            self.spaces[spacename] = bpp.function_space(self.grid,
                                                        origspace.kind,
                                                        origspace.order,
                                                        origspace.domains,
                                                        origspace.closed,
                                                        origspace.gridname)

    def __setstate__(self, state):
        gridstate = state['grid']
        vertices = gridstate['vertices'],
        elements = gridstate['elements'],
        domain_indices = gridstate['domain_indices']
        vertices = vertices[0]
        elements = elements[0]
        self.grid = bpp.grid_from_element_data(vertices,elements,domain_indices)
        self.spaces = {}
        for spacename in state['spaces'].keys():
            spacestate = state['spaces'][spacename]
            self.spaces[spacename] =  bpp.function_space(self.grid,
                                                         spacestate['kind'],
                                                         spacestate['order'],
                                                         spacestate['domains'],
                                                         spacestate['closed'],
                                                         spacestate['gridname'])

    def __getstate__(self):
        state = dict()
        state['grid'] = self.grid.__getstate__()
        state['spaces'] = dict()
        for spacename in self.spaces.keys():
            state['spaces'][spacename] = self.spaces[spacename].__getstate__()
        return state
