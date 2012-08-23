// Docstrings ------------------------------------------------------------------

%define GridFactory_docstring
"Grid factory.

This class provides static member functions to construct grids on the
fly and to import grids from existing files."
%enddef 

%define GridFactory_createStructuredGrid_autodoc_docstring
"createStructuredGrid(topology, lowerLeft, upperRight, nElements) -> Grid"
%enddef

%define GridFactory_createStructuredGrid_docstring
"Construct a regular structured grid.

*Arguments:*
    topology (string)
        Topology of the grid to be constructed (one of `linear',
        `triangular', `quadrilateral', `hybrid').
    lowerLeft (tuple)
        Coordinates of the lower left corner of the grid.
    upperRight (tuple)
        Coordinates of the upper right corner of the grid.
    nElements (tuple)
        Number of grid subdivisions in each direction.

This function constructs a regular structured grid. Its dimension,
dimGrid, and the dimension of the surrounding space, dimWorld, are
determined from the parameter `topology'. The constructed grid covers
the \p dimGrid-dimensional cube

    [lowerLeft(0) upperRight(0)] x [lowerLeft(1) upperRight(1)] 
        x ... x [lowerLeft(dimGrid-1), upperRight(dimGrid-1)].

The last (dimWorld - dimGrid) dimensions of all grid points are set to zero.

Each side of the cube parallel to the nth coordinate axis is subdivided into nElements(n) segments.

*Note:* Currently only grids with triangular topology are supported."
%enddef 

%define GridFactory_importGmshGrid_autodoc_docstring
"importGmshGrid(topology, fileName, verbose = True, 
    insertBoundarySegments = False) -> Grid"
%enddef

%define GridFactory_importGmshGrid_docstring
"Import grid from a file in Gmsh format.

*Arguments:*
    topology (string)
        Topology of the grid to be constructed (one of `linear', 
        `triangular', `quadrilateral', `hybrid').
    fileName (string)
        Name of the Gmsh file.
    verbose (bool, default: True)
        Output diagnostic information.
    insertBoundarySegments (bool, default: False)
        Insert boundary segments.

See http://geuz.org/gmsh for information about the Gmsh file format.
See Dune::GmshReader documentation for information about the supported
Gmsh features."
%enddef 

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (GridFactory);
DECLARE_METHOD_DOCSTRING(GridFactory, createStructuredGrid, 0);
DECLARE_METHOD_DOCSTRING(GridFactory, importGmshGrid, 0);

} // namespace Bempp
