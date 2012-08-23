// Docstrings ------------------------------------------------------------------

%define GeometryType_docstring
"Unique label for each type of entities that can occur in grids."
%enddef 

%define GeometryType_GeometryType_autodoc_docstring
"__init__(*args) -> GeometryType"
%enddef

%define GeometryType_GeometryType_docstring
"Constructor.

*Variants:*
    __init__(self) -> GeometryType
        Construct the geometry type `none'.

    __init__(self, dim) -> GeometryType
        Construct a vertex or segment.
        Arguments:
            dim (int)
                Dimension, must be either 0 (vertex) or 1 (segment).

    __init__(self, basicType, dim) -> GeometryType
        Construct a general geometry type.
        Arguments:
            basicType (string)
                Basic geometry type, one of `simplex', `cube', `pyramid', 
                `prism', `extended' or `none'.
            dim (int)
                Dimension.
      
"
%enddef

%define GeometryType_makeVertex_docstring
"Make a vertex"
%enddef 

%define GeometryType_makeLine_docstring
"Make a line segment"
%enddef 

%define GeometryType_makeTriangle_docstring
"Make a triangle"
%enddef 

%define GeometryType_makeQuadrilateral_docstring
"Make a quadrilateral"
%enddef 

%define GeometryType_makeTetrahedron_docstring
"Make a tetrahedron"
%enddef 

%define GeometryType_makePyramid_docstring
"Make a pyramid"
%enddef 

%define GeometryType_makePrism_docstring
"Make a prism"
%enddef 

%define GeometryType_makeHexahedron_docstring
"Make a hexahedron"
%enddef 

%define GeometryType_makeSimplex_docstring
"Make a simplex of dimension dim."
%enddef 

%define GeometryType_makeCube_docstring
"Make a hypercube of dimension dim."
%enddef 

%define GeometryType_makeNone_docstring
"Make a singular of dimension dim."
%enddef 

%define GeometryType_isVertex_docstring
"Return true if entity is a vertex"
%enddef 

%define GeometryType_isLine_docstring
"Return true if entity is a line segment"
%enddef 

%define GeometryType_isTriangle_docstring
"Return true if entity is a triangle"
%enddef 

%define GeometryType_isQuadrilateral_docstring
"Return true if entity is a quadrilateral."
%enddef 

%define GeometryType_isTetrahedron_docstring
"Return true if entity is a tetrahedron."
%enddef 

%define GeometryType_isPyramid_docstring
"Return true if entity is a pyramid."
%enddef 

%define GeometryType_isPrism_docstring
"Return true if entity is a prism."
%enddef 

%define GeometryType_isHexahedron_docstring
"Return true if entity is a hexahedron."
%enddef 

%define GeometryType_isSimplex_docstring
"Return true if entity is a simplex of any dimension."
%enddef 

%define GeometryType_isCube_docstring
"Return true if entity is a cube of any dimension"
%enddef 

%define GeometryType_isNone_docstring
"Return true if entity is a singular of any dimension"
%enddef 

%define GeometryType_dim_docstring
"Return dimension of the type."
%enddef

// Declarations ----------------------------------------------------------------

namespace Dune
{

DECLARE_CLASS_DOCSTRING(GeometryType);
DECLARE_METHOD_DOCSTRING(GeometryType, GeometryType, 0);
DECLARE_METHOD_DOCSTRING(GeometryType, makeVertex, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeLine, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeTriangle, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeQuadrilateral, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeTetrahedron, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makePyramid, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makePrism, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeHexahedron, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeSimplex, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeCube, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, makeNone, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isVertex, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isLine, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isTriangle, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isQuadrilateral, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isTetrahedron, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isPyramid, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isPrism, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isHexahedron, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isSimplex, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isCube, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, isNone, 1);
DECLARE_METHOD_DOCSTRING(GeometryType, dim, 1);

} // namespace Dune
