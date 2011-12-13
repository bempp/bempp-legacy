// Docstrings ------------------------------------------------------------------

%define GridView_docstring
"Grid view.

This class provides a means to access a specific subset of entities,
for example the leaf entities or the entities on a certain refinement
level."
%enddef 

%define GridView_indexSet_docstring
"The index set."
%enddef 

%define GridView_entityCount_docstring
"Number of entities with certain characteristics.

*Arguments:*
    codim (int)
        codimension
    type (GeometryType)
        geometry type
"
%enddef 

%define GridView_containsEntity_docstring
"containsEntity(self, e) -> bool

True if the entity e belongs to this grid view.

If e is not an element of the grid, then the result of containsEntity() is 
undefined."
%enddef 

%define GridView_entities_docstring
"entities(self, codim) -> EntityIterator

Iterator over entities of codimension codim contained in this view."
%enddef 

// TODO: handle the optional param properly 
%define GridView_vtkWriter_docstring
"vtkWriter(self, dm = 'conforming') -> VtkWriter

VtkWriter for this grid view.

*Arguments:*
    dm (string, optional)
        data mode, `conforming' (default) or `nonconforming'

See the documentation of Dune::VTK::DataMode for more information about the 
implications of the choice of data mode."
%enddef 

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING(GridView);
DECLARE_METHOD_DOCSTRING(GridView, indexSet, 1);
DECLARE_METHOD_DOCSTRING(GridView, entityCount, 1);
DECLARE_METHOD_DOCSTRING(GridView, containsEntity, 0);
DECLARE_METHOD_DOCSTRING(GridView, entities, 0);
DECLARE_METHOD_DOCSTRING(GridView, vtkWriter, 0);

} // namespace Bempp
