// Docstrings ------------------------------------------------------------------

%define IndexSet_docstring
"Index set."
%enddef 

%define IndexSet_entityIndex_docstring
"entityIndex(self, e) -> int

Index of the entity e.

The result of calling this method with an entity that is not in the
index set is undefined.

Returns an integer in the range 0 ... (max number of entities in set - 1)."
%enddef

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (IndexSet);
DECLARE_METHOD_DOCSTRING(IndexSet, entityIndex, 0);

} // namespace Bempp
