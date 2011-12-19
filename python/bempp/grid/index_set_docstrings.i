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

%define IndexSet_subEntityIndex_docstring
"subEntityIndex(self, e, i, codimSub) -> int

Index of i'th subentity of codimension codimSub of entity e of codimension 0."
%enddef


// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (IndexSet);
DECLARE_METHOD_DOCSTRING(IndexSet, entityIndex, 0);
DECLARE_METHOD_DOCSTRING(IndexSet, subEntityIndex, 0);

} // namespace Bempp
