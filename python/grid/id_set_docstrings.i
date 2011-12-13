// Docstrings ------------------------------------------------------------------

%define IdSet_docstring
"Id set."
%enddef 

%define IdSet_entityId_docstring
"entityId(self, e) -> int

Id of the entity e."
%enddef

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (IdSet);
DECLARE_METHOD_DOCSTRING(IdSet, entityId, 0);

} // namespace Bempp
