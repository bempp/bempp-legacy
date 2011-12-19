// Docstrings ------------------------------------------------------------------

%define IdSet_docstring
"Id set."
%enddef 

%define IdSet_entityId_docstring
"entityId(self, e) -> int

Id of the entity e."
%enddef

%define IdSet_subEntityId_docstring
"subEntityId(self, e, i, codimSub) -> int

Id of i'th subentity of codimension codimSub of entity e of codimension 0."
%enddef


// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (IdSet);
DECLARE_METHOD_DOCSTRING(IdSet, entityId, 0);
DECLARE_METHOD_DOCSTRING(IdSet, subEntityId, 0);

} // namespace Bempp
