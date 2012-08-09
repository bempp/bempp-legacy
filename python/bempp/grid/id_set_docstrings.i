// Docstrings ------------------------------------------------------------------

%define IdSet_docstring
"Id set."
%enddef 

%define IdSet_entityId_autodoc_docstring
"entityId(self, e) -> int"
%enddef

%define IdSet_entityId_docstring
"Id of the entity e."
%enddef

%define IdSet_subEntityId_autodoc_docstring
"subEntityId(self, e, i, codimSub) -> int"
%enddef

%define IdSet_subEntityId_docstring
"Id of i'th subentity of codimension codimSub of entity e of codimension 0."
%enddef


// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (IdSet);
DECLARE_METHOD_DOCSTRING(IdSet, entityId, 0);
DECLARE_METHOD_DOCSTRING(IdSet, subEntityId, 0);

} // namespace Bempp
