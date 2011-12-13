// Docstrings ------------------------------------------------------------------

%define EntityIteratorCodim0_docstring
"Iterator over entities of codimension 0."
%enddef 

%define EntityIteratorCodim0_next_docstring
"next(self) -> Entity

Return the current entity and increment the iterator."
%enddef 

%define EntityIteratorCodim1_docstring
"Iterator over entities of codimension 1."
%enddef 
#define EntityIteratorCodim1_next_docstring EntityIteratorCodim0_next_docstring

%define EntityIteratorCodim2_docstring
"Iterator over entities of codimension 2."
%enddef 
#define EntityIteratorCodim2_next_docstring EntityIteratorCodim0_next_docstring

%define EntityIteratorCodim3_docstring
"Iterator over entities of codimension 3."
%enddef 
#define EntityIteratorCodim3_next_docstring EntityIteratorCodim0_next_docstring

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityIteratorCodim0, EntityIterator<0>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityIteratorCodim0, EntityIterator<0>, next, 0);

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityIteratorCodim1, EntityIterator<1>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityIteratorCodim1, EntityIterator<1>, next, 0);

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityIteratorCodim2, EntityIterator<2>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityIteratorCodim2, EntityIterator<2>, next, 0);

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityIteratorCodim3, EntityIterator<3>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityIteratorCodim3, EntityIterator<3>, next, 0);

} // namespace Bempp
