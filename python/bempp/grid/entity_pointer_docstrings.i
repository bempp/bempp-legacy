// Docstrings ------------------------------------------------------------------

%define EntityCodim0_docstring
"Entity of codimension 0."
%enddef

%define EntityCodim0_level_docstring
"Entity level."
%enddef

%define EntityCodim0_geometry_docstring
"Geometry of this entity.

This object gives, among other things, the map from a reference element to
world coordinates."
%enddef

%define EntityCodim0_type_docstring
"Type of the reference element."
%enddef

%define EntityCodim0_subEntityCount_docstring
"Number of subentities of codimension codimSub."
%enddef

%define EntityCodim0_subEntities_autodoc_docstring
"subEntities(self, codimSub) -> EntityIterator"
%enddef

%define EntityCodim0_subEntities_docstring
"Iterator over subentities of codimension codimSub.

*Note:* codimSub must be greater than 0 and less than the dimension of
the grid."
%enddef

%define EntityCodim0_father_autodoc_docstring
"father(self) -> Entity"
%enddef

%define EntityCodim0_father_docstring
"Inter-level access to father entity on the next-coarser grid.

The given entity resulted directly from a subdivision of its father
entity. For macro (level-0) elements father() returns None."
%enddef

%define EntityCodim0_hasFather_docstring
"True if entity has a father entity which can be accessed using the father()
method."
%enddef

%define EntityCodim0_isLeaf_docstring
"True if the entity is contained in the leaf grid."
%enddef

%define EntityCodim0_isRegular_docstring
"True if the element is of regular type in red/green type refinement.

In bisection or hanging node refinement this is always true."
%enddef

%define EntityCodim0_sons_autodoc_docstring
"sons(self, maxLevel) -> EntityIterator"
%enddef

%define EntityCodim0_sons_docstring
"Iterator over elements that resulted from (recursive) subdivision of this element.

*Parameters:*
   - maxlevel (int)
        Iterator does not stop at elements with level greater than maxLevel."
%enddef

%define EntityCodim0_isNew_docstring
"True if the entity has been created during the last call to adapt()."
%enddef

%define EntityCodim0_mightVanish_docstring
"True if the entity might disappear during the next call to adapt().

If the method returns False, the entity is guaranteed to still be present after
adaptation."
%enddef


%define EntityCodim1_docstring
"Entity of codimension 1."
%enddef

#define EntityCodim1_level_docstring EntityCodim0_level_docstring
#define EntityCodim1_geometry_docstring EntityCodim0_geometry_docstring
#define EntityCodim1_type_docstring EntityCodim0_type_docstring

%define EntityCodim2_docstring
"Entity of codimension 2."
%enddef

#define EntityCodim2_level_docstring EntityCodim0_level_docstring
#define EntityCodim2_geometry_docstring EntityCodim0_geometry_docstring
#define EntityCodim2_type_docstring EntityCodim0_type_docstring

%define EntityCodim3_docstring
"Entity of codimension 3."
%enddef

#define EntityCodim3_level_docstring EntityCodim0_level_docstring
#define EntityCodim3_geometry_docstring EntityCodim0_geometry_docstring
#define EntityCodim3_type_docstring EntityCodim0_type_docstring


// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_TEMPLATE_CLASS_DOCSTRING(EntityCodim0, EntityPointer<0>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, level, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, geometry, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, type, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, subEntityCount, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, subEntities, 0);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, father, 0);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, hasFather, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, isLeaf, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, isRegular, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, sons, 0);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, isNew, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim0, EntityPointer<0>, mightVanish, 1);

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityCodim1, EntityPointer<1>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim1, EntityPointer<1>, level, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim1, EntityPointer<1>, geometry, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim1, EntityPointer<1>, type, 1);

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityCodim2, EntityPointer<2>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim2, EntityPointer<2>, level, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim2, EntityPointer<2>, geometry, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim2, EntityPointer<2>, type, 1);

DECLARE_TEMPLATE_CLASS_DOCSTRING (EntityCodim3, EntityPointer<3>);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim3, EntityPointer<3>, level, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim3, EntityPointer<3>, geometry, 1);
DECLARE_TEMPLATE_METHOD_DOCSTRING(EntityCodim3, EntityPointer<3>, type, 1);

} // namespace Bempp
