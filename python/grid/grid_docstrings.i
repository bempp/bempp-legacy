// Docstrings ------------------------------------------------------------------

%define Grid_docstring
"Grid."
%enddef 

%define Grid_dim_docstring
"Dimension of the grid."
%enddef 

%define Grid_dimWorld_docstring
"Dimension of the space containing the grid."
%enddef 

%define Grid_maxLevel_docstring
"Maximum level defined in this grid.

Levels are numbered 0 ... maxLevel() with 0 the coarsest level."
%enddef 

%define Grid_boundarySegmentCount_docstring
"boundarySegmentCount(self) -> int

Number of boundary segments within the macro (level-0) grid."
%enddef 

%define Grid_levelView_docstring
"levelView(self, level) -> GridView

View of the entities on grid level `level'."
%enddef 

%define Grid_leafView_docstring
"leafView(self) -> GridView

View of the leaf entities."
%enddef 

%define Grid_globalIdSet_docstring
"Grid's global id set."
%enddef 

%define Grid_refineGlobally_docstring
"Refine the grid refCount times using the default refinement rule.

This behaves like marking all elements for refinement and then calling
preAdapt(), adapt() and postAdapt(). The state after refineGlobally()
is comparable to the state after postAdapt()."
%enddef 

%define Grid_mark_docstring
"mark(self, refCount, e) -> bool

Mark an entity to be refined/coarsened in a subsequent adapt().

*Arguments:*
    refCount (int)
        Number of subdivisions that should be applied. Negative value
        means coarsening.
    e (Entity)
        Entity of codimension 0 that should be marked.

Returns True if e was marked, False otherwise."
%enddef 

%define Grid_getMark_docstring
"getMark(self, e) -> int

Adaptation mark for entity e of codimension 0."
%enddef 

%define Grid_preAdapt_docstring
"To be called after marking entities, but before calling adapt().

This sets the mightVanish flags of the elements for the next adapt() call.

Returns True if an entity may be coarsened during a subsequent
adapt(), False otherwise."
%enddef 

%define Grid_adapt_docstring
"Refine all positive marked leaf entities, coarsen all negative marked
entities if possible.

Returns True if a least one entity was refined, false otherwise.

The complete adaptation process works as follows:
   - mark entities with the mark() method
   - call preAdapt()
   - if preAdapt() returned true: possibly save current solution
   - call adapt()
   - if adapt() returned true: possibly interpolate the (saved) solution
   - call postAdapt()."
%enddef 

%define Grid_postAdapt_docstring
"To be called when the grid has been adapted and information left over
by the adaptation has been processed.

This removes the isNew flags of the elements from the last adapt()
call."
%enddef 

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (Grid);
DECLARE_METHOD_DOCSTRING(Grid, dim, 1);
DECLARE_METHOD_DOCSTRING(Grid, dimWorld, 1);
DECLARE_METHOD_DOCSTRING(Grid, maxLevel, 1);
DECLARE_METHOD_DOCSTRING(Grid, boundarySegmentCount, 0);
DECLARE_METHOD_DOCSTRING(Grid, levelView, 0);
DECLARE_METHOD_DOCSTRING(Grid, leafView, 0);
DECLARE_METHOD_DOCSTRING(Grid, globalIdSet, 1);
DECLARE_METHOD_DOCSTRING(Grid, refineGlobally, 1);
DECLARE_METHOD_DOCSTRING(Grid, mark, 0);
DECLARE_METHOD_DOCSTRING(Grid, getMark, 0);
DECLARE_METHOD_DOCSTRING(Grid, preAdapt, 1);
DECLARE_METHOD_DOCSTRING(Grid, adapt, 1);
DECLARE_METHOD_DOCSTRING(Grid, postAdapt, 1);

} // namespace Bempp
