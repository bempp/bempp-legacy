// Docstrings ------------------------------------------------------------------

%define Geometry_docstring
"Geometry of an entity."
%enddef

%define Geometry_type_docstring
"Type of the reference element."
%enddef

%define Geometry_affine_docstring
"True if the geometry mapping is affine and false otherwise."
%enddef

%define Geometry_cornerCount_docstring
"Number of corners of the reference element."
%enddef

%define Geometry_getCorners_autodoc_docstring
"corners(self) -> ndarray"
%enddef

%define Geometry_getCorners_docstring
"Positions of the geometry corners.

Returns a 2D array whose ith column contains the coordinates of
the ith corner. The numbering of corners follows the conventions
of the generic reference element."
%enddef

%define Geometry_local2global_autodoc_docstring
"local2global(self, local) -> ndarray"
%enddef

%define Geometry_local2global_docstring
"Convert local (logical) to global (physical) coordinates.

*Parameters:*
   - local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i

Returns a 2D array whose ith column contains the global coordinates of x_i."
%enddef

%define Geometry_global2local_autodoc_docstring
"global2local(self, global_) -> ndarray"
%enddef

%define Geometry_global2local_docstring
"Convert global (physical) to local (logical) coordinates.

*Parameters:*
   - global (ndarray)
        2D array whose ith column contains the global coordinates of a point x_i.

Returns a 2D whose ith column contains the local coordinates of x_i."
%enddef

%define Geometry_getIntegrationElements_autodoc_docstring
"integrationElement(self, local) -> ndarray"
%enddef

%define Geometry_getIntegrationElements_docstring
"The factor mu appearing in the integral transformation formula.

See the documentation of the C++ interface for the definition of mu.

*Parameters:*
   - local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i.

Returns a vector whose ith entry contains the integration element mu(x_i)."
%enddef

%define Geometry_volume_autodoc_docstring
"volume(self) -> float"
%enddef

%define Geometry_volume_docstring
"Volume of geometry."
%enddef

%define Geometry_getCenter_autodoc_docstring
"center(self) -> ndarray"
%enddef

%define Geometry_getCenter_docstring
"Center of geometry.

Note that this method is still subject to a change of name and
semantics. At the moment, the center is not required to be the centroid
of the geometry, or even the centroid of its corners. This makes
acceptable the current default implementation, which maps the centroid
of the reference element to the geometry.

We may change the name (and semantics) of the method to centroid() if
Dune's developers find reasonably efficient ways to implement it
properly.

Returns a vector containing the coordinates of the center of geometry."
%enddef

%define Geometry_getJacobiansTransposed_autodoc_docstring
"jacobiansTransposed(self, local) -> ndarray"
%enddef

%define Geometry_getJacobiansTransposed_docstring
"Transposed Jacobian matrices.

See the documentation of the C++ interface for the definition of the
Jacobian matrix.

*Parameters:*
   - local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i.

Returns a 3D array whose ith slice (i.e. ...(:,:,i)) contains the
transposed Jacobian matrix at x_i."
%enddef

%define Geometry_getJacobianInversesTransposed_autodoc_docstring
"jacobianInversesTransposed(self, local) -> ndarray"
%enddef

%define Geometry_getJacobianInversesTransposed_docstring
"Inverses of the transposed Jacobian matrices.

See the documentation of the C++ interface for the definition of the
Jacobian matrix.

*Parameters:*
   - local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i.

Returns a 3D array whose ith slice (i.e. ...(:,:,i)) contains the
inverse of the transposed Jacobian matrix at x_i.

*Note:* In the non-symmetric case dimGrid != dimWorld the
pseudoinverse of the transposed Jacobian matrix is returned. This
means that it is an inverse for all vectors tangential to the grid
while mapping all normal vectors to zero."
%enddef

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (Geometry);
DECLARE_METHOD_DOCSTRING(Geometry, type, 1);
DECLARE_METHOD_DOCSTRING(Geometry, affine, 1);
DECLARE_METHOD_DOCSTRING(Geometry, cornerCount, 1);
DECLARE_METHOD_DOCSTRING(Geometry, getCorners, 0);
DECLARE_METHOD_DOCSTRING(Geometry, local2global, 0);
DECLARE_METHOD_DOCSTRING(Geometry, global2local, 0);
DECLARE_METHOD_DOCSTRING(Geometry, getIntegrationElements, 0);
DECLARE_METHOD_DOCSTRING(Geometry, volume, 0);
DECLARE_METHOD_DOCSTRING(Geometry, getCenter, 0);
DECLARE_METHOD_DOCSTRING(Geometry, getJacobiansTransposed, 0);
DECLARE_METHOD_DOCSTRING(Geometry, getJacobianInversesTransposed, 0);

} // namespace Bempp
