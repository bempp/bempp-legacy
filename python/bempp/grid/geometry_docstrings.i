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

%define Geometry_corners_docstring
"corners(self) -> ndarray

Positions of the geometry corners.

Returns a 2D array whose ith column contains the coordinates of
the ith corner. The numbering of corners follows the conventions
of the generic reference element."
%enddef

%define Geometry_local2global_docstring
"local2global(self, local) -> ndarray

Convert local (logical) to global (physical) coordinates.

*Arguments:*
    local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i

Returns a 2D array whose ith column contains the global coordinates of x_i."
%enddef

%define Geometry_global2local_docstring
"global2local(self, global_) -> ndarray

Convert global (physical) to local (logical) coordinates.

*Arguments:*
    global (ndarray)
        2D array whose ith column contains the global coordinates of a point x_i.
      
Returns a 2D whose ith column contains the local coordinates of x_i."
%enddef

%define Geometry_integrationElements_docstring
"integrationElement(self, local) -> ndarray

The factor mu appearing in the integral transformation formula.

See the documentation of the C++ interface for the definition of mu.

*Arguments:*
    local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i.

Returns a vector whose ith entry contains the integration element mu(x_i)."
%enddef

%define Geometry_volume_docstring
"volume(self) -> float

Volume of geometry."
%enddef

%define Geometry_center_docstring
"center(self) -> ndarray

Center of geometry.

Note that this method is still subject to a change of name and
semantics. At the moment, the center is not required to be the centroid
of the geometry, or even the centroid of its corners. This makes
acceptable the current default implementation, which maps the centroid
of the reference element to the geometry.

We may change the name (and semantic) of the method to centroid() if
Dune's developers find reasonably efficient ways to implement it
properly.

Returns a vector containing the coordinates of the center of geometry."
%enddef

%define Geometry_jacobiansTransposed_docstring
"jacobiansTransposed(self, local) -> ndarray

Transposed Jacobian matrices.

See the documentation of the C++ interface for the definition of the
Jacobian matrix.

*Arguments:*
    local (ndarray)
        2D array whose ith column contains the local coordinates of a point x_i.

Returns a 3D array whose ith slice (i.e. ...(:,:,i)) contains the
transposed Jacobian matrix at x_i."
%enddef

%define Geometry_jacobianInversesTransposed_docstring
"jacobianInversesTransposed(self, local) -> ndarray

Inverses of the transposed Jacobian matrices.

See the documentation of the C++ interface for the definition of the
Jacobian matrix.

*Arguments:*
    local (ndarray)
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
DECLARE_METHOD_DOCSTRING(Geometry, corners, 0);
DECLARE_METHOD_DOCSTRING(Geometry, local2global, 0);
DECLARE_METHOD_DOCSTRING(Geometry, global2local, 0);
DECLARE_METHOD_DOCSTRING(Geometry, integrationElements, 0);
DECLARE_METHOD_DOCSTRING(Geometry, volume, 0);
DECLARE_METHOD_DOCSTRING(Geometry, center, 0);
DECLARE_METHOD_DOCSTRING(Geometry, jacobiansTransposed, 0);
DECLARE_METHOD_DOCSTRING(Geometry, jacobianInversesTransposed, 0);

} // namespace Bempp
