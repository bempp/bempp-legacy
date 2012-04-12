// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef bempp_geometry_hpp
#define bempp_geometry_hpp

#include "geometry_type.hpp"
#include "common.hpp"

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <armadillo>

namespace Fiber
{

template <typename ValueType> struct GeometricalData;

} // namespace Fiber

namespace Bempp
{

/** \brief Abstract wrapper of a geometry. */
class Geometry
{
public:
    /** \brief Destructor. */
    virtual ~Geometry() {}

    /** \brief Dimension of the geometry. */
    virtual int dim() const = 0;

    /** \brief Dimension of the space containing the geometry. */
    virtual int dimWorld() const = 0;

    /** \brief Set up geometry of an entity.

      \param[in] corners  Coordinates of the entity's vertices,
                          stored columnwise.

      \param[in] auxData  Auxiliary data necessary for the description of the
                          entity. Interpretation of these data in subclasses of
                          Geometry may vary. They can be used for example to
                          define a *curvilinear* element. */
    virtual void setup(const arma::Mat<ctype>& corners,
                       const arma::Col<char>& auxData) = 0;

    /** \brief Type of the reference element.

     The type can be used to access the <tt>Dune::GenericReferenceElement</tt>.
     */
    virtual GeometryType type() const = 0;

    /** \brief True if the geometry mapping is affine and false otherwise. */
    virtual bool affine() const = 0;

    /** \brief Number of corners of the reference element.

      Since a geometry is a convex polytope the number of corners is a
      well-defined concept. The method is redundant because this
      information is also available via the reference element. It is
      here for efficiency and ease of use.
     */
    virtual int cornerCount() const = 0;

    /** \brief Positions of the geometry corners.
     *
     *  \param[out] c Matrix whose \f$i\f$th column contains the coordinates of
     *  the \f$i\f$th corner. The numbering of corners follows the conventions
     *  of the generic reference element.
     */
    virtual void corners(arma::Mat<ctype>& c) const = 0;

    /** \brief Convert local (logical) to global (physical) coordinates.

      \param[in] local Matrix whose \f$i\f$th column contains the local coordinates of a point \f$x_i \in D\f$.
      \param[out] global Matrix whose \f$i\f$th column contains the global coordinates of \f$x_i\f$, i.e. \f$g(x_i)\f$.
    */
    virtual void local2global(const arma::Mat<ctype>& local,
                              arma::Mat<ctype>& global) const = 0;

    /** \brief Convert global (physical) to local (logical) coordinates.

      \param[in] global Matrix whose \f$i\f$th column contains the global coordinates of a point \f$x_i \in W\f$.
      \param[out] local Matrix whose \f$i\f$th column contains the local coordinates of \f$x_i\f$, i.e. \f$g^{-1}(x_i)\f$.

      \fixme This is going to be tricky to implement for dimGrid < dimWorld.
      Maybe the docstring should say that we convert some sort of *projection*
      of global to local.
    */
    virtual void global2local(const arma::Mat<ctype>& global,
                              arma::Mat<ctype>& local) const = 0;

    /** \brief The factor appearing in the integral transformation formula.

      Let \f$ g : D \to W\f$ denote the transformation described by the Geometry.
      Then the jacobian of the transformation is defined as the
      \f$\textrm{cdim}\times\textrm{mydim}\f$ matrix
      \f[ J_g(x) = \left( \begin{array}{ccc} \frac{\partial g_0}{\partial x_0} &
      \cdots & \frac{\partial g_0}{\partial x_{n-1}} \\
      \vdots & \ddots & \vdots \\ \frac{\partial g_{m-1}}{\partial x_0} &
      \cdots & \frac{\partial g_{m-1}}{\partial x_{n-1}}
      \end{array} \right).\f]
      Here we abbreviated \f$m=\textrm{cdim}\f$ and \f$n=\textrm{mydim}\f$ for ease of
      readability.

      The integration element \f$\mu(x)\f$ for any \f$x\in D\f$ is then defined as
      \f[ \mu(x) = \sqrt{|\det J_g^T(x)J_g(x)|}.\f]

      \param[in]   local       Matrix whose \f$i\f$th column contains the local coordinates of a point \f$x_i \in D\f$.
      \param[out]  int_element  Row vector whose \f$i\f$th entry contains the integration element \f$\mu(x_i)\f$.

      \note Each implementation computes the integration element with optimal
      efficiency. For example in an equidistant structured mesh it may be as
      simple as \f$h^\textrm{mydim}\f$.
    */
    virtual void integrationElement(const arma::Mat<ctype>& local,
                                    arma::Row<ctype>& int_element) const = 0;

    /** \brief Volume of geometry. */
    virtual ctype volume() const = 0;

    /** \brief Center of geometry.
     *
     *  Note that this method is still subject to a change of name and
     *  semantics. At the moment, the center is not required to be the centroid
     *  of the geometry, or even the centroid of its corners. This makes
     *  acceptable the current default implementation, which maps the centroid
     *  of the reference element to the geometry.
     *
     *  We may change the name (and semantic) of the method to centroid() if
     *  Dune's developers find reasonably efficient ways to implement it
     *  properly.
     *
     * \param[out]  c  Coordinates of the center of geometry.
     */
    virtual void center(arma::Col<ctype>& c) const = 0;

    /** \brief Transpose of the Jacobian matrix.
     *
     *  The Jacobian matrix is defined in the documentation of
     *  integrationElement().
     *
     *  \param[in]  local
     *    Matrix whose \f$i\f$th column contains the local coordinates of a point \f$x_i \in D\f$.
     *  \param[out]  jacobian_t
     *    3D array whose \f$i\f$th slice (i.e. jacobian_t(:,:,i)) contains the
     *    transposed Jacobian matrix at \f$x_i\f$, i.e. \f$J_g^T(x_i)\f$.
     */
    virtual void jacobianTransposed(const arma::Mat<ctype>& local,
                                    arma::Cube<ctype>& jacobian_t) const = 0;

    /** \brief Inverse of the transposed Jacobian matrix.
     *
     *  The Jacobian matrix is defined in the documentation of
     *  integrationElement().
     *
     *  \param[in]  local
     *    Matrix whose \f$i\f$th column contains the local coordinates of a point \f$x_i \in D\f$.
     *  \param[out]  jacobian_inv_t
     *    3D array whose \f$i\f$th slice (i.e. jacobian_inv_t(:,:,i)) contains the
     *    inverse of the transposed Jacobian matrix at \f$x_i\f$, i.e. \f$J_g^{-T}(x_i)\f$.
     *
     *  The use of this function is to compute the gradient of some function
     *  \f$f : W \to \textbf{R}\f$ at some position \f$y=g(x)\f$, where
     *  \f$x\in D\f$ and \f$g\f$ the transformation of the Geometry.
     *  When we set \f$\hat{f}(x) = f(g(x))\f$ and apply the chain rule we obtain
     *  \f[\nabla f(g(x)) = J_g^{-T}(x) \nabla \hat{f}(x).\f]
     *
     *  \note In the non-symmetric case \f$\textrm{cdim} \neq \textrm{mydim}\f$, the
     *        pseudoinverse of \f$J_g^T(x)\f$ is returned.
     *        This means that it is inverse for all tangential vectors in
     *        \f$g(x)\f$ while mapping all normal vectors to zero.
     */
    virtual void jacobianInverseTransposed(const arma::Mat<ctype>& local,
                                           arma::Cube<ctype>& jacobian_inv_t) const = 0;

    /** \brief Unit vector normal to the entity.
     *
     *  An exception is thrown if dim() != dimWorld() - 1.
     *
     *  \param[in]  local
     *    Matrix whose \f$i\f$th column contains the local coordinates of a point \f$x_i \in D\f$.
     *  \param[out]  jacobian_t
     *    Matrix whose \f$i\f$th column containts components of a unit vector
     *    normal to the entity at \f$x_i\f$.
     */
    virtual void normal(const arma::Mat<ctype>& local,
                        arma::Mat<ctype>& normal) const = 0;

    virtual void getData(int what, const arma::Mat<ctype>& local,
                         Fiber::GeometricalData<ctype>& data) const = 0;
};

} // namespace Bempp

#endif
