// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_geometry_decl_hpp
#define bempp_geometry_decl_hpp

#include "geometry_type_decl.hpp"
#include "common.hpp"

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <armadillo>

namespace Bempp
{

/** Abstract wrapper of a geometry */
class Geometry
{
public:
    /** Destructor */
    virtual ~Geometry() {}

    /** \brief Type of the reference element. The type can
      be used to access the Dune::GenericReferenceElement.
     */
    virtual GeometryType type() const = 0;

    /** \brief True if the geometry mapping is affine and false otherwise */
    virtual bool affine() const = 0;

    /** \brief Number of corners of the reference element.
     *
      Since a geometry is a convex polytope the number of corners is a
      well-defined concept. The method is redundant because this
      information is also available via the reference element. It is
      here for efficiency and ease of use.
     */
    virtual int n_corners() const = 0;

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
     *
     *  \return \f$J_g^T(x)\f$
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
};

/** Wrapper of a Dune geometry of type DuneGeometry */

template<typename DuneGeometry>
class ConcreteGeometry : public Geometry
{
private:
    const DuneGeometry* m_dune_geometry;

    /** \brief Default constructor

    \internal Should be used only by friend classes that call setDuneGeometry() later on. */
    ConcreteGeometry() : m_dune_geometry(0) {
    }

    void setDuneGeometry(const DuneGeometry* dune_geometry) {
        m_dune_geometry = dune_geometry;
    }

    template<int codim, typename DuneEntity> friend class ConcreteEntity;

public:
    /** Constructor from a pointer to DuneGeometry */
    explicit ConcreteGeometry(const DuneGeometry* dune_geometry) :
        m_dune_geometry(dune_geometry) {}

    const DuneGeometry& duneGeometry() const {
        return *m_dune_geometry;
    }

    virtual GeometryType type() const {
        return m_dune_geometry->type();
    }

    virtual bool affine() const {
        return m_dune_geometry->affine();
    }

    virtual int n_corners() const {
        return m_dune_geometry->corners();
    }

    virtual void corners(arma::Mat<ctype>& c) const {
        const int cdim = DuneGeometry::dimensionworld;
        const int n = m_dune_geometry->corners();
        c.set_size(cdim, n);

        /** \fixme In future this copying should be optimised away by casting
        appropriate columns of c to Dune field vectors. But this
        can't be done until unit tests are in place. */
        typename DuneGeometry::GlobalCoordinate g;
        for (int j = 0; j < n; ++j) {
            g = m_dune_geometry->corner(j);
            for (int i = 0; i < cdim; ++i)
                c(i,j) = g[i];
        }
    }

    virtual void local2global(const arma::Mat<ctype>& local,
                              arma::Mat<ctype>& global) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if (local.n_rows != mdim)
            throw std::invalid_argument("Geometry::local2global(): invalid dimensions of the 'local' array");
#endif
        const int n = local.n_cols;
        global.set_size(cdim, n);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::GlobalCoordinate g;
        typename DuneGeometry::LocalCoordinate l;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,j);
            g = m_dune_geometry->global(l);
            for (int i = 0; i < cdim; ++i)
                global(i,j) = g[i];
        }
    }

    virtual void global2local(const arma::Mat<ctype>& global,
                              arma::Mat<ctype>& local) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if (global.n_rows != cdim)
            throw std::invalid_argument("Geometry::global2local(): invalid dimensions of the 'global' array");
#endif
        const int n = global.n_cols;
        local.set_size(mdim, n);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::GlobalCoordinate g;
        typename DuneGeometry::LocalCoordinate l;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < cdim; ++i)
                g[i] = global(i,j);
            l = m_dune_geometry->local(l);
            for (int i = 0; i < cdim; ++i)
                local(i,j) = l[i];
        }
    }

    virtual void integrationElement(const arma::Mat<ctype>& local,
                                    arma::Row<ctype>& int_element) const {
        const int mdim = DuneGeometry::mydimension;
#ifndef NDEBUG
        if (local.n_rows != mdim)
            throw std::invalid_argument("Geometry::local2global(): invalid dimensions of the 'local' array");
#endif
        const int n = local.n_cols;
        int_element.set_size(n);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::LocalCoordinate l;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,j);
            ctype ie = m_dune_geometry->integrationElement(l);
            int_element(j) = ie;
        }
    }

    virtual ctype volume() const {
        return m_dune_geometry->volume();
    }

    virtual void center(arma::Col<ctype>& c) const {
        const int cdim = DuneGeometry::coorddimension;
        c.set_size(cdim);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::GlobalCoordinate g = m_dune_geometry->center();
        for (int i = 0; i < cdim; ++i)
            c(i) = g[i];
    }

    virtual void jacobianTransposed(const arma::Mat<ctype>& local,
                                    arma::Cube<ctype>& jacobian_t) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if (local.n_rows != mdim)
            throw std::invalid_argument("Geometry::jacobianTransposed(): invalid dimensions of the 'local' array");
#endif
        const int n = local.n_cols;
        jacobian_t.set_size(mdim, cdim, n);

        /** \bug Unfortunately Dune::FieldMatrix (the underlying type of
        JacobianTransposed) stores elements rowwise, while Armadillo does it
        columnwise. Hence element-by-element filling of jacobian_t seems
        unavoidable). */
        typename DuneGeometry::JacobianTransposed j_t;
        typename DuneGeometry::LocalCoordinate l;
        for (int k = 0; k < n; ++k) {
            /** \fixme However, this bit of data copying could be avoided. */
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,k);
            j_t = m_dune_geometry->jacobianTransposed(l);
            for (int j = 0; j < cdim; ++j)
                for (int i = 0; i < mdim; ++i)
                    jacobian_t(i,j,k) = j_t[i][j];
        }
    }

    virtual void jacobianInverseTransposed(const arma::Mat<ctype>& local,
                                           arma::Cube<ctype>& jacobian_inv_t) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if (local.n_rows != mdim)
            throw std::invalid_argument("Geometry::jacobianInverseTransposed(): invalid dimensions of the 'local' array");
#endif
        const int n = local.n_cols;
        jacobian_inv_t.set_size(cdim, mdim, n);

        /** \bug Unfortunately Dune::FieldMatrix (the underlying type of
        Jacobian) stores elements rowwise, while Armadillo does it
        columnwise. Hence element-by-element filling of jacobian_t seems
        unavoidable). */
        typename DuneGeometry::Jacobian j_inv_t;
        typename DuneGeometry::LocalCoordinate l;
        for (int k = 0; k < n; ++k) {
            /** \fixme However, this bit of data copying could be avoided. */
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,k);
            j_inv_t = m_dune_geometry->jacobianTransposed(l);
            for (int j = 0; j < mdim; ++j)
                for (int i = 0; i < cdim; ++i)
                    jacobian_inv_t(i,j,k) = j_inv_t[i][j];
        }
    }
};

} // namespace Bempp

#endif
