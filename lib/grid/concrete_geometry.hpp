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

#ifndef bempp_concrete_geometry_hpp
#define bempp_concrete_geometry_hpp

#include "../common/common.hpp"

#include "geometry.hpp"
#include "dune.hpp"
#include "geometry_type.hpp"

#include "../common/not_implemented_error.hpp"
#include "../fiber/geometrical_data.hpp"
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/static_assert.hh>
#include <dune/grid/common/grid.hh>

#include "../common/armadillo_fwd.hpp"
#include <memory>

namespace Bempp
{

/** \brief Wrapper of a Dune geometry of type \p DuneGeometry */

template <typename DuneGeometry> class ConcreteGeometryFactory;

// Internal helper function for setting up geometries manually (independently
// from the grid).

template <typename DuneGeometry>
void setupDuneGeometry(
        typename Dune::MakeableInterfaceObject<DuneGeometry>::ImplementationType& duneGeometry,
        const GeometryType& type,
        const std::vector<Dune::FieldVector<double, DuneGeometry::dimensionworld> >& corners)
{
    throw std::runtime_error("setupDuneGeometry() is only supported for a "
                             "small number of Dune geometry classes");
}

template <>
void setupDuneGeometry<Default2dIn3dDuneGrid::Codim<0>::Geometry>(
        Dune::MakeableInterfaceObject<Default2dIn3dDuneGrid::Codim<0>::Geometry>::ImplementationType& duneGeometry,
        const GeometryType& type,
        const std::vector<Dune::FieldVector<double, Default2dIn3dDuneGrid::dimensionworld> >& corners)
{
    duneGeometry.setup(type, corners);
}

template <typename DuneGeometry>
class ConcreteGeometry : public Geometry
{
    dune_static_assert((int)DuneGeometry::coorddimension ==
                       (int)DuneGeometry::dimensionworld,
                       "ConcreteGeometry: world dimension does not agree with "
                       "number of coordinates");

private:
    const DuneGeometry* m_dune_geometry;
    bool m_owns_dune_geometry;

    /** \brief Default constructor.

    Should be followed by a call to setDuneGeometry(). */
    ConcreteGeometry() : m_dune_geometry(0), m_owns_dune_geometry(false) {
    }

    void setDuneGeometry(const DuneGeometry* dune_geometry, bool owns) {
        if (m_owns_dune_geometry)
            delete m_dune_geometry;
        m_dune_geometry = dune_geometry;
        m_owns_dune_geometry = owns;
    }

    template<int codim, typename DuneEntity> friend class ConcreteEntity;
    friend class ConcreteGeometryFactory<DuneGeometry>;

public:
    /** \brief Constructor from a pointer to DuneGeometry.

      This object does not assume ownership of \p *dune_geometry.
    */
    explicit ConcreteGeometry(const DuneGeometry* dune_geometry) :
        m_dune_geometry(dune_geometry) {}

    ~ConcreteGeometry() {
        if (m_owns_dune_geometry)
            delete m_dune_geometry;
    }

    /** \brief Read-only access to the underlying Dune geometry object. */
    const DuneGeometry& duneGeometry() const {
        return *m_dune_geometry;
    }

    virtual int dim() const {
        return DuneGeometry::mydimension;
    }

    virtual int dimWorld() const {
        return DuneGeometry::coorddimension;
    }

    virtual void setupImpl(const arma::Mat<double>& corners,
                           const arma::Col<char>& auxData) {
        const int dimWorld = DuneGeometry::coorddimension;
        const int cornerCount = corners.n_cols;
        assert((int)corners.n_rows == dimWorld);

        GeometryType type;
        if (DuneGeometry::mydimension == 0) {
            assert(cornerCount == 1);
            type.makeVertex();
        }
        else if (DuneGeometry::mydimension == 1) {
            assert(cornerCount == 2);
            type.makeLine();
        }
        else if (DuneGeometry::mydimension == 2) {
            assert (cornerCount == 3 || cornerCount == 4);
            if (cornerCount == 3)
                type.makeTriangle();
            else
                type.makeQuadrilateral();
        }
        else
            throw NotImplementedError("ConcreteGeometry::setup(): "
                                      "not implemented yet for 3D entities");

        std::vector<Dune::FieldVector<double, dimWorld> > duneCorners(cornerCount);
        for (size_t i = 0; i < corners.n_cols; ++i)
            for (int j = 0; j < dimWorld; ++j)
                duneCorners[i][j] = corners(j, i);

        typedef Dune::MakeableInterfaceObject<DuneGeometry> DuneMakeableGeometry;
        typedef typename DuneMakeableGeometry::ImplementationType DuneGeometryImp;

        DuneGeometryImp newDuneGeometry;
        setupDuneGeometry<DuneGeometry>(newDuneGeometry, type, duneCorners);
        // newDuneGeometry.setup(type, duneCorners);
        // I wish I could avoid this heap allocation...
        // unfortunately I can't find a way to obtain access from
        // Dune's Geometry<... GeometryImp> to the underlying GeometryImp object
        setDuneGeometry(new DuneMakeableGeometry(newDuneGeometry), true /* owns */);
    }

    virtual GeometryType type() const {
        return m_dune_geometry->type();
    }

    virtual bool affine() const {
        return m_dune_geometry->affine();
    }

    virtual int cornerCount() const {
        return m_dune_geometry->corners();
    }

    virtual void getCornersImpl(arma::Mat<double>& c) const {
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

    virtual void local2globalImpl(const arma::Mat<double>& local,
                                  arma::Mat<double>& global) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if ((int)local.n_rows != mdim)
            throw std::invalid_argument("Geometry::local2global(): invalid dimensions of the 'local' array");
#endif
        const size_t n = local.n_cols;
        global.set_size(cdim, n);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::GlobalCoordinate g;
        typename DuneGeometry::LocalCoordinate l;
        for (size_t j = 0; j < n; ++j) {
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,j);
            g = m_dune_geometry->global(l);
            for (int i = 0; i < cdim; ++i)
                global(i,j) = g[i];
        }
    }

    virtual void global2localImpl(const arma::Mat<double>& global,
                                  arma::Mat<double>& local) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if ((int)global.n_rows != cdim)
            throw std::invalid_argument("Geometry::global2local(): invalid dimensions of the 'global' array");
#endif
        const size_t n = global.n_cols;
        local.set_size(mdim, n);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::GlobalCoordinate g;
        typename DuneGeometry::LocalCoordinate l;
        for (size_t j = 0; j < n; ++j) {
            for (int i = 0; i < cdim; ++i)
                g[i] = global(i,j);
            l = m_dune_geometry->local(g);
            for (int i = 0; i < mdim; ++i)
                local(i,j) = l[i];
        }
    }

    virtual void getIntegrationElementsImpl(const arma::Mat<double>& local,
                                            arma::Row<double>& int_element) const {
        const int mdim = DuneGeometry::mydimension;
#ifndef NDEBUG
        if ((int)local.n_rows != mdim)
            throw std::invalid_argument("Geometry::local2global(): invalid dimensions of the 'local' array");
#endif
        const size_t n = local.n_cols;
        int_element.set_size(n);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::LocalCoordinate l;
        for (size_t j = 0; j < n; ++j) {
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,j);
            double ie = m_dune_geometry->integrationElement(l);
            int_element(j) = ie;
        }
    }

    virtual double volume() const {
        return m_dune_geometry->volume();
    }

    virtual void getCenterImpl(arma::Col<double>& c) const {
        const int cdim = DuneGeometry::coorddimension;
        c.set_size(cdim);

        /** \fixme Optimise (get rid of data copying). */
        typename DuneGeometry::GlobalCoordinate g = m_dune_geometry->center();
        for (int i = 0; i < cdim; ++i)
            c(i) = g[i];
    }

    virtual void getJacobiansTransposedImpl(const arma::Mat<double>& local,
                                            arma::Cube<double>& jacobian_t) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if ((int)local.n_rows != mdim)
            throw std::invalid_argument("Geometry::getJacobiansTransposed(): "
                                        "invalid dimensions of the 'local' array");
#endif
        const size_t n = local.n_cols;
        jacobian_t.set_size(mdim, cdim, n);

        /** \bug Unfortunately Dune::FieldMatrix (the underlying type of
        JacobianTransposed) stores elements rowwise, while Armadillo does it
        columnwise. Hence element-by-element filling of jacobian_t seems
        unavoidable). */
        typename DuneGeometry::JacobianTransposed j_t;
        typename DuneGeometry::LocalCoordinate l;
        for (size_t k = 0; k < n; ++k) {
            /** \fixme However, this bit of data copying could be avoided. */
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,k);
            j_t = m_dune_geometry->jacobianTransposed(l);
            for (int j = 0; j < cdim; ++j)
                for (int i = 0; i < mdim; ++i)
                    jacobian_t(i,j,k) = j_t[i][j];
        }
    }

    virtual void getJacobianInversesTransposedImpl(
            const arma::Mat<double>& local,
            arma::Cube<double>& jacobian_inv_t) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;
#ifndef NDEBUG
        if ((int)local.n_rows != mdim)
            throw std::invalid_argument("Geometry::getJacobianInversesTransposed(): "
                                        "invalid dimensions of the 'local' array");
#endif
        const size_t n = local.n_cols;
        jacobian_inv_t.set_size(cdim, mdim, n);

        /** \bug Unfortunately Dune::FieldMatrix (the underlying type of
        Jacobian) stores elements rowwise, while Armadillo does it
        columnwise. Hence element-by-element filling of jacobian_t seems
        unavoidable). */
        typename DuneGeometry::Jacobian j_inv_t;
        typename DuneGeometry::LocalCoordinate l;
        for (size_t k = 0; k < n; ++k) {
            /** \fixme However, this bit of data copying could be avoided. */
            for (int i = 0; i < mdim; ++i)
                l[i] = local(i,k);
            j_inv_t = m_dune_geometry->jacobianInverseTransposed(l);
            for (int j = 0; j < mdim; ++j)
                for (int i = 0; i < cdim; ++i)
                    jacobian_inv_t(i,j,k) = j_inv_t[i][j];
        }
    }

    virtual void getNormalsImpl(const arma::Mat<double>& local,
                                arma::Mat<double>& normal) const {
        arma::Cube<double> jacobian_t;
        getJacobiansTransposed(local, jacobian_t);
        calculateNormals(jacobian_t, normal);
    }

    virtual void getDataImpl(size_t what, const arma::Mat<double>& local,
                             Fiber::GeometricalData<double>& data) const {
        // In this first implementation we call the above virtual functions as required.
        // In future some optimisations (elimination of redundant calculations)
        // might be possible.

        typedef ConcreteGeometry<DuneGeometry> This; // to avoid virtual function calls

        if (what & Fiber::GLOBALS)
            This::local2global(local, data.globals);
        if (what & Fiber::INTEGRATION_ELEMENTS)
            This::getIntegrationElements(local, data.integrationElements);
        if (what & Fiber::JACOBIANS_TRANSPOSED || what & Fiber::NORMALS)
            This::getJacobiansTransposed(local, data.jacobiansTransposed);
        if (what & Fiber::JACOBIAN_INVERSES_TRANSPOSED)
            This::getJacobianInversesTransposed(
                        local, data.jacobianInversesTransposed);
        if (what & Fiber::NORMALS)
            calculateNormals(data.jacobiansTransposed, data.normals);
    }

private:
    void calculateNormals(const arma::Cube<double>& jt,
                          arma::Mat<double>& normals) const {
        const int mdim = DuneGeometry::mydimension;
        const int cdim = DuneGeometry::coorddimension;

        if (mdim != cdim - 1)
            throw std::logic_error("ConcreteGeometry::calculateNormals(): "
                                   "normal vectors are defined only for "
                                   "entities of dimension (worldDimension - 1)");

        const size_t pointCount = jt.n_slices;
        normals.set_size(cdim, pointCount);

        // First calculate normal vectors of arbitrary length

        // Compile-time if
        if (cdim == 3)
            for (size_t i = 0; i < pointCount; ++i) {
                normals(0,i) = jt(0,1,i) * jt(1,2,i) - jt(0,2,i) * jt(1,1,i);
                normals(1,i) = jt(0,2,i) * jt(1,0,i) - jt(0,0,i) * jt(1,2,i);
                normals(2,i) = jt(0,0,i) * jt(1,1,i) - jt(0,1,i) * jt(1,0,i);
            }
        else if (cdim == 2)
            for (size_t i = 0; i < pointCount; ++i) {
                normals(0,i) = jt(0,1,i);
                normals(1,i) = jt(0,0,i);
            }
        else if (cdim == 1) // probably unnecessary
            for (size_t i = 0; i < pointCount; ++i)
                normals(0,i) = 1.;
        else
            throw std::runtime_error("ConcreteGeometry::calculateNormals(): "
                                     "Normal vector is not defined for "
                                     "zero-dimensional space");

        // Now set vector length to 1.

        for (size_t i = 0; i < pointCount; ++i) {
            double sum = 0.;
            for (int dim = 0; dim < cdim; ++dim)
                sum += normals(dim, i) * normals(dim, i);
            normals.col(i) /= sqrt(sum);
        }
    }
};

} // namespace Bempp

#endif
