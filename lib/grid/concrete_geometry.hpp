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

#ifndef bempp_concrete_geometry_hpp
#define bempp_concrete_geometry_hpp

#include "geometry.hpp"
#include "geometry_type.hpp"
#include "common.hpp"

#include "../common/not_implemented_error.hpp"
#include "../fiber/geometrical_data.hpp"

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/static_assert.hh>
#include <dune/grid/common/grid.hh>

#include <armadillo>
#include <memory>

namespace Bempp
{

/** \brief Wrapper of a Dune geometry of type \p DuneGeometry */

template <typename DuneGeometry> class ConcreteGeometryFactory;

template<typename DuneGeometry>
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

    virtual void setup(const arma::Mat<ctype>& corners,
                       const arma::Col<char>& auxData) {
        const int dimWorld = DuneGeometry::coorddimension;
        const int cornerCount = corners.n_cols;
        assert(corners.n_rows == dimWorld);

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

        std::vector<Dune::FieldVector<ctype, dimWorld> > duneCorners(cornerCount);
        for (int i = 0; i < corners.n_cols; ++i)
            for (int j = 0; j < dimWorld; ++j)
                duneCorners[i][j] = corners(j, i);

        typedef Dune::MakeableInterfaceObject<DuneGeometry> DuneMakeableGeometry;
        typedef typename DuneMakeableGeometry::ImplementationType DuneGeometryImp;

        DuneGeometryImp newDuneGeometry;
        newDuneGeometry.setup(type, duneCorners);
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
            l = m_dune_geometry->local(g);
            for (int i = 0; i < mdim; ++i)
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
            j_inv_t = m_dune_geometry->jacobianInverseTransposed(l);
            for (int j = 0; j < mdim; ++j)
                for (int i = 0; i < cdim; ++i)
                    jacobian_inv_t(i,j,k) = j_inv_t[i][j];
        }
    }

    virtual void getData(int what, const arma::Mat<ctype>& local,
                         Fiber::GeometricalData<ctype>& data) const {
        // In this first implementation we call the above virtual functions as required.
        // In future some optimisations (elimination of redundant calculations)
        // might be possible.

        typedef ConcreteGeometry<DuneGeometry> This; // to avoid virtual function calls

        if (what & Fiber::GLOBALS)
            This::local2global(local, data.globals);
        if (what & Fiber::INTEGRATION_ELEMENTS)
            This::integrationElement(local, data.integrationElements);
        if (what & Fiber::JACOBIANS_TRANSPOSED)
            This::jacobianTransposed(local, data.jacobiansTransposed);
        if (what & Fiber::JACOBIAN_INVERSES_TRANSPOSED)
            This::jacobianInverseTransposed(local, data.jacobianInversesTransposed);
        if (what & Fiber::NORMALS)
            throw std::runtime_error("Geometry::getData(): calculation of"
                                     "normals not implemented yet");
    }

};

} // namespace Bempp

#endif
