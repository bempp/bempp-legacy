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

#ifndef bempp_piecewise_linear_discontinuous_scalar_space_hpp
#define bempp_piecewise_linear_discontinuous_scalar_space_hpp

#include "../common/common.hpp"

#include "../grid/grid_view.hpp"
#include "piecewise_linear_scalar_space.hpp"
#include "../common/types.hpp"
// The name is absurd. Change to linear_scalar_basis.hpp
#include "../fiber/piecewise_linear_continuous_scalar_basis.hpp"

#include <map>
#include <memory>
#include <tbb/mutex.h>

namespace Bempp
{

/** \cond FORWARD_DECL */
class GridSegment;
class GridView;
/** \endcond */

/** \ingroup space
 *  \brief Space of piecewise linear, not necessarily continuous, scalar
 *  functions. */
template <typename BasisFunctionType>
class PiecewiseLinearDiscontinuousScalarSpace : public PiecewiseLinearScalarSpace<BasisFunctionType>
{
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;
    typedef typename Space<BasisFunctionType>::ComplexType ComplexType;

    /** \brief Constructor.
     *
     *  Construct a space of piecewise linear, not necessarily continuous,
     *  scalar functions defined on the grid \p grid with level \p level.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    explicit PiecewiseLinearDiscontinuousScalarSpace(
            const shared_ptr<const Grid>& grid, unsigned int level=0);

    /** \brief Constructor.
     *
     *  Construct a space of piecewise linear, not necessarily continuous,
     *  scalar functions defined on the segment \p segment of the grid \p grid.
     *  If \p strictlyOnSegment is set to \c false (default), the space will
     *  include all basis functions associated with vertices belonging to \p
     *  segment, regardless of whether the elements on which these functions
     *  are defined belong themselves to \p segment. In consequence, the
     *  resulting space will be (in the mathematical sense) a superset of a
     *  PiecewiseLinearContinuousScalarSpace defined on the same segment. If \p
     *  strictlyOnSegment is set to \c true, the space will only include basis
     *  functions defined on elements belonging to \p segment.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    PiecewiseLinearDiscontinuousScalarSpace(
            const shared_ptr<const Grid>& grid,
            const GridSegment& segment,
            bool strictlyOnSegment = false,
            unsigned int level=0);
    virtual ~PiecewiseLinearDiscontinuousScalarSpace();

    virtual shared_ptr<const Space<BasisFunctionType> > discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType> >& self) const;
    virtual bool isDiscontinuous() const;

    virtual bool isBarycentric() const {
        return false;
    }

    virtual size_t globalDofCount() const;
    virtual size_t flatLocalDofCount() const;
    virtual void getGlobalDofs(const Entity<0>& element,
                            std::vector<GlobalDofIndex>& dofs) const;
    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs) const;
    virtual void flatLocal2localDofs(
            const std::vector<FlatLocalDofIndex>& flatLocalDofs,
            std::vector<LocalDof>& localDofs) const;

    virtual void getGlobalDofPositions(std::vector<Point3D<CoordinateType> >& positions) const;
    virtual void getFlatLocalDofPositions(std::vector<Point3D<CoordinateType> >& positions) const;

    virtual void getGlobalDofBoundingBoxes(
            std::vector<BoundingBox<CoordinateType> >& bboxes) const;
    virtual void getFlatLocalDofBoundingBoxes(
            std::vector<BoundingBox<CoordinateType> >& bboxes) const;

    virtual void getGlobalDofNormals(std::vector<Point3D<CoordinateType> >& normals) const;
    virtual void getFlatLocalDofNormals(std::vector<Point3D<CoordinateType> >& normals) const;

    virtual void dumpClusterIds(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs) const;
    virtual void dumpClusterIdsEx(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
            DofType dofType) const;

private:
    void initialize(const GridSegment& segment, bool strictlyOnSegment = false);
    void assignDofsImpl(const GridSegment& segment,
                        bool strictlyOnSegment = false);

private:
    /** \cond PRIVATE */
    std::auto_ptr<GridView> m_view;
    Fiber::PiecewiseLinearContinuousScalarBasis<2, BasisFunctionType> m_lineBasis;
    Fiber::PiecewiseLinearContinuousScalarBasis<3, BasisFunctionType> m_triangleBasis;
    Fiber::PiecewiseLinearContinuousScalarBasis<4, BasisFunctionType> m_quadrilateralBasis;
    std::vector<std::vector<GlobalDofIndex> > m_local2globalDofs;
    std::vector<std::vector<LocalDof> > m_global2localDofs;
    std::vector<LocalDof> m_flatLocal2localDofs;
    /** \endcond */
};

} // namespace Bempp

#endif
