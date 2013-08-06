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

#ifndef bempp_piecewise_polynomial_discontinuous_scalar_space_hpp
#define bempp_piecewise_polynomial_discontinuous_scalar_space_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../fiber/lagrange_scalar_basis.hpp"

#include "scalar_space.hpp"
#include "dof_assignment_mode.hpp"

#include <map>
#include <memory>
#include <tbb/mutex.h>

namespace Bempp
{

/** \cond FORWARD_DECL */
class GridSegment;
class GridView;
template <typename CoordinateType> class BoundingBox;
/** \endcond */

/** \ingroup space
 *  \brief Space of discontinuous, piecewise polynomial scalar functions. */
template <typename BasisFunctionType>
class PiecewisePolynomialDiscontinuousScalarSpace :
        public ScalarSpace<BasisFunctionType>
{
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;
    typedef typename Space<BasisFunctionType>::ComplexType ComplexType;


    /** \brief Constructor.
     *
     *  Construct a space of functions whose restrictions to
     *  elements of the grid \p grid with level \p level will be polynomials of order at most \p
     *  polynomialOrder. */
    PiecewisePolynomialDiscontinuousScalarSpace(
            const shared_ptr<const Grid>& grid, int polynomialOrder, unsigned int level=0);
    /** \brief Constructor.
     *
     *  Construct a space of continuous functions whose restrictions to
     *  elements of the grid \p grid will be polynomials of order at most \p
     *  polynomialOrder. The space will contain only the basis functions deemed
     *  to belong to the segment \p segment. The precise way in which this is
     *  determined is controlled by the parameter \p dofMode.
     *
     *  \p dofMode can be set to REFERENCE_POINT_ON_SEGMENT,
     *  ELEMENT_ON_SEGMENT, or their (OR-ed) combination. If \p dofMode
     *  contains REFERENCE_POINT_ON_SEGMENT, the constructed space will include
     *  only the vertex basis functions associated with vertices belonging to
     *  \p segment, edge functions associated with edges belonging to \p
     *  segment and bubble function associated with elements belonging to \p
     *  segment. Note that, as a result, the space will be (in the mathematical
     *  sense) a superset of a PiecewisePolynomialContinuousScalarSpace, of the
     *  same order, defined on the same segment. If \p dofMode contains
     *  ELEMENT_ON_SEGMENT, the space will contain only the basis function
     *  whose supports are the elements belonging to \p segment.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    PiecewisePolynomialDiscontinuousScalarSpace(
            const shared_ptr<const Grid>& grid,
            int polynomialOrder,
            const GridSegment& segment,
            int dofMode = REFERENCE_POINT_ON_SEGMENT,
            unsigned int level=0);
    virtual ~PiecewisePolynomialDiscontinuousScalarSpace();

    virtual int domainDimension() const;
    virtual int codomainDimension() const;

    virtual bool isBarycentric() const {
        return false;
    }


    /** \brief Return the variant of element \p element.
     *
     *  Possible return values:
     *    - 2: one-dimensional segment,
     *    - 3: triangular element,
     *    - 4: quadrilateral element. */
    virtual ElementVariant elementVariant(const Entity<0>& element) const;
    virtual void setElementVariant(const Entity<0>& element,
                                   ElementVariant variant);

    virtual const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;

    virtual shared_ptr<const Space<BasisFunctionType> > discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType> >& self) const;
    virtual bool isDiscontinuous() const;

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

    virtual void getGlobalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const;
    virtual void getFlatLocalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const;

    virtual void getGlobalDofBoundingBoxes(
            std::vector<BoundingBox<CoordinateType> >& bboxes) const;
    virtual void getFlatLocalDofBoundingBoxes(
            std::vector<BoundingBox<CoordinateType> >& bboxes) const;

    virtual void getGlobalDofNormals(
            std::vector<Point3D<CoordinateType> >& normals) const;
    virtual void getFlatLocalDofNormals(
            std::vector<Point3D<CoordinateType> >& normals) const;

    virtual void dumpClusterIds(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs) const;
    virtual void dumpClusterIdsEx(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
            DofType dofType) const;

private:
    void initialize(const GridSegment& segment,
                    int dofMode = REFERENCE_POINT_ON_SEGMENT);
    void assignDofsImpl(const GridSegment& segment,
                        int dofMode = REFERENCE_POINT_ON_SEGMENT);

private:
    /** \cond PRIVATE */
    int m_polynomialOrder;
    boost::scoped_ptr<Fiber::Basis<BasisFunctionType> > m_triangleBasis;
    std::auto_ptr<GridView> m_view;
    std::vector<std::vector<GlobalDofIndex> > m_local2globalDofs;
    std::vector<std::vector<LocalDof> > m_global2localDofs;
    std::vector<LocalDof> m_flatLocal2localDofs;
    size_t m_flatLocalDofCount;
    std::vector<BoundingBox<CoordinateType> > m_globalDofBoundingBoxes;
    /** \endcond */
};

} // namespace Bempp

#endif
