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

#ifndef bempp_raviart_thomas_0_vector_space_hpp
#define bempp_raviart_thomas_0_vector_space_hpp

#include "../common/common.hpp"

#include "space.hpp"
#include "../grid/grid_view.hpp"
#include "../common/types.hpp"
#include "../fiber/raviart_thomas_order_0_basis.hpp"

#include <boost/scoped_ptr.hpp>
#include <map>
#include <memory>

namespace Bempp
{

/** \cond FORWARD_DECL */
class GridView;
/** \endcond */

/** \ingroup space
 *  \brief Space of continuous, piecewise linear scalar functions. */
template <typename BasisFunctionType>
class RaviartThomas0VectorSpace : public Space<BasisFunctionType>
{
    typedef Space<BasisFunctionType> Base;
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;
    typedef typename Space<BasisFunctionType>::ComplexType ComplexType;
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;

    explicit RaviartThomas0VectorSpace(
        const shared_ptr<const Grid>& grid);
    virtual ~RaviartThomas0VectorSpace();

    virtual const CollectionOfBasisTransformations& shapeFunctionValue() const;


    virtual int domainDimension() const;
    virtual int codomainDimension() const;

    /** \brief Return the variant of element \p element.
     *
     *  Possible return values:
     *    - 3: triangular element,
     *    - 4: quadrilateral element. */
    virtual ElementVariant elementVariant(const Entity<0>& element) const;
    virtual void setElementVariant(const Entity<0>& element,
                                   ElementVariant variant);

    virtual const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;

    virtual size_t globalDofCount() const;
    virtual size_t flatLocalDofCount() const;
    virtual void getGlobalDofs(const Entity<0>& element,
                               std::vector<GlobalDofIndex>& dofs,
                               std::vector<BasisFunctionType>& dofWeights) const;
    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs,
            std::vector<std::vector<BasisFunctionType> >& localDofWeights) const;
    virtual void flatLocal2localDofs(
            const std::vector<FlatLocalDofIndex>& flatLocalDofs,
            std::vector<LocalDof>& localDofs) const;

    virtual void getGlobalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const;
    virtual void getFlatLocalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const;
    virtual void dumpClusterIds(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs) const;
    virtual void dumpClusterIdsEx(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
            DofType dofType) const;

private:
    void assignDofsImpl();

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
    std::auto_ptr<GridView> m_view;
    Fiber::RaviartThomasOrder0Basis<3, BasisFunctionType> m_triangleBasis;
    std::vector<std::vector<GlobalDofIndex> > m_local2globalDofs;
    std::vector<std::vector<BasisFunctionType> > m_local2globalDofWeights;
    std::vector<std::vector<LocalDof> > m_global2localDofs;
    std::vector<std::vector<BasisFunctionType> > m_global2localDofWeights;
    std::vector<LocalDof> m_flatLocal2localDofs;
    size_t m_flatLocalDofCount;
};

} // namespace Bempp

#endif
