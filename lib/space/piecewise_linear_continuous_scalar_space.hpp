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

#ifndef bempp_piecewise_linear_continuous_scalar_space_hpp
#define bempp_piecewise_linear_continuous_scalar_space_hpp

#include "space.hpp"
#include "../common/types.hpp"
#include "../fiber/piecewise_linear_continuous_scalar_basis.hpp"

#include <map>
#include <memory>

namespace Bempp
{

class GridView;

// Element variants: 2 (line element), 3 (triangular element), 4 (quadrilateral element)

template <typename ValueType>
class PiecewiseLinearContinuousScalarSpace : public Space<ValueType>
{
public:
    PiecewiseLinearContinuousScalarSpace(Grid& grid);

    virtual int domainDimension() const;
    virtual int codomainDimension() const;

    virtual ElementVariant elementVariant(const Entity<0>& element) const;
    virtual void setElementVariant(const Entity<0>& element,
                                   ElementVariant variant);

    virtual void getBases(const std::vector<const EntityPointer<0>*>& elements,
                          std::vector<const Fiber::Basis<ValueType>*>& bases) const;

    virtual const Fiber::Basis<ValueType>& basis(const EntityPointer<0>& element) const;

    virtual void assignDofs();
    virtual bool dofsAssigned() const;
    virtual int globalDofCount() const;
    virtual void globalDofs(const Entity<0>& element,
                            std::vector<GlobalDofIndex>& dofs) const;    
    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs) const;

    virtual void globalDofPositions(std::vector<Point3D>& positions) const;

private:
    std::auto_ptr<GridView> m_view;
    Fiber::PiecewiseLinearContinuousScalarBasis<2, ValueType> m_lineBasis;
    Fiber::PiecewiseLinearContinuousScalarBasis<3, ValueType> m_triangleBasis;
    Fiber::PiecewiseLinearContinuousScalarBasis<4, ValueType> m_quadrilateralBasis;
    typedef std::map<EntityIndex, std::vector<GlobalDofIndex> > GlobalDofMap;
    GlobalDofMap m_local2globalDofs;
    std::vector<std::vector<LocalDof> > m_global2localDofs;
};

} // namespace Bempp

#endif
