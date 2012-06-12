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

#ifndef fiber_dune_basis_helper_hpp
#define fiber_dune_basis_helper_hpp

#include "../common/common.hpp"

#include "types.hpp"
#include "_4d_array.hpp"

#include "../common/armadillo_fwd.hpp"
#include <vector>

namespace Fiber
{

template <typename CoordinateType, typename ValueType, typename DuneBasis>
void evaluateBasisFunctionsWithDune(
        const arma::Mat<CoordinateType>& local,
        LocalDofIndex localDofIndex,
        arma::Cube<ValueType>& result)
{
    typedef typename DuneBasis::Traits Traits;

    DuneBasis basis;
    const int functionCount = localDofIndex == ALL_DOFS ? basis.size() : 1;
    const int pointCount = local.n_cols;

    typename Traits::DomainType point;
    std::vector<typename Traits::RangeType> values;
    result.set_size(Traits::dimRange, functionCount, pointCount);

    for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex)
    {
        for (int dim = 0; dim < Traits::dimDomain; ++dim)
            point[dim] = local(dim, pointIndex);
        basis.evaluateFunction(point, values);
        if (localDofIndex == ALL_DOFS)
            for (int functionIndex = 0; functionIndex < functionCount; ++functionIndex)
                for (int dim = 0; dim < Traits::dimRange; ++dim)
                    result(dim, functionIndex, pointIndex) = values[functionIndex][dim];
        else
            for (int dim = 0; dim < Traits::dimRange; ++dim)
                result(dim, 0, pointIndex) = values[localDofIndex][dim];
    }
}

template <typename CoordinateType, typename ValueType, typename DuneBasis>
void evaluateBasisFunctionDerivativesWithDune(
        const arma::Mat<CoordinateType>& local,
        LocalDofIndex localDofIndex,
        _4dArray<ValueType>& result)
{
    typedef typename DuneBasis::Traits Traits;

    DuneBasis basis;
    const int functionCount = localDofIndex == ALL_DOFS ? basis.size() : 1;
    const int pointCount = local.n_cols;

    typename Traits::DomainType point;
    std::vector<typename Traits::JacobianType> jacobians;
    result.set_size(Traits::dimRange, Traits::dimDomain, functionCount, pointCount);

    for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex)
    {
        for (int dim = 0; dim < Traits::dimDomain; ++dim)
            point[dim] = local(dim, pointIndex);
        basis.evaluateJacobian(point, jacobians);
        if (localDofIndex == ALL_DOFS)
            for (int functionIndex = 0; functionIndex < functionCount; ++functionIndex)
                for (int dimD = 0; dimD < Traits::dimDomain; ++dimD)
                    for (int dimR = 0; dimR < Traits::dimRange; ++dimR)
                        result(dimR, dimD, functionIndex, pointIndex) =
                                jacobians[functionIndex][dimR][dimD];
        else
            for (int dimD = 0; dimD < Traits::dimDomain; ++dimD)
                for (int dimR = 0; dimR < Traits::dimRange; ++dimR)
                    result(dimR, dimD, 0, pointIndex) =
                            jacobians[localDofIndex][dimR][dimD];
    }
}

}

#endif // DUNE_BASIS_HELPER_HPP
