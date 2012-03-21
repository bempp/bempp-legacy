#ifndef fiber_dune_basis_helper_hpp
#define fiber_dune_basis_helper_hpp

#include "types.hpp"
#include "array_4d.hpp"

#include <armadillo>
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
        Array4D<ValueType>& result)
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
