#ifndef fiber_piecewise_constant_scalar_basis_hpp
#define fiber_piecewise_constant_scalar_basis_hpp

#include "basis.hpp"
#include "basis_data.hpp"

#include <algorithm>

namespace Fiber
{

template <typename ValueType>
class PiecewiseConstantScalarBasis : public Basis<ValueType>
{
public:
    virtual int size() const {
        return 1;
    }

    virtual int order() const {
        return 0;
    }

    virtual void evaluate(int what,
                          const arma::Mat<ValueType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const {
        if (localDofIndex != ALL_DOFS && localDofIndex != 0)
            throw std::invalid_argument("PiecewiseConstantScalarBasis::evaluate(): "
                                        "Invalid localDofIndex");
        // Since there is only one basis function, there is no difference
        // between calculating all basis functions and just one.

        const int componentCount = 1;
        const int functionCount = 1;
        const int pointCount = points.n_cols;
        if (what & VALUES)
        {
            data.values.set_size(componentCount, functionCount, pointCount);
            data.values.fill(1.);
        }
        if (what & DERIVATIVES)
        {
            const int coordCount = points.n_rows;
            data.derivatives.set_size(componentCount, coordCount,
                                 functionCount, pointCount);
            std::fill(data.derivatives.begin(), data.derivatives.end(), 0.);
        }
    }
};

} // namespace Fiber

#endif
