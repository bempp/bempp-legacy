#ifndef fiber_basis_data_hpp
#define fiber_basis_data_hpp

#include <armadillo>
#include "../common/multidimensional_arrays.hpp"

namespace Fiber
{

enum BasisDataType
{
    VALUES = 0x0001,
    DERIVATIVES = 0x0002
};

template <typename ValueType>
struct BasisData
{
    // values(i,j,l) = (f_k)_i(x_l) ->
    // ith component of kth basis function at lth point
    arma::Cube<ValueType> values;
    // derivatives(i,j,k,l) = (d_j (f_k)_i)(x_l) ->
    // derivative in direction j of ith component of kth basis function at lth point
    Array4D<ValueType> derivatives;
};

} // namespace Fiber

#endif // BASIS_DATA_TYPES_HPP
