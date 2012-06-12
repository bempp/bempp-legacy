#ifndef bempp_mass_matrix_container_hpp
#define bempp_mass_matrix_container_hpp

#include "../common/common.hpp"

#include <memory>

namespace Bempp
{

template <typename ValueType> class DiscreteLinearOperator;

template <typename BasisFunctionType>
struct MassMatrixContainer
{
    std::auto_ptr<DiscreteLinearOperator<BasisFunctionType> > massMatrix;
    std::auto_ptr<DiscreteLinearOperator<BasisFunctionType> > inverseMassMatrix;
};

} //namespace Bempp

#endif
