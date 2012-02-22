#ifndef fiber_basis_hpp
#define fiber_basis_hpp

#include "types.hpp"

namespace Fiber
{

template <typename ValueType> struct BasisData;

template <typename ValueType>
class Basis
{
public:
    virtual int size() const = 0;
    /** \brief Maximum polynomial order of basis elements. */
    virtual int order() const = 0;
    virtual void evaluate(int what,
                          const arma::Mat<ValueType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const = 0;
};

} // namespace Fiber

#endif
