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
    virtual ~Basis() {}

    virtual int size() const = 0;
    /** \brief Maximum polynomial order of basis elements. */
    virtual int order() const = 0;
    virtual void evaluate(int what,
                          const arma::Mat<ValueType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const = 0;

    /**
     * \brief Returns an OpenCL code snippet for basis function evaluation
     * \note The code snippet must provide device function devBasisEval
     */
    virtual std::string clCodeString (const std::string modifier) const {
        throw std::runtime_error("Basis: clCodeString not implemented yet");
    }
};

} // namespace Fiber

#endif
