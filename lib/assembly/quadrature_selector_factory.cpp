#include "quadrature_selector_factory.hpp"

#include "../common/not_implemented_error.hpp"

namespace Bempp
{

template <typename ValueType>
std::auto_ptr<QuadratureSelector<ValueType> >
QuadratureSelectorFactory<ValueType>::make(
        const AssemblyOptions& options) {
    throw NotImplementedError("QuadratureSelectorFactory::make(): "
                              "not implemented yet");
}

#ifdef COMPILE_FOR_FLOAT
template class QuadratureSelectorFactory<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class QuadratureSelectorFactory<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class QuadratureSelectorFactory<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class QuadratureSelectorFactory<std::complex<double> >;
#endif

} // namespace Bempp
