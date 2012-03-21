#include "hypersingular_operator_3d.hpp"

namespace Bempp
{

#ifdef COMPILE_FOR_FLOAT
template class HypersingularOperator3D<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class HypersingularOperator3D<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class HypersingularOperator3D<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class HypersingularOperator3D<std::complex<double> >;
#endif

} // namespace Bempp
