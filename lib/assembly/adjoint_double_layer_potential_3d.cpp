#include "adjoint_double_layer_potential_3d.hpp"

namespace Bempp
{

#ifdef COMPILE_FOR_FLOAT
template class AdjointDoubleLayerPotential3D<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class AdjointDoubleLayerPotential3D<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class AdjointDoubleLayerPotential3D<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class AdjointDoubleLayerPotential3D<std::complex<double> >;
#endif

} // namespace Bempp
