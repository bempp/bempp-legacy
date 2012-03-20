#include "double_layer_potential_3d.hpp"

namespace Bempp
{

#ifdef COMPILE_FOR_FLOAT
template class DoubleLayerPotential3D<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class DoubleLayerPotential3D<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class DoubleLayerPotential3D<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class DoubleLayerPotential3D<std::complex<double> >;
#endif

} // namespace Bempp
