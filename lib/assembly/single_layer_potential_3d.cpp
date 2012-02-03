#include "single_layer_potential_3d.hpp"

namespace Bempp
{

#ifdef COMPILE_FOR_FLOAT
template class SingleLayerPotential3D<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class SingleLayerPotential3D<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class SingleLayerPotential3D<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class SingleLayerPotential3D<std::complex<double> >;
#endif

} // namespace Bempp
