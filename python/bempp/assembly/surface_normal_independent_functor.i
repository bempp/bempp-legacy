%{
#include "assembly/surface_normal_independent_functor.hpp"
  %}




%define SURFACE_NORMAL_INDEPENDENT_FUNCTOR(TYPE,NPYNAME)

%include "assembly/surface_normal_independent_functor.hpp"

namespace Bempp{
  %template(SurfaceNormalIndependentFunctor_## NPYNAME) SurfaceNormalIndependentFunctor< TYPE >;
    }

%enddef // SURFACE_NORMAL_INDEPENDENT_FUNCTOR_REAL


SURFACE_NORMAL_INDEPENDENT_FUNCTOR(float,float32)
SURFACE_NORMAL_INDEPENDENT_FUNCTOR(double,float64)
SURFACE_NORMAL_INDEPENDENT_FUNCTOR(std::complex<float>,complex32)
SURFACE_NORMAL_INDEPENDENT_FUNCTOR(std::complex<double>,complex64)
