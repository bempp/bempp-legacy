%{
#include "assembly/surface_normal_dependent_functor.hpp"
  %}




%define SURFACE_NORMAL_DEPENDENT_FUNCTOR(TYPE,NPYNAME)

namespace Bempp
{
  %extend  SurfaceNormalDependentFunctor< TYPE >
  {
    %ignore evaluate;
  }
}

%include "assembly/surface_normal_dependent_functor.hpp"

namespace Bempp{
  %template(SurfaceNormalDependentFunctor_## NPYNAME) SurfaceNormalDependentFunctor< TYPE >;
    }

%enddef // SURFACE_NORMAL_DEPENDENT_FUNCTOR_REAL


SURFACE_NORMAL_DEPENDENT_FUNCTOR(float,float32)
SURFACE_NORMAL_DEPENDENT_FUNCTOR(double,float64)
SURFACE_NORMAL_DEPENDENT_FUNCTOR(std::complex<float>,complex32)
SURFACE_NORMAL_DEPENDENT_FUNCTOR(std::complex<double>,complex64)









