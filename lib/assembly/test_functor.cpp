
#include "test_functor.hpp"

namespace Bempp
{

  int testFunctor(const SurfaceNormalIndependentFunctor<double>& functor){
    return functor.argumentDimension();
  }

}
