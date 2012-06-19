
#include "test_functor.hpp"
#include <armadillo>

namespace Bempp
{

  double testFunctor(const SurfaceNormalIndependentFunctor<double>& functor){
      arma::Col<double> point(3);
      point(0)=1.0; point(1)=2.0; point(2)=3.0;
      arma::Col<double> result(functor.resultDimension());
      functor.evaluate(point,result);
      for (size_t i=0;i<functor.resultDimension();i++) std::cout <<  result(i) << " ";
      std::cout << std::endl;
      return result(0);

  }

}
