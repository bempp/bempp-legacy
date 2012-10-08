#ifdef WITH_AHMED

%{
#include "assembly/aca_approximate_lu_inverse.hpp"
#include "assembly/discrete_aca_boundary_operator.hpp"
#include "fiber/scalar_traits.hpp"
#include "common/shared_ptr.hpp"
%}

%define NUMERICAL_ACA_LU_APPROXIMATE_INVERSE(VALUE,PY_VALUE)
     %inline %{
     namespace Bempp{
     boost::shared_ptr<const DiscreteBoundaryOperator< VALUE > > createAcaApproximateLuInverse_## PY_VALUE (
             const boost::shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op,
             Fiber::ScalarTraits< VALUE >::RealType delta){

         const DiscreteAcaBoundaryOperator< VALUE >& acaOp =
                 DiscreteAcaBoundaryOperator< VALUE >::castToAca(*op);
         boost::shared_ptr<const DiscreteBoundaryOperator< VALUE > > acaLu(new AcaApproximateLuInverse< VALUE >(acaOp,delta));
         return acaLu;

     }
 }
     %}
%enddef
BEMPP_ITERATE_OVER_BASIS_TYPES(NUMERICAL_ACA_LU_APPROXIMATE_INVERSE)

#endif // WITH_AHMED
