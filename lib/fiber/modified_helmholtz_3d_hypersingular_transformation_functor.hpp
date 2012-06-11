#ifndef fiber_modified_helmholtz_3d_hypersingular_transformation_functor_hpp
#define fiber_modified_helmholtz_3d_hypersingular_transformation_functor_hpp

#include "basis_transformation_functor_wrappers.hpp"
#include "scalar_function_value_functor.hpp"
#include "surface_curl_3d_functor.hpp"

namespace Fiber
{

/** \brief Functor calculating basis-function transformations necessary for
 *  the implementation of the hypersingular operator for the modified Helmholtz
 *  equation in 3D.
 *
 *  Two following quantities are calculated:
 *  * shape-function value
 *  * shape-function surface curl. */
template <typename CoordinateType_>
class ModifiedHelmholtz3dHypersingularTransformationFunctor :
        public ElementaryBasisTransformationFunctorPairWrapper<
        ScalarFunctionValueElementaryFunctor<CoordinateType_>,
        SurfaceCurl3dElementaryFunctor<CoordinateType_> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
