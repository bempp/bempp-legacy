#ifndef fiber_surface_curl_3d_functor_hpp
#define fiber_surface_curl_3d_functor_hpp

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"
#include "basis_transformation_functor_wrappers.hpp"

#include <boost/array.hpp>

namespace Fiber
{

template <typename CoordinateType_>
class SurfaceCurl3dElementaryFunctor
{
public:
    typedef CoordinateType_ CoordinateType;

    int argumentDimension() const { return 1.; }
    int resultDimension() const { return 3; }

    void addDependencies(int& basisDeps, int& geomDeps) const {
        basisDeps |= DERIVATIVES;
        geomDeps |= NORMALS | JACOBIAN_INVERSES_TRANSPOSED;
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            _1dSliceOf3dArray<ValueType>& result) const {
        assert(basisData.componentCount() == 1);
        assert(geomData.dimWorld() == 3);

        // vec := gradient of the basis function extended outside
        // the surface so that its normal derivative on the surf. is zero
        boost::array<ValueType, 3> vec;
        vec[0] = basisData.derivatives(0, 0) * geomData.jacobianInverseTransposed(0, 0) +
                basisData.derivatives(0, 1) * geomData.jacobianInverseTransposed(0, 1);
        vec[1] = basisData.derivatives(0, 0) * geomData.jacobianInverseTransposed(1, 0) +
                basisData.derivatives(0, 1) * geomData.jacobianInverseTransposed(1, 1);
        vec[2] = basisData.derivatives(0, 0) * geomData.jacobianInverseTransposed(2, 0) +
                basisData.derivatives(0, 1) * geomData.jacobianInverseTransposed(2, 1);

        // result := n \times vec
        result(0) = geomData.normal(1) * vec[2] - geomData.normal(2) * vec[1];
        result(1) = geomData.normal(2) * vec[0] - geomData.normal(0) * vec[2];
        result(2) = geomData.normal(0) * vec[1] - geomData.normal(1) * vec[0];
    }
};

// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
template <typename CoordinateType_>
class SurfaceCurl3dFunctor :
        public ElementaryBasisTransformationFunctorWrapper<
        SurfaceCurl3dElementaryFunctor<CoordinateType_> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
