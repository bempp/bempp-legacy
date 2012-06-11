#ifndef fiber_scalar_function_value_functor_hpp
#define fiber_scalar_function_value_functor_hpp

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"
#include "basis_transformation_functor_wrappers.hpp"

namespace Fiber
{

template <typename CoordinateType_>
class ScalarFunctionValueElementaryFunctor
{
public:
    typedef CoordinateType_ CoordinateType;

    int argumentDimension() const { return 1; }
    int resultDimension() const { return 1; }

    void addDependencies(int& basisDeps, int& geomDeps) const {
        basisDeps |= VALUES;
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            _1dSliceOf3dArray<ValueType>& result) const {
        assert(basisData.componentCount() == 1);
        result(0) = basisData.values(0);
    }
};

// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
template <typename CoordinateType_>
class ScalarFunctionValueFunctor :
        public ElementaryBasisTransformationFunctorWrapper<
        ScalarFunctionValueElementaryFunctor<CoordinateType_> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
