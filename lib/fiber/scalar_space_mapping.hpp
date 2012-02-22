#ifndef fiber_scalar_space_mapping_hpp
#define fiber_scalar_space_mapping_hpp

#include <stdexcept>
#include "basis_data.hpp"
#include "geometrical_data.hpp"

namespace Fiber
{

/** \brief Mapping of functions defined on reference elements to ones defined
  on physical elements. */
template <typename ValueType>
class ScalarSpaceMapping
{
public:
    static void addShapeFunctionDependencies(int& basisDeps, int& geomDeps) {
        basisDeps |= VALUES;
    }

    static void addSurfaceCurlDependencies(int& basisDeps, int& geomDeps) {
        throw std::logic_error("ScalarSpaceMapping::"
                               "addSurfaceCurlDependencies(): "
                               "Not implemented yet");
    }

    static void evaluateShapeFunctions(const BasisData<ValueType>& basisData,
                                       const GeometricalData<ValueType>& geomData,
                                       arma::Cube<ValueType>& result) {
        result = basisData.values;
    }

    static void evaluateSurfaceCurls(const BasisData<ValueType>& basisData,
                                     const GeometricalData<ValueType>& geomData,
                                     arma::Mat<ValueType>& result) {
        throw std::logic_error("ScalarSpaceMapping::"
                               "evaluateSurfaceCurls(): "
                               "Not implemented yet");
    }
};

} // namespace Fiber

#endif // SCALAR_SPACE_MAPPING_HPP
