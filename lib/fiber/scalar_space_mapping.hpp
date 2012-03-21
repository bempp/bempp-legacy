#ifndef fiber_scalar_space_mapping_hpp
#define fiber_scalar_space_mapping_hpp

#include <cassert>
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
        basisDeps |= DERIVATIVES;
        geomDeps |= NORMALS | JACOBIAN_INVERSES_TRANSPOSED;
    }

    static void evaluateShapeFunctions(const BasisData<ValueType>& basisData,
                                       const GeometricalData<ValueType>& geomData,
                                       arma::Cube<ValueType>& result) {
        result = basisData.values;
    }

    /**
      \param[out] result
        dimensions: (worldDim, functionCount, pointCount)
    */
    static void evaluateSurfaceCurls3D(const BasisData<ValueType>& basisData,
                                     const GeometricalData<ValueType>& geomData,
                                     arma::Cube<ValueType>& result) {
        const arma::Mat<ValueType>& n = geomData.normals;
        // jt(i, j): dx_j/dq_i
        const arma::Cube<ValueType>& jit = geomData.jacobianInversesTransposed;
        const Array4D<ValueType>& d = basisData.derivatives;

        assert(d.extent(0) == 1); // scalar functions
        const int worldDim = 3;
        assert(d.extent(1) + 1 == worldDim);
        const int functionCount = d.extent(2);
        const int pointCount = d.extent(3);

        result.set_size(worldDim, functionCount, pointCount);

        for (int p = 0; p < pointCount; ++p)
            for (int f = 0; f < functionCount; ++f) {
                // vec := gradient of the basis function extended outside
                // the surface so that its normal derivative on the surf. is zero
                arma::Col<ValueType> vec(3);
                vec(0) = d(0, 0, f, p) * jit(0, 0, p) + d(0, 1, f, p) * jit(0, 1, p);
                vec(1) = d(0, 0, f, p) * jit(1, 0, p) + d(0, 1, f, p) * jit(1, 1, p);
                vec(2) = d(0, 0, f, p) * jit(2, 0, p) + d(0, 1, f, p) * jit(2, 1, p);

                // result := n \times vec
                result(0, f, p) = n(1, p) * vec(2) - n(2, p) * vec(1);
                result(1, f, p) = n(2, p) * vec(0) - n(0, p) * vec(2);
                result(2, f, p) = n(0, p) * vec(1) - n(1, p) * vec(0);
            }
    }

};

} // namespace Fiber

#endif // SCALAR_SPACE_MAPPING_HPP
