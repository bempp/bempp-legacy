#ifndef fiber_modified_helmholtz_3d_hypersingular_integrand_functor_hpp
#define fiber_modified_helmholtz_3d_hypersingular_integrand_functor_hpp

#include <cassert>
#include "collection_of_3d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class ModifiedHelmholtz3dHypersingularIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    ModifiedHelmholtz3dHypersingularIntegrandFunctor(KernelType waveNumber) :
        m_waveNumber(waveNumber)
    {}

    void addGeometricalDependencies(int& testGeomDeps, int& trialGeomDeps) const {
        testGeomDeps |= NORMALS;
        trialGeomDeps |= NORMALS;
    }

    KernelType waveNumber() const { return m_waveNumber; }
    void setWaveNumber(KernelType k) { m_waveNumber = k; }

    // It is possible that this function could be generalised to
    // multiple basis transformations or kernels and that the additional
    // loops could be optimised away by the compiler.
    template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& testTransfValues,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& trialTransfValues,
            const CollectionOf2dSlicesOfConstNdArrays<KernelType>& kernelValues) const {
        const int dimWorld = 3;

        // Assert that there is a single scalar-valued kernel
        assert(kernelValues.size() == 1);
        assert(kernelValues[0].extent(0) == 1);
        assert(kernelValues[0].extent(1) == 1);

        // Assert that there are two test and trial transformations
        // (function value and surface curl) of correct dimensions
        assert(testTransfValues.size() == 2);
        assert(trialTransfValues.size() == 2);
        _1dSliceOfConst3dArray<BasisFunctionType> testValues = testTransfValues[0];
        _1dSliceOfConst3dArray<BasisFunctionType> trialValues = trialTransfValues[0];
        _1dSliceOfConst3dArray<BasisFunctionType> testSurfaceCurls = testTransfValues[1];
        _1dSliceOfConst3dArray<BasisFunctionType> trialSurfaceCurls = trialTransfValues[1];
        assert(testValues.extent(0) == 1);
        assert(trialValues.extent(0) == 1);
        assert(testSurfaceCurls.extent(0) == dimWorld);
        assert(trialSurfaceCurls.extent(0) == dimWorld);

        // K(x, y) [kappa^2 u*(x) v(y) n(x) . n(y) + curl u*(x) . curl v(y)]

        ResultType result = 0.;
        for (int dim = 0; dim < dimWorld; ++dim)
            result += testGeomData.normal(dim) * trialGeomData.normal(dim);
        result *= m_waveNumber * m_waveNumber *
                conjugate(testValues(0)) * trialValues(0);
        for (int dim = 0; dim < dimWorld; ++dim)
            result += conjugate(testSurfaceCurls(dim)) *
                    trialSurfaceCurls(dim);
        result *= kernelValues[0](0, 0);
        return result;
    }

private:
    KernelType m_waveNumber;
};

} // namespace Fiber

#endif
