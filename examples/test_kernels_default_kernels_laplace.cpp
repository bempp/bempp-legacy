#include <stdexcept>

#include "fiber/default_collection_of_kernels_imp.hpp"
#include "fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "common/auto_timer.hpp"

using namespace Fiber;

typedef double RT;

int main()
{
    typedef Laplace3dSingleLayerPotentialKernelFunctor<RT> Functor;
    typedef DefaultCollectionOfKernels<Functor> KernelCollection;  
    KernelCollection col((Functor()));

    const int pointCount = 4, dimWorld = 3;
    const int testPointCount = pointCount, trialPointCount = pointCount;
    GeometricalData<double> testGeomData, trialGeomData;
    testGeomData.globals.set_size(dimWorld, testPointCount);
    int i = 1;
    for (int d = 0; d < dimWorld; ++d)
        for (int p = 0; p < testPointCount; ++p)
            testGeomData.globals(d, p) = (double)(i++);
    trialGeomData.globals.set_size(dimWorld, trialPointCount);
    i = 1;
    for (int d = 0; d < dimWorld; ++d)
        for (int p = 0; p < trialPointCount; ++p)
            trialGeomData.globals(d, p) = -(double)(i++);
    
    CollectionOf4dArrays<RT > result(1);
    result[0].set_size(0, 0, testPointCount, trialPointCount);

    {
    Bempp::AutoTimer timer("Kernel evaluation took ");
    for (size_t i = 0; i < 10000000; ++i)
        col.evaluateOnGrid(testGeomData, trialGeomData, result);
    }
    for (size_t i = 0; i < testPointCount; ++i) {
        for (size_t j = 0; j < trialPointCount; ++j)
            std::cout << result[0](0, 0, i, j) << " ";
        std::cout << std::endl;
    }
}

