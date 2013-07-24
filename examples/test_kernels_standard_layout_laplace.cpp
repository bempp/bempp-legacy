#include <stdexcept>

#include "fiber/collection_of_4d_arrays.hpp"
#include "fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "fiber/laplace_3d_single_layer_potential_collection_of_kernels.hpp"
#include "common/auto_timer.hpp"

using namespace Fiber;

int main()
{
    typedef Laplace3dSingleLayerPotentialCollectionOfKernels<double> 
        KernelCollection;  
    KernelCollection col;

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
    
    CollectionOf4dArrays<double> result(1);
    result[0].set_size(1, 1, testPointCount, trialPointCount);

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

