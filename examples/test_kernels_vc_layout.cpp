#include <stdexcept>

#include "common/auto_timer.hpp"
#include "vc_multicomponent_kernel_values.hpp"
#include "vc_geometrical_data.hpp"

#include "vc_modified_helmholtz_3d_single_layer_potential_collection_of_kernels.hpp"

using namespace Fiber;

// typedef std::complex<double> RT;
typedef double CT;
typedef double RT;

int main()
{
    typedef VcModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels<RT> KernelCollection;  
    KernelCollection col(1., 0. /*RT(0., -1.)*/);

    const int pointCount = 4, dimWorld = 3;
    const int testPointCount = pointCount, trialPointCount = pointCount;
    const int paddedTestPointCount = ((testPointCount + Vc::Vector<CT>::Size - 1) / 
                                      Vc::Vector<CT>::Size) * Vc::Vector<CT>::Size;
    const int paddedTrialPointCount = ((testPointCount + Vc::Vector<CT>::Size - 1) / 
                                       Vc::Vector<CT>::Size) * Vc::Vector<CT>::Size;
    VcGeometricalData<double> testGeomData(testPointCount, dimWorld),
        trialGeomData(trialPointCount, dimWorld);
    int i = 1;
    for (int d = 0; d < dimWorld; ++d)
        for (int p = 0; p < testPointCount; ++p)
            testGeomData.globals[p + d * paddedTestPointCount] = (double)(i++);
    i = 1;
    for (int d = 0; d < dimWorld; ++d)
        for (int p = 0; p < trialPointCount; ++p)
            trialGeomData.globals[p + d * paddedTrialPointCount] = -(double)(i++);
    
    Fiber::CollectionOf4dArrays<Vc::Vector<CT> > result(1);
    result[0].set_size(testPointCount / Vc::Vector<CT>::Size, trialPointCount, 1, 1);
    Fiber::CollectionOf4dArrays<Vc::Vector<CT> > resultImag(1);
    resultImag[0].set_size(testPointCount / Vc::Vector<CT>::Size, trialPointCount, 1, 1);

    {
    Bempp::AutoTimer timer("Kernel evaluation took ");
    for (size_t i = 0; i < 5000000; ++i)
        col.evaluateOnGrid(testGeomData, trialGeomData, result, resultImag);
    }

    const int testPointChunkCount = testPointCount / Vc::Vector<CT>::Size;
    for (size_t i = 0; i < testPointChunkCount; ++i)
        for (size_t k = 0; k < Vc::Vector<CT>::Size && i * Vc::Vector<CT>::Size + k < testPointCount; ++k) {
            for (size_t j = 0; j < trialPointCount; ++j)
                std::cout << "(" << result[0](i, j, 0, 0)[k] << ", " 
                          << resultImag[0](i, j, 0, 0)[k] << 
                    ") ";
        std::cout << std::endl;
    }

}

