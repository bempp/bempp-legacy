#include "../common/common.hpp"

#include "../assembly/weak_form_hmat_assembly_helper.hpp"

#include <iostream>
#include <chrono>

namespace Bempp {

void test_cuda() {

    std::cout << "Hello! This is the CUDA auto-tuner to optimize the CUDA block size and the threshold size for computing H-Matrix blocks on the GPU." << std::endl;
    std::cout << "Info: The tuner is based on single-precision CUDA computations of Laplace SLP weak-form matrix coefficients." << std::endl;

    unsigned int cudaBlockSizeOpt = 0;
    unsigned int thresholdSizeOpt = 10000;

    typedef std::chrono::steady_clock clock;
    typedef double                    BasisFunctionType;
    typedef double                    KernelType;
    typedef double                    ResultType;
    typedef float                     CudaBasisFunctionType;
    typedef float                     CudaKernelType;
    typedef float                     CudaResultType;

//    const Space<BasisFunctionType>                                 testSpace;
//    const Space<BasisFunctionType>                                 trialSpace;
//    const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> kernel;
//    const Fiber::Shapeset<BasisFunctionType>                       testShapeset;
//    const Fiber::Shapeset<BasisFunctionType>                       trialShapeset;
//    const shared_ptr<hmat::DefaultBlockClusterTreeType>            blockClusterTree;
//    const std::vector<LocalAssembler*>                             assemblers;
//    const std::vector<const DiscreteLinOp*>                        sparseTermsToAdd;
//    const std::vector<ResultType>                                  denseTermsMultipliers;
//    const std::vector<ResultType>                                  sparseTermsMultipliers;
//    const Context<BasisFunctionType, ResultType>                   context;
//
//    typedef WeakFormHMatAssemblyHelper<
//            BasisFunctionType,     KernelType,     ResultType,
//        CudaBasisFunctionType, CudaKernelType, CudaResultType>
//    AssemblyHelper;
//
//    AssemblyHelper helper(testSpace, trialSpace, kernel, testShapeset,
//                          trialShapeset, blockClusterTree, assemblers,
//                          sparseTermsToAdd, denseTermsMultipliers,
//                          sparseTermsMultipliers, context);
//
//    // Loop over CUDA block sizes between 32 and 512
//    for (int i = 32; i <= 512; i += 32) {
//
//      helper.setCudaMinBlockSize(i);
//
//      // Loop over increasing H-Matrix block sizes
//      for (int N = 1; N < 10000; ++N) {
//
//        // CPU H-Matrix block computation
//        helper.setIsCudaEnabled(false);
//        double cpuTimer = 0.0;
//        clock::time_point cpuBegin = clock::now();
//        helper.computeMatrixBlock();
//        clock::time_point cpuEnd = clock::now();
//        cpuTimer = std::chrono::duration_cast<std::chrono::microseconds>(cpuEnd-cpuBegin).count();
//
//        // GPU H-Matrix block computation
//        helper.setIsCudaEnabled(true);
//        double gpuTimer = 0.0;
//        clock::time_point gpuBegin = clock::now();
//        helper.computeMatrixBlock();
//        clock::time_point gpuEnd = clock::now();
//        gpuTimer = std::chrono::duration_cast<std::chrono::microseconds>(gpuEnd-gpuBegin).count();
//
//        if (gpuTimer < cpuTimer) {
//          printf("CUDA block size = %d: optimal threshold size = %d.\n", i, N);
//          if (N < thresholdSizeOpt) {
//            thresholdSizeOpt = N;
//            cudaBlockSizeOpt = i;
//          }
//          break;
//        }
//      }
//    }

    std::cout << "Result: auto-tuner recommends to use a CUDA block size of "
              << cudaBlockSizeOpt << " along with a threshold size of "
              << thresholdSizeOpt << " to achieve best performance on the current machine." << std::endl;
}

}
