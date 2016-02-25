// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "global_parameters.hpp"

namespace Bempp {

ParameterList GlobalParameters::parameterList() {

  ParameterList parameters;

  // Specifies the maximum number of threads to use
  // -1 denotes automatic selection of number of threads via TBB.
  parameters.put("options.global.maxThreadCount", static_cast<int>(-1));

  // Default Verbosity of BEM++. Supported values are
  // -5 (low verbosity), 0 (default), 5 (high verbosity)
  parameters.put("options.global.verbosityLevel", static_cast<int>(5));

  // Default assembly type for boundary operators. Allowed values are
  // "dense" and "hmat".
  parameters.put("options.assembly.boundaryOperatorAssemblyType",
                 std::string("hmat"));

  // Default assembly type for potential oeprators.
  // Allowed values are "dense" and "hmat".
  parameters.put("options.assembly.potentialOperatorAssemblyType",
                 std::string("hmat"));

  // If true then singular integrals are pre-calculated and cached
  // before the boundary oeprator assembly.

  parameters.put("options.assembly.enableSingularIntegralCaching", true);

  // Use polynomial interpolation instead of exponentials to assemble
  // Helmholtz or Maxwell type kernels.
  parameters.put("options.assembly.enableInterpolationForOscillatoryKernels", true);

  // Number of interpolation points per wavelength for oscillatory kernels.
  parameters.put("options.assembly.interpolationPointsPerWavelength", static_cast<int>(5000));
   

  // Order for singular double integrals.
  parameters.put("options.quadrature.doubleSingular", static_cast<int>(6));

  auto createQuadratureOptions = [&parameters](const std::string name,
                                               double relDist, int singleOrder,
                                               int doubleOrder) {

    // Relative distance of quadrature point to element.
    parameters.put(
        (std::string("options.quadrature.") + name + std::string(".maxRelDist"))
            .c_str(),
        static_cast<double>(relDist));

    // Order of single regular integrals.
    parameters.put((std::string("options.quadrature.") + name +
                    std::string(".singleOrder"))
                       .c_str(),
                   static_cast<int>(singleOrder));

    // Order of double regular integrals.
    parameters.put((std::string("options.quadrature.") + name +
                    std::string(".doubleOrder"))
                       .c_str(),
                   static_cast<int>(doubleOrder));

  };

  createQuadratureOptions("near", 2, 4, 4);
  createQuadratureOptions("medium", 4, 3, 3);
  createQuadratureOptions("far", std::numeric_limits<double>::infinity(), 2, 2);

  parameters.erase("options.quadrature.far.maxRelDist");

  // Specifies the minimum block size below which blocks are assumed to be dense
  parameters.put("options.hmat.minBlockSize", static_cast<int>(20));

  // Specifies the maximum size of an admissible block
  parameters.put("options.hmat.maxBlockSize", static_cast<int>(1000000));

  // Specifies the block separation parameter eta.
  parameters.put("options.hmat.eta", static_cast<double>(1.2));

  // Specifies the type of admissibility function ('strong' or 'weak')
  parameters.put("options.hmat.admissibility", std::string("weak"));

  // Specifies the accuracy of low-rank approximations.
  parameters.put("options.hmat.eps", static_cast<double>(1E-3));

  // Maximum rank of a low rank subblock
  parameters.put("options.hmat.maxRank", static_cast<int>(30));

  // Compression algorithm
  parameters.put("options.hmat.compressionAlgorithm", std::string("aca"));

  // Enable coarsening
  parameters.put("options.hmat.coarsening", true);

  // Accuracy for coarsening
  // 0: Use same as options.hmat.eps
  parameters.put("options.hmat.coarseningAccuracy", static_cast<double>(0));

  // Number of levels for matvec parallelisation
  // The total number of tasks is 4^matVecParallelLevels
  parameters.put("options.hmat.matVecParallelLevels", static_cast<int>(5));

  return parameters;
}
}
